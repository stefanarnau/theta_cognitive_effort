clear all;

% PATH VARS
PATH_EEGLAB        = 'insert_path_here';
PATH_LOGFILES      = 'insert_path_here';
PATH_RAW           = 'insert_path_here';
PATH_ICSET         = 'insert_path_here';
PATH_AUTOCLEANED   = 'insert_path_here';

% Subjects
subject_list = {'VP09', 'VP17', 'VP25', 'VP10', 'VP11', 'VP13', 'VP14', 'VP15', 'VP16', 'VP18',...
                'VP19', 'VP20', 'VP21', 'VP22', 'VP23', 'VP08', 'VP24', 'VP26', 'VP27', 'VP28',...
                'VP29', 'VP30', 'VP31', 'VP32', 'VP33', 'VP34'};

% Init eeglab
addpath(PATH_EEGLAB);
eeglab;
channel_location_file = which('dipplot.m');
channel_location_file = channel_location_file(1 : end - length('dipplot.m'));
channel_location_file = [channel_location_file, 'standard_BESA/standard-10-5-cap385.elp'];

% Iterate subjects
for s = 1 : length(subject_list)

    % participant identifiers
    subject = subject_list{s};
    id = str2num(subject(3 : 4));

    % Load
    EEG = pop_loadbv(PATH_RAW, [subject, '.vhdr'], [], []);

    % Repair subject 13 (first block start marker missing)
    if id == 13
        EEG = pop_editeventvals(EEG, 'insert',...
                                {1, [], [], [], [], [], [], [], [], []},...
                                'changefield', {1, 'latency', 0.5},...
                                'changefield', {1, 'duration', 0.001},...
                                'changefield', {1, 'channel', 0},...
                                'changefield', {1, 'bvtime', []},...
                                'changefield', {1, 'visible', []},...
                                'changefield', {1, 'bvmknum', 3733},...
                                'changefield', {1, 'type', 'S121'},...
                                'changefield', {1, 'code', 'Stimulus'});
    end

    % Subject 30 has been restarted after a couple of trials.
    % Remove all events until second start of block 1 (second occurence of 'S121'),
    % which is event number 36...
    if id == 30
        EEG.event(1 : 35) = [];
    end

    % Fork response button channels
    RESPS = pop_select(EEG, 'channel', [65, 66]);
    EEG = pop_select(EEG, 'nochannel', [65, 66]);

    % Open log file
    fid = fopen([PATH_LOGFILES, subject, '_degreeLog.txt'], 'r');

    % Extract lines as strings
    logcell = {};
    tline = fgetl(fid);
    while ischar(tline)
        logcell{end + 1} = tline;
        tline = fgetl(fid);
    end

    % Delete header
    logcell(1 : 3) = [];

    % Get color and tilt positions in probe display (numbers 1-8)
    positions = [];
    for l = 1 : length(logcell)
        line_values = split(logcell{l}, ' ');
        positions(l, 1) = str2num(line_values{8});
        positions(l, 2) = str2num(line_values{10});
    end

    % Open trial log file
    fid = fopen([PATH_LOGFILES, subject, '_trials.txt'], 'r');

    % Extract lines as strings
    logcell = {};
    tline = fgetl(fid);
    while ischar(tline)
        logcell{end + 1} = tline;
        tline = fgetl(fid);
    end

    % Delete header
    logcell(1 : 3) = [];

    % Get response side, accuracy and rt from log file
    trial_log = [];
    for l = 1 : length(logcell)
        line_values = split(logcell{l}, '|');
        trial_log(l, 1) = str2num(line_values{5});
        trial_log(l, 2) = str2num(line_values{6});
        trial_log(l, 3) = str2num(line_values{7});
    end

    % Get version of task
    if id < 8
        error("Preprocessing invalid for id < 8.");
    elseif id == 8
        EEG.task_version = 1;
    else
        EEG.task_version = mod(id, 8);
        if EEG.task_version == 0
            EEG.task_version = 8;
        end
    end

     % Open log file
     fid = fopen([PATH_LOGFILES, subject, '_degreeLog.txt'], 'r');

     % Extract lines as strings
     logcell = {};
     tline = fgetl(fid);
     while ischar(tline)
         logcell{end + 1} = tline;
         tline = fgetl(fid);
     end

     % Delete header
     logcell(1 : 3) = [];

     % Iterate last 100 trials and extract rt thresholds
     rt_threshs = [];
     for l = 1 : 100
         line_values = split(logcell{length(logcell) - l}, ' ');
         rt_threshs(l, 1) = str2num(line_values{5});
         rt_threshs(l, 2) = str2num(line_values{13});
     end
     rt_thresh_color = mean(rt_threshs(rt_threshs(:, 1) == 2, 2));
     rt_thresh_tilt = mean(rt_threshs(rt_threshs(:, 1) == 1, 2));

    % Event coding
    EEG = event_coding(EEG, RESPS, positions, trial_log, rt_thresh_color, rt_thresh_tilt);

    % Add FCz as empty channel
    EEG.data(end + 1, :) = 0;
    EEG.nbchan = size(EEG.data, 1);
    EEG.chanlocs(end + 1).labels = 'FCz';

    % Add channel locations
    EEG = pop_chanedit(EEG, 'lookup', channel_location_file);

    % Save original channel locations (for later interpolation)
    EEG.chanlocs_original = EEG.chanlocs;

    % Reref to CPz, so that FCz obtains non-interpolated data
    EEG = pop_reref(EEG, 'CPz');

    % Remove data at boundaries
    EEG = pop_rmdat(EEG, {'boundary'}, [0, 1], 1);

    % Resample data
    ERP = pop_resample(EEG, 500);
    EEG = pop_resample(EEG, 200);

    % Reject continuous data
    [ERP, selected_regions] = pop_rejcont(ERP, 'freqlimit', [20, 40], 'taper', 'hamming');
    ERP.rejcont_regions = selected_regions;
    [EEG, selected_regions] = pop_rejcont(EEG, 'freqlimit', [20, 40], 'taper', 'hamming');
    EEG.rejcont_regions = selected_regions;

    % Filter
    ERP = pop_basicfilter(ERP, [1 : ERP.nbchan], 'Cutoff', [0.01, 40], 'Design', 'butter', 'Filter', 'bandpass', 'Order', 4, 'RemoveDC', 'on', 'Boundary', 'boundary'); 
    EEG = pop_basicfilter(EEG, [1 : EEG.nbchan], 'Cutoff', [   1, 40], 'Design', 'butter', 'Filter', 'bandpass', 'Order', 4, 'RemoveDC', 'on', 'Boundary', 'boundary');
        
    % Bad channel detection
    [ERP, ERP.chans_rejected] = pop_rejchan(ERP, 'elec', [1 : ERP.nbchan], 'threshold', 5, 'norm', 'on', 'measure', 'kurt');
    [EEG, EEG.chans_rejected] = pop_rejchan(EEG, 'elec', [1 : EEG.nbchan], 'threshold', 5, 'norm', 'on', 'measure', 'kurt');

    % Interpolate channels
    ERP = pop_interp(ERP, ERP.chanlocs_original, 'spherical');
    EEG = pop_interp(EEG, EEG.chanlocs_original, 'spherical');

    % Reref common average
    ERP = pop_reref(ERP, []);
    EEG = pop_reref(EEG, []);

    % Determine rank of data
    dataRank = sum(eig(cov(double(EEG.data'))) > 1e-6); 

    % Epoch data
    ERP = pop_epoch(ERP, {'trial'}, [-1, 2.6], 'newname', [subject '_epoched'], 'epochinfo', 'yes');
    ERP = pop_rmbase(ERP, [-200, 0]);
    EEG = pop_epoch(EEG, {'trial'}, [-1, 2.6], 'newname', [subject '_epoched'], 'epochinfo', 'yes');
    EEG = pop_rmbase(EEG, [-200, 0]);

    % Autoreject trials
    [ERP, ERP.rejected_epochs] = pop_autorej(ERP, 'nogui', 'on');
    [EEG, EEG.rejected_epochs] = pop_autorej(EEG, 'nogui', 'on');

    % Find standard latency of event in epoch
    lats = [];
    for e = 1 : length(ERP.event)
        lats(end+1) = mod(ERP.event(e).latency, ERP.pnts);
    end
    lat_mode = mode(lats);
    
    % Compile a trialinfo matrix
    trialinfo = [];
    counter = 0;
    for e = 1 : length(ERP.event)
        if strcmpi(ERP.event(e).type, 'trial') & (mod(ERP.event(e).latency, ERP.pnts) == lat_mode)

            counter = counter + 1;

            % Compile table
            trialinfo(counter, :) = [id,...
                                        ERP.event(e).block_nr,...
                                        ERP.event(e).trial_nr,...
                                        ERP.event(e).bonustrial,...
                                        ERP.event(e).tilt_task,...
                                        ERP.event(e).cue_ax,...
                                        ERP.event(e).target_red_left,...
                                        ERP.event(e).distractor_red_left,...
                                        ERP.event(e).response_interference,...
                                        ERP.event(e).task_switch,...
                                        ERP.event(e).prev_switch,...
                                        ERP.event(e).prev_accuracy,...
                                        ERP.event(e).correct_response,...
                                        ERP.event(e).response_side,...
                                        ERP.event(e).rt,...
                                        ERP.event(e).rt_thresh_color,...
                                        ERP.event(e).rt_thresh_tilt,...
                                        ERP.event(e).accuracy,...
                                        ERP.event(e).position_color,...
                                        ERP.event(e).position_tilt,...
                                        ERP.event(e).position_target,...
                                        ERP.event(e).position_distractor,...    
                                        ERP.event(e).sequence_position,...                  
                                        ];

        end
    end

    % Save trialinfo
    ERP.trialinfo = trialinfo;
    writematrix(trialinfo, [PATH_AUTOCLEANED, subject, '_trialinfo_erp.csv']);

    % Find standard latency of event in epoch
    lats = [];
    for e = 1 : length(EEG.event)
        lats(end+1) = mod(EEG.event(e).latency, EEG.pnts);
    end
    lat_mode = mode(lats);
    
    
    % Compile a trialinfo matrix
    trialinfo = [];
    counter = 0;
    for e = 1 : length(EEG.event)
        if strcmpi(EEG.event(e).type, 'trial') & (mod(EEG.event(e).latency, EEG.pnts) == lat_mode)

            counter = counter + 1;

            % Compile table
            trialinfo(counter, :) = [id,...
                                     EEG.event(e).block_nr,...
                                     EEG.event(e).trial_nr,...
                                     EEG.event(e).bonustrial,...
                                     EEG.event(e).tilt_task,...
                                     EEG.event(e).cue_ax,...
                                     EEG.event(e).target_red_left,...
                                     EEG.event(e).distractor_red_left,...
                                     EEG.event(e).response_interference,...
                                     EEG.event(e).task_switch,...
                                     EEG.event(e).prev_switch,...
                                     EEG.event(e).prev_accuracy,...
                                     EEG.event(e).correct_response,...
                                     EEG.event(e).response_side,...
                                     EEG.event(e).rt,...
                                     EEG.event(e).rt_thresh_color,...
                                     EEG.event(e).rt_thresh_tilt,...
                                     EEG.event(e).accuracy,...
                                     EEG.event(e).position_color,...
                                     EEG.event(e).position_tilt,...
                                     EEG.event(e).position_target,...
                                     EEG.event(e).position_distractor,...    
                                     EEG.event(e).sequence_position,...                  
                                     ];

        end
    end

    % Save trialinfo
    EEG.trialinfo = trialinfo;
    writematrix(trialinfo, [PATH_AUTOCLEANED, subject, '_trialinfo.csv']);

    % Runica & ICLabel
    EEG = pop_runica(EEG, 'extended', 1, 'interrupt', 'on', 'PCA', dataRank);
    EEG = iclabel(EEG);

    % Find nobrainer
    EEG.nobrainer = find(EEG.etc.ic_classification.ICLabel.classifications(:, 1) < 0.3 | EEG.etc.ic_classification.ICLabel.classifications(:, 3) > 0.3);

    % Copy ICs to erpset
    ERP = pop_editset(ERP, 'icachansind', 'EEG.icachansind', 'icaweights', 'EEG.icaweights', 'icasphere', 'EEG.icasphere');
    ERP.etc = EEG.etc;
    ERP.nobrainer = EEG.nobrainer;

    % Save IC set
    pop_saveset(ERP, 'filename', [subject, '_icset_erp.set'], 'filepath', PATH_ICSET, 'check', 'on');
    pop_saveset(EEG, 'filename', [subject, '_icset.set'], 'filepath', PATH_ICSET, 'check', 'on');

    % Remove components
    ERP = pop_subcomp(ERP, ERP.nobrainer, 0);
    EEG = pop_subcomp(EEG, EEG.nobrainer, 0);

    % Save clean data
    pop_saveset(ERP, 'filename', [subject, '_cleaned_erp.set'], 'filepath', PATH_AUTOCLEANED, 'check', 'on');
    pop_saveset(EEG, 'filename', [subject, '_cleaned.set'], 'filepath', PATH_AUTOCLEANED, 'check', 'on');

    % Save channel label in order for creating 10-20 montage in mne
    channel_labels = '';
    for ch = 1 : EEG.nbchan
        channel_labels = [channel_labels, ' ', EEG.chanlocs(ch).labels];
    end
    save([PATH_AUTOCLEANED, 'channel_labels.mat'], 'channel_labels');
    
end