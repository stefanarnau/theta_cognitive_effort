function[EEG] = event_coding(EEG, RESPS, positions, trial_log, rt_thresh_color, rt_thresh_tilt)

    % Struct for new events
    nec = 0;
    new_events = struct();

    % Prepare response channels
    idx_left = find(strcmp({RESPS.chanlocs.labels}, 'LeftKey'));
    idx_right = find(strcmp({RESPS.chanlocs.labels}, 'RightKey'));
    resps_left = rescale(RESPS.data(idx_left, :));
    resps_right = rescale(RESPS.data(idx_right, :));
    critval_responses = 0.5;
    maxrt = 1500;

    % Look up correct color responses
    if ismember(EEG.task_version, [1, 2, 5, 6])
        correct_green = 0;
        correct_red = 1;
    else
        correct_green = 1;
        correct_red = 0;
    end

    % Iterate events recoding rt and accuracy info
    trial_nr = 0;
    for e = 1 : length(EEG.event)

        % Check for block start markers
        if (strcmpi(EEG.event(e).type(1), {'S'}) & ismember(str2num(EEG.event(e).type(2 : 4)), [120 : 130]))

            % Get block number
            block_nr = str2num(EEG.event(e).type(2 : 4)) - 120;

            % Rest sequential marker
            previous_task = -1;

        end

        % If trial
        if (strcmpi(EEG.event(e).type(1), {'S'}) & ismember(str2num(EEG.event(e).type(2 : 4)) - 40, [0 : 32]))
    
            trial_nr = trial_nr + 1;

            % Get event code
            ecode = str2num(EEG.event(e).type(2 : 4)) - 40;

            % Decode bonustrial
            if ecode <= 16
                bonustrial = 1;
            else
                bonustrial = 0;
            end

            % Decode task
            if ismember(ecode, [1 : 8, 17 : 24])
                tilt_task = 1;
            else
                tilt_task = 0;
            end

            % Decode cue
            if ismember(mod(ecode, 8), [1, 2, 3, 4]) 
                cue_ax = 1;
            else
                cue_ax = 0;
            end

            % Decode target
            if ismember(mod(ecode, 4), [1, 2]) 
                target_red_left = 1;
            else
                target_red_left = 0;
            end

            % Decode distractor
            if mod(ecode, 2) == 1
                distractor_red_left = 1;
            else
                distractor_red_left = 0;
            end

            % Check response interference
            if target_red_left == distractor_red_left
                response_interference = 0;
            else
                response_interference = 1;
            end

            % Code task sequence
            current_task = tilt_task;
            if previous_task == -1
                task_switch = -1;
            elseif current_task == previous_task
                task_switch = 0;
            else
                task_switch = 1;
            end
            previous_task = current_task;

            % Decode correct response side

            % If tilt task target left
            if tilt_task == 1 & target_red_left == 1
                correct_response = 0;
            % If tilt task target right
            elseif tilt_task == 1 & target_red_left == 0
                correct_response = 1;
            % If color task target red 
            elseif tilt_task == 0 & target_red_left == 1
                correct_response = correct_red;
            % If color task target green
            elseif tilt_task == 0 & target_red_left == 0
                correct_response = correct_green;
            end

            % Create event
            nec = nec + 1;
            new_events(nec).latency = EEG.event(e).latency;
            new_events(nec).duration = 1;
            new_events(nec).type = "trial";
            new_events(nec).code = "trial";
            new_events(nec).urevent = EEG.event(e).urevent;
            new_events(nec).block_nr = block_nr;
            new_events(nec).trial_nr = trial_nr;
            new_events(nec).bonustrial = bonustrial;
            new_events(nec).tilt_task = tilt_task;
            new_events(nec).cue_ax = cue_ax;
            new_events(nec).target_red_left = target_red_left;
            new_events(nec).distractor_red_left = distractor_red_left;
            new_events(nec).response_interference = response_interference;
            new_events(nec).task_switch = task_switch;
            new_events(nec).correct_response = correct_response;
            new_events(nec).position_color = positions(trial_nr, 1);
            new_events(nec).position_tilt = positions(trial_nr, 2);
            new_events(nec).rt_thresh_color = rt_thresh_color;
            new_events(nec).rt_thresh_tilt = rt_thresh_tilt;

            % Code log response side
            if trial_log(trial_nr, 1) == 0
                response_side = 2;
            elseif trial_log(trial_nr, 1) == 1
                response_side = 0;
            elseif trial_log(trial_nr, 1) == 2
                response_side = 1;
            end
            new_events(nec).response_side = response_side;

            % Code log rt
            rt = trial_log(trial_nr, 3);
            if rt == -1
                rt = NaN;
            end
            if rt > maxrt
                rt = NaN;
            end
            new_events(nec).rt = rt;

            % Code accuracy
            accuracy = trial_log(trial_nr, 2);
            if accuracy == 3
                accuracy = 2;
            end
            if isnan(rt)
                accuracy = 2;
            end
            new_events(nec).accuracy = accuracy;
            
            if tilt_task == 0
                new_events(nec).position_target = new_events(nec).position_color;
                new_events(nec).position_distractor = new_events(nec).position_tilt;
            else
                new_events(nec).position_target = new_events(nec).position_tilt;
                new_events(nec).position_distractor = new_events(nec).position_color; 
            end
            
            % Mark bugged trials
            if ecode == 0
                new_events(nec).type = "zerocoded";
                new_events(nec).code = "zerocoded";
            end
        
        end

        % If event is boundary event...
        if strcmpi(EEG.event(e).type, 'boundary')
            stimulus_type = 'boundary';
        end

    end

    % Get indices of trials
    trial_idx = [];
    for e = 1 : length(new_events)
        if strcmpi(new_events(e).type, "trial")
            trial_idx(end + 1) = e;
        end
    end

    % Code sequences
    old_block = 0;
    n_nontrial = 0;
    for tidx = 1 : length(trial_idx)

        % Get trial index in event structure
        e = trial_idx(tidx);

        % If block changes
        current_block = new_events(e).block_nr;
        if current_block ~= old_block

            old_block = current_block;
            sequence_position = 1;
            new_events(e).task_switch = -1;
            sequence_start = e;
            previous_type = new_events(e).bonustrial;
            current_type = new_events(e).bonustrial;

            switch_prev = -1;
            switch_before = new_events(e).task_switch;

            acc_prev = -1;
            acc_before = new_events(e).accuracy;

        % If block does not change
        else

            % If sequence continues
            current_type = new_events(e).bonustrial;

            if current_type == previous_type

                sequence_position = sequence_position + 1;
                previous_type = new_events(e).bonustrial;

                switch_prev = switch_before;
                switch_before = new_events(e).task_switch;

                acc_prev = acc_before;
                acc_before = new_events(e).accuracy;

            % If sequence does not continue
            else

                % Reset sequence parameters
                sequence_position = 1;
                new_events(e).task_switch = -1;
                sequence_start = e;
                previous_type = new_events(e).bonustrial;

                switch_prev = -1;
                switch_before = new_events(e).task_switch;

                acc_prev = -1;
                acc_before = new_events(e).accuracy;
            end
            
        end

        new_events(e).sequence_position = sequence_position;
        new_events(e).prev_accuracy = acc_prev;
        new_events(e).prev_switch = switch_prev;

    end

    % Check if number of trials match
    if trial_nr ~= size(positions, 1)
        error('number of trials in logfile and number of trials in event markers do not match...')
    end

    % Replace events
    EEG.event = new_events;
    EEG = eeg_checkset(EEG, 'eventconsistency');

end
