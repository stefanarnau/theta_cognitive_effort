% Compute beamformer in the frequency domain

clear all 
close all
clc


%% Params
PlotHeadmodel=0; 
% normally is 0, for Bocotilt: if 1, then time zero of each trials's epoch is response onset. Epoch go from -2000 till response onset. StartSec=-0.4 would mean 400 ms prior to response onset
% normally bocotilt epochs have there time zero at task cue onset
ResponseLockedEpochs=0; 

% ProjectName='Auvibo';  
%     CondName='All'; Version='_V1_';    % normally it is V1, only if you use different versions within the same PRxxx number change it
%         CondLabel={'BreakAll','NoBreakAll','Bon','Stan','Au','Vi','Prcue','Imprcue','Aucue','Vicue',... % 1-10
%         'AuBon','AuStan','ViBon','ViStan',...% 11-14
%         'PrcueBon','PrcueStan','ImprcueBon','ImprcueStan','ImprcueAuBon','ImprcueAuStan',... % 15-20
%         'AucueBon','AucueStan','VicueBon','VicueStan',... % 21-24
%         'PrcueAuBon','PrcueAuStan','PrcueViBon','PrcueViStan','ImprcueViBon','ImprcueViStan',... % 27-30
%         'PrcueVi','ImprcueVi','PrcueAu','ImprcueAu',...%31-34,,all conditions until here are without break trials! This line is different from AllOld and behavioral data
%         'BreakBon','BreakStan',...
%         };
%    SubjList={'005','006','007','008','009','010','011','012','013','014','015','016','019','020','021','022','023','024','025','027','028','029','030','031','032','033','034','035','036','042','049','050'};  %'026',
ProjectName='Bocotilt';  
    CondName='LikeSensorLevel'; Version='_V1_';    % normally it is V1, only if you use different versions within the same PRxxx number change it
        CondLabel={'All','Bon','Stan','Switch','Repeat','LeftAnswer','RightAnswer'}; 
   %SubjList={'008','010','011','013','014','015','016','017','018','019','020','021','023','024','025','026','027','028','029','030','031','032','033','034'}; % '022', erst ab DI501
   SubjList={'010','011','013','014','015','016','017','018','019','020','021','022','023','024','025','026','027','028','029','030','031','032','034'}; % '008','033'(nodipoli headmodel, bad scan)
%    %SubjList={'010','011'};

if isunix, PathMain='/media/sf_nathaExtern'; else PathMain = 'I:'; end  % alternatively ispc for windows
    PathEeglab = fullfile(PathMain,'Toolboxes','eeglab2021.1'); addpath(PathEeglab); eeglab; 
    PathFT = fullfile(PathMain,'Toolboxes','fieldtrip-20210413'); addpath(PathFT); ft_defaults;
    PathLib1 = fullfile(PathMain,'Toolboxes','FunctionLib'); addpath(genpath(PathLib1));
% Pathes belonging to MRI preprocessing:
HeadmodelNum='EL001HM001BN001SG001AL002'; SourcemodelNum='SM500AL002'; 
    PathElec=fullfile(PathMain,[ProjectName 'DataBase'],HeadmodelNum);
    PathHeadmodel=fullfile(PathMain,[ProjectName 'DataBase'],HeadmodelNum(end-19:end));
    PathSourcemodel=fullfile(PathMain,[ProjectName 'DataBase'],SourcemodelNum);
% Pathes belonging to EEG preprocessing and saving path:
AnalysisName='DI500PR010'; % analysis code of outcome of this script, if ='C004P103'; then P001 is loaded
    PathPreProcMat = fullfile(PathMain,[ProjectName AnalysisName(6:end)],'OutputMats');
    PathPreProc = fullfile(PathMain,[ProjectName AnalysisName(6:end)],'3 preprocessed'); 
    PathTrialInfo = fullfile(PathMain,[ProjectName AnalysisName(6:end)], 'Trialinfo' );
    PathFilter = fullfile(PathMain,[ProjectName AnalysisName], 'OutputMats' ); mkdir(PathFilter);


Cell={{'Switch','Repeat',9,2,0.7,1.2,0.7,1.2,'All',-0.8,2.2,'V1'},{'Switch','Repeat',5.5,2,1.5,2.2,1.5,2.2,'All',-0.8,2.2,'V1'},{'Bon','Stan',6,2,0,0.6,0,0.6,'All',-0.8,2.2,'V1'},{'Bon','Stan',12,2,0.3,1.3,0.3,1.3,'All',-0.8,2.2,'V1'}}; 


%% Params end
for analysis=1:length(Cell)
    AnalysisCell=Cell{1,analysis};
    Conditions{1}=AnalysisCell{1,1};
    Conditions{2}=AnalysisCell{1,2};
    Freq=AnalysisCell{1,3};
    FreqRange=AnalysisCell{1,4};
    StartSec1=AnalysisCell{1,5};
    EndSec1=AnalysisCell{1,6};
    StartSec2=AnalysisCell{1,7};
    EndSec2=AnalysisCell{1,8};
    Conditions{3}=AnalysisCell{1,9};
    StartSecFil=AnalysisCell{1,10};
    EndSecFil=AnalysisCell{1,11};
    AnalysisVersion=AnalysisCell{1,12};  
    
    for subj=1:length(SubjList)
        %% Load MRI related data
        load(fullfile(PathSourcemodel,SubjList{subj})) % loaded is 'sourcemodel'
        load(fullfile(PathHeadmodel,SubjList{subj})) % loaded is 'headmodel'
        load(fullfile(PathElec,SubjList{subj})); % loaded is elec
       
    
        %% Load data and adjust elec to data structure
        if strcmpi(ProjectName,'Bocotilt')
            EEG=pop_loadset('filepath',PathPreProc,'filename',['VP' SubjList{1}(end-1:end) '_cleaned.set']);
            % Rename faulty electrode names
            elec.label(3)={'FCz'}; 
        else
            EEG=pop_loadset('filepath',PathPreProc,'filename',[SubjList{1} '.set']);
        end
        
        % Extract idx of chan order in data. idx can be applied on elec, strIdx makes sure the label (e.g. 'Fz') appears only once in elec.label :
        label={EEG.chanlocs.labels}; 
        chanIdx=nan(length(label),1);
        for chan=1:length(label)
            chanIdx(chan)=strIdx(elec.label,label{chan});
        end
        % shorten elec structure so that elec only includes the channels present in the data
        % Bocotilt elec has fields: unit ('mm', coordsys ('ctf'), label(cell), elecpos(double), chanpos(double), tra(127x127 double), cfg(struct):
        elec.chanpos=elec.chanpos(chanIdx,:); % double
        elec.elecpos=elec.elecpos(chanIdx,:); % double
        %elec.chantype=elec.chantype(chanIdx); % cell
        %elec.chanunit=elec.chanunit(chanIdx); % cell
        elec.label=elec.label(chanIdx); % cell
        elec.tra=elec.tra(chanIdx,chanIdx);
        if PlotHeadmodel
            figure
            ft_plot_mesh(headmodel.bnd(1), 'edgecolor','none','facealpha',0.8,'facecolor',[0.6 0.6 0.8]);
            hold on;
            ft_plot_sens(elec);
            savefig(fullfile(PathFilter,[SubjList{subj} '.fig']))
        end
    
        
        
        %% EEGLab 2 Fieldtrip
        % Convert EEGLAB data to fieldtrip data
        if strcmpi(ProjectName,'Bocotilt')
            EEG=pop_loadset('filepath',PathPreProc,'filename',['VP' SubjList{subj}(end-1:end) '_cleaned.set']);
        else
            EEG=pop_loadset('filepath',PathPreProc,'filename',[SubjList{subj} '.set']);
        end
        %EEG=pop_epoch(EEG,{'trialcue_real'},[0,2.6]);  % concentrate on cue interval
        dataRaw=eeglab2fieldtrip( EEG, 'raw', 'none' );  % all information here that I need for the leadfield
        
        % Load trialinfo 
        load( fullfile( PathTrialInfo,[CondName 'CondIdx2D' Version SubjList{subj}] ) ); % loaded is CondIdx2D: trial x cond, logical indices
        
        if ResponseLockedEpochs % Response-locked epochs and time intervals,this loop is Bocotil-specific!!
            for trial=1:size(dataRaw.trialinfo,1)
                % Step 1:
                % Find new time Point zero (time point of response onset): 
                %Rts=dataRaw.trialinfo{:,17}; % 17 = Rts
                Rts=EEG.trialinfo(:,17);  % better read from EEGlab data than from fieldtrip in case anything wnet wrong duirng transformation from eeglab 2fieldtrip structure
                Tp0=Rts(trial)/1000+0.8; 
                if or(isnan(Tp0), Tp0>2.0)  % missing answer or RT exceeds data range /trial goes until 2.2, I want to look at timings until 200 after response onset
                    disp([ SubjList{subj} ' trial:' num2str(trial) ' no RT in range'])
                    Tp0=0.5+0.8;  % take standard RT to have comparable time range
                end
                % Change time vector of that trial
                TimeOld=dataRaw.time{trial};
                TimeNew=TimeOld-Tp0; 
                % Step2 : 
                % % Alternative1
                %dataRaw.time{trial}=TimeNew; 
                %%%%%%%%%%%%%%%%%%%%% Alternative 2: %%%%%%%%%%%%%%%%%%%%%%%%%%
                % Cut out a time interval from -2 seconds until 200 ms after RTonset from data and time vector 
                Idx=dsearchn(TimeNew',-2):dsearchn(TimeNew',Tp0+0.2);
                dataRaw.time{1,trial}=TimeNew(Idx); 
                DataTemp=dataRaw.trial{1,trial}; % DataTemp: channel x time points
                dataRaw.trial{1,trial}=DataTemp(:,Idx);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
        end
    
    
        cfg = [];   cfg.trials = CondIdx2D(:,strIdx(CondLabel,Conditions{3}));
            cfg.toilim = [StartSecFil EndSecFil];
            data = ft_redefinetrial(cfg, dataRaw);
        cfg = []; cfg.trials = CondIdx2D(:,strIdx(CondLabel,Conditions{1})); 
            cfg.toilim = [StartSec1 EndSec1];
            cond1SgTr = ft_redefinetrial(cfg, dataRaw);
        cfg = []; cfg.trials = CondIdx2D(:,strIdx(CondLabel,Conditions{2})); 
            cfg.toilim = [StartSec2 EndSec2];
            cond2SgTr = ft_redefinetrial(cfg, dataRaw);
        
        % make sure to have the right electrode positions : 
        data.elec=elec;
        cond1SgTrial.elec=elec; 
        cond2SgTrial.elec=elec; 
        
    
        %% Calculate cross power spectral density
    
        cfg = [];
        cfg.method    = 'fft'; % fast fourier transformation
        cfg.output    = 'powandcsd';
        cfg.tapsmofrq = FreqRange;  
        cfg.foilim    = [Freq Freq]; 
        freqCond1 = ft_freqanalysis(cfg, cond1SgTr);
        
        cfg = [];
        cfg.method    = 'fft';
        cfg.output    = 'powandcsd';
        cfg.tapsmofrq = FreqRange;
        cfg.foilim    = [Freq Freq];
        freqCond2 = ft_freqanalysis(cfg, cond2SgTr);
    
        % dataAll = ft_appenddata([], dataPre, dataPost);
        cfg = [];
        cfg.method    = 'fft';
        cfg.output    = 'powandcsd';
        cfg.tapsmofrq = FreqRange;
        cfg.foilim    = [Freq Freq];
        freqAll = ft_freqanalysis(cfg, data);
    

        %% Source analysis
    
        % Compute filter
        cfg              = [];
        cfg.method       = 'dics';  % 'dics'
        cfg.frequency    = Freq;
        cfg.sourcemodel  = sourcemodel;
        cfg.headmodel    = headmodel;
        cfg.dics.projectnoise = 'yes'; % yes=is actually used when we do not compare conditions to reduce centre of head bias
        cfg.dics.lambda       = '5%';
        cfg.dics.keepfilter   = 'yes'; % filter must be saved of course, in order to apply it later to single conditions. *beamformerTut
        cfg.dics.realfilter   = 'yes';  
        sourceAll = ft_sourceanalysis(cfg, freqAll);  
            
        % Apply filter to conditions
        cfg.sourcemodel.filter = sourceAll.avg.filter;
    
        sourceCond1  = ft_sourceanalysis(cfg, freqCond1 );
        sourceCond2 = ft_sourceanalysis(cfg, freqCond2);
        

        % StanBonFreq3Range2Time0to800AllTime0to800V2
        PathSaving= fullfile(PathFilter,...
            [Conditions{1} Conditions{2} 'Freq' num2str(Freq*10) 'FreqRange' num2str(FreqRange*10) ...
            'Time1' num2str(StartSec1*1000) 'to' num2str(EndSec1*1000) 'Time2' num2str(StartSec2*1000) 'to' num2str(EndSec2*1000) ...
            Conditions{3} 'Time' num2str(StartSecFil*1000) 'to' num2str(EndSecFil*1000) AnalysisVersion]); 
        mkdir(PathSaving);
        save(fullfile(PathSaving,SubjList{subj}),'sourceAll','sourceCond1','sourceCond2');  % '-v7.3'
    end 
end
