% Make conditions for Bocotilt
clear all
close all
clc

% info from Stefan: 
% 1: id
% 2: block_nr
% 3: trial_nr
% 4: bonustrial
% 5: tilt_task
% 6: cue_ax (0 or 1)
% 7: target_red_left (0 or 1)
% 8: distractor_red_left (0 or 1)
% 9: response_interference
% 10: task_switch (-1,0,1) irgendeine Zuordnung zu switch , repeat, erster trial + deshalb unklar?
% 11: prev_switch (-1,0,1)
% 12: prev_accuracy (-1,0,1 or 2) (missing answer, no previous trial, correct, incorrect)
% 13: correct_response
% 14: response_side(0,1 or 2) 2 muss missing sein
% 15: rt
% 16: rt_thresh_color
% 17: rt_thresh_tilt
% 18: accuracy
% 19: position_color
% 20: position_tilt
% 21: position_target
% 22: position_distractor
% 23: sequence_position (Zahlen zwischen 1 + 8) --> sequenz innerhalb 8erBlock: trial 1-8


ProjectName='Bocotilt'; 
if isunix, PathMain='/media/sf_nathaExtern'; else PathMain = 'H:'; end
    PathEeglab = fullfile(PathMain,'Toolboxes','eeglab2020_0'); addpath(PathEeglab); eeglab;
    PathLib1 = fullfile(PathMain,'Toolboxes','FunctionLib'); addpath(genpath(PathLib1));
AnalysisName='PR010'; % analysis code of outcome of this script, if ='C004P103'; then P001 is loaded
    PathPreProcMat = fullfile(PathMain,[ProjectName AnalysisName],'OutputMats');
    PathPreProc = fullfile(PathMain,[ProjectName AnalysisName],'3 preprocessed'); 
    PathTrialInfo = fullfile(PathMain,[ProjectName AnalysisName], 'Trialinfo' ); 
%SubjList={'008','010','011','013','014','015','016','017','018','019','020','021','022','023','024','025','026','027','028','029','030','031','032','033','034'};
SubjList={'008','009','010','011','013','014','015','016','017','018','019','020','021','022','023','024','025','026','027','028','029','030','031','032','033','034'}; 


%%
NumTrialsPerBinV1=nan(56,4,length(SubjList));
NumTrialsPerBinV2=nan(56,4,length(SubjList));

for subj=1:length(SubjList)

% load trialinfo matrix: 
if strcmpi(AnalysisName(1:2),'PR')
    EEG=pop_loadset(['VP' SubjList{subj}(end-1:end) '_cleaned.set'],PathPreProc);
    Ti=EEG.trialinfo;  
    disp(['loading' SubjList{subj}]);
else % normally not used
    EEG=pop_loadset(['VP' SubjList{subj}(end-1:end) '_cleaned.set'],fullfile(PathMain,[ProjectName AnalysisName(6:end)],'3 preprocessed')); 
    Ti=EEG.trialinfo;  
end





CondName='ForDecoding'; Version='_V1_'; 
        % Calculate num trials per condition
        CondLabel={'BonTask1','BonTask2','StanTask1','StanTask2'};
     
        NoBreak=Ti(:, 2) > 4; 
        BonTask1=and( NoBreak, and( Ti(:,4)==1, Ti(:,5)==0) ) ;
        BonTask2=and( NoBreak, and( Ti(:,4)==1, Ti(:,5)==1) ) ;
        StanTask1=and( NoBreak, and( Ti(:,4)==0, Ti(:,5)==0) ) ;
        StanTask2=and( NoBreak, and( Ti(:,4)==0, Ti(:,5)==1) ) ;

        
        clear CondIdx2D;
        CondIdx2D(:,1)=BonTask1; CondIdx2D(:,2)=BonTask2; CondIdx2D(:,3)=StanTask1; CondIdx2D(:,4)=StanTask2; 
        save(fullfile(PathTrialInfo,[CondName 'CondIdx2D' Version SubjList{subj}]),'CondIdx2D');
        NumTrialsPerBinV1(:,:,subj)=calcNumTrials(CondIdx2D,CondLabel,Ti);


CondName='ForDecoding'; Version='_V2_';  % gleich viele Trials in allen 4 Bedingungen
        % Calculate num trials per condition
        CondLabel={'BonTask1','BonTask2','StanTask1','StanTask2'};
        NoBreak=Ti(:, 2) > 4; 
        BonTask1=and( NoBreak, and( Ti(:,4)==1, Ti(:,5)==0) ) ;
        BonTask2=and( NoBreak, and( Ti(:,4)==1, Ti(:,5)==1) ) ;
        StanTask1=and( NoBreak, and( Ti(:,4)==0, Ti(:,5)==0) ) ;
        StanTask2=and( NoBreak, and( Ti(:,4)==0, Ti(:,5)==1) ) ;
              
        clear CondIdx2D;
        CondIdx2D(:,1)=BonTask1; CondIdx2D(:,2)=BonTask2; CondIdx2D(:,3)=StanTask1; CondIdx2D(:,4)=StanTask2; 
        
        [TrialIndsReduced, CondIdx2D]=getSameAmountTrials(CondIdx2D, CondLabel);
        
        save(fullfile(PathTrialInfo,[CondName 'CondIdx2D' Version SubjList{subj}]),'CondIdx2D');
        NumTrialsPerBinV2(:,:,subj)=calcNumTrials(CondIdx2D,CondLabel,Ti);

CondName='LikeSensorLevel'; Version='_V1_'; 
        CondLabel={'Bon','Stan','Switch','Repeat','LeftAnswer','RightAnswer','Color','Ori'}; 
        clear CondIdx2D;
        %EEG=pop_loadset(['VP' SubjList{subj}(end-1:end) '_cleaned.set'],PathPreProc);
        %load(fullfile(PathPreProcMat, SubjList{subj}));  % loaded is EEG
        NoBreak=Ti(:, 2) > 4;
%         Task1=Ti(:,)
%         Task2=
        
        % Exclude trials
        % 1: id % 2: block_nr % 3: trial_nr % 4: bonustrial % 5: tilt_task % 6: cue_ax % 7: target_red_left % 8: distractor_red_left
         % 9: response_interference % 10: task_switch % 11: prev_switch % 12: prev_accuracy % 13: correct_response % 14: response_side (0,1 or2)
        % 15: rt % 16: rt_thresh_color % 17: rt_thresh_tilt % 18: accuracy % 19: position_color % 20: position_tilt % 21: position_target
        % 22: position_distractor % 23: sequence_position
        CondIdx2D(:,1)=and(Ti(:, 2) > 4 ,Ti(:, 23) > 1)  ; %  break trials and first trials of block excluded
        CondIdx2D(:,2)=and( Ti(:, 4) ==1 , and(Ti(:, 2) > 4 ,Ti(:, 23) > 1) ) ;
        CondIdx2D(:,3)=and( Ti(:, 4) ==0 , and(Ti(:, 2) > 4 ,Ti(:, 23) > 1) ) ;
        CondIdx2D(:,4)=and( Ti(:, 10) ==1 , and(Ti(:, 2) > 4 ,Ti(:, 23) > 1) ) ;
        CondIdx2D(:,5)=and( Ti(:, 10) ==0 , and(Ti(:, 2) > 4 ,Ti(:, 23) > 1) ) ;
        CondIdx2D(:,6)=and( Ti(:, 14) ==0 , and(Ti(:, 2) > 4 ,Ti(:, 23) > 1) ) ; % Stefan: 0=left, 1=right, 2 = no answer
        CondIdx2D(:,7)=and( Ti(:, 14) ==1 , and(Ti(:, 2) > 4 ,Ti(:, 23) > 1) ) ;
        save(fullfile(PathTrialInfo,[CondName 'CondIdx2D' Version SubjList{subj}]),'CondIdx2D');

% CondName='LikeSensorLevel'; Version='_V1_'; 
%         CondLabel={'All','Bon','Stan','Switch','Repeat'}; 
%         clear CondIdx2D;
%         %EEG=pop_loadset(['VP' SubjList{subj}(end-1:end) '_cleaned.set'],PathPreProc);
%             load(fullfile(PathPreProcMat, SubjList{subj}));  % loaded is EEG
%         disp(['loading' SubjList{subj}]);
%         % Exclude trials
%         Ti=EEG.trialinfo; 
%         % 1: id % 2: block_nr % 3: trial_nr % 4: bonustrial % 5: tilt_task % 6: cue_ax % 7: target_red_left % 8: distractor_red_left
%          % 9: response_interference % 10: task_switch % 11: prev_switch % 12: prev_accuracy % 13: correct_response % 14: response_side (0,1 or2)
%         % 15: rt % 16: rt_thresh_color % 17: rt_thresh_tilt % 18: accuracy % 19: position_color % 20: position_tilt % 21: position_target
%         % 22: position_distractor % 23: sequence_position
%         CondIdx2D(:,1)=and(Ti(:, 2) > 4 ,Ti(:, 23) > 1)  ; %  break trials and first trials of block excluded
%         CondIdx2D(:,2)=and( Ti(:, 4) ==1 , and(Ti(:, 2) > 4 ,Ti(:, 23) > 1) ) ;
%         CondIdx2D(:,3)=and( Ti(:, 4) ==0 , and(Ti(:, 2) > 4 ,Ti(:, 23) > 1) ) ;
%         CondIdx2D(:,4)=and( Ti(:, 10) ==1 , and(Ti(:, 2) > 4 ,Ti(:, 23) > 1) ) ;
%         CondIdx2D(:,5)=and( Ti(:, 10) ==0 , and(Ti(:, 2) > 4 ,Ti(:, 23) > 1) ) ;
%         save(fullfile(PathTrialInfo,[CondName 'CondIdx2D' Version SubjList{subj}]),'CondIdx2D');
CondName='Reward'; Version='_V1_';
        CondLabel={'Bon','Stan'}; 
        clear CondIdx2D;
        %EEG=pop_loadset(['VP' SubjList{subj}(end-1:end) '_cleaned.set'],PathPreProc);
        load(fullfile(PathPreProcMat, SubjList{subj}));  % loaded is EEG
        disp(['loading' SubjList{subj}]);
        % Exclude trials
        % 1: id % 2: block_nr % 3: trial_nr % 4: bonustrial % 5: tilt_task % 6: cue_ax % 7: target_red_left % 8: distractor_red_left
         % 9: response_interference % 10: task_switch % 11: prev_switch % 12: prev_accuracy % 13: correct_response % 14: response_side
        % 15: rt % 16: rt_thresh_color % 17: rt_thresh_tilt % 18: accuracy % 19: position_color % 20: position_tilt % 21: position_target
        % 22: position_distractor % 23: sequence_position
        CondIdx2D(:,1)=and( Ti(:, 4) ==1 , and(Ti(:, 2) > 4 ,Ti(:, 23) > 1) ) ;
        CondIdx2D(:,2)=and( Ti(:, 4) ==0 , and(Ti(:, 2) > 4 ,Ti(:, 23) > 1) ) ;
        TrialIndsReduced=getSameAmountTrials(CondIdx2D , CondLabel);
        save(fullfile(PathTrialInfo,[CondName 'TrialIndsReduced' Version SubjList{subj}]),'TrialIndsReduced');
        save(fullfile(PathTrialInfo,[CondName 'CondIdx2D' Version SubjList{subj}]),'CondIdx2D');

CondName='Transition'; Version='_V1_'; 
        CondLabel={'Switch','Repeat'}; 
        clear CondIdx2D;
        %EEG=pop_loadset(['VP' SubjList{subj}(end-1:end) '_cleaned.set'],PathPreProc);
        load(fullfile(PathPreProcMat, SubjList{subj}));  % loaded is EEG
        disp(['loading' SubjList{subj}]);
        % Exclude trials
        % 1: id % 2: block_nr % 3: trial_nr % 4: bonustrial % 5: tilt_task % 6: cue_ax % 7: target_red_left % 8: distractor_red_left
         % 9: response_interference % 10: task_switch % 11: prev_switch % 12: prev_accuracy % 13: correct_response % 14: response_side
        % 15: rt % 16: rt_thresh_color % 17: rt_thresh_tilt % 18: accuracy % 19: position_color % 20: position_tilt % 21: position_target
        % 22: position_distractor % 23: sequence_position
        CondIdx2D(:,1)=and( Ti(:, 10) ==1 , and(Ti(:, 2) > 4 ,Ti(:, 23) > 1) ) ;
        CondIdx2D(:,2)=and( Ti(:, 10) ==0 , and(Ti(:, 2) > 4 ,Ti(:, 23) > 1) ) ;
        TrialIndsReduced=getSameAmountTrials(CondIdx2D , CondLabel);
        save(fullfile(PathTrialInfo,[CondName 'TrialIndsReduced' Version SubjList{subj}]),'TrialIndsReduced');
        save(fullfile(PathTrialInfo,[CondName 'CondIdx2D' Version SubjList{subj}]),'CondIdx2D');

  

   
end


%% Function definitions
% function - make sure all conditions have the same amount of trials 
function [TrialIndsReduced, CondIdx2DNew]=getSameAmountTrials(CondIdx2D, CondLabel)
    % Calculate number of trials per condition 
    nTrialsPerCond=nan(1,length(CondLabel));
    for cond=1:length(CondLabel)
        nTrialsPerCond(cond)=sum(CondIdx2D(:,cond));
    end

%     % Randomly draw min(nTrialsPerCond) trials for each condition (calulcate output, whereas output has integers)
%     TrialIndsReduced = nan( min(nTrialsPerCond) , length(CondLabel) );
%     for cond=1:length(CondLabel)
%         TrialInds=find(CondIdx2D(:,cond));
%         % p = randperm(n,k) returns a row vector containing k unique integers selected randomly from 1 to n
%         TrialIndsReduced(:,cond) = sort(TrialInds(randperm( nTrialsPerCond(cond) , min(nTrialsPerCond) )));
%     end

    % Randomly draw trials that shall be removed (calculate output, whereas output has logicals)
    HowManyTrialsToRemove = nan(1,length(CondLabel));
    for cond=1:length(CondLabel)
        HowManyTrialsToRemove(cond)=nTrialsPerCond(cond) - min(nTrialsPerCond);
    end
    % Randomly draw those trials that shall be removed
    %IdxRemovedTrials = nan(1,length(CondLabel));
    TrialIndsReduced = nan( min(nTrialsPerCond) , length(CondLabel) );
    CondIdx2DNew = CondIdx2D;
    for cond=1:length(CondLabel)
        TrialInds=find(CondIdx2D(:,cond));
        % p = randperm(n,k) returns a row vector containing k unique integers selected randomly from 1 to n
        IdxRemovedTrials = sort(TrialInds(randperm( nTrialsPerCond(cond) , HowManyTrialsToRemove(cond) ) ));
        CondIdx2DNew( IdxRemovedTrials ,cond)=false(1,length(IdxRemovedTrials));
        TrialIndsReduced(:,cond) = find(CondIdx2DNew(:,cond));
    end


end



function NumTrialsPerBin=calcNumTrials(CondIdx2D,CondLabel,Ti)  % specific for bocotilt!!!
    % how many trials per cond? copy and paste in Excel
        NumTrialsPerBin=nan(56,length(CondLabel));
        for cond=1:length(CondLabel)
             % 1: id % 2: block_nr % 3: trial_nr % 4: bonustrial % 5: tilt_task % 6: cue_ax % 7: target_red_left % 8: distractor_red_left
             % 9: response_interference % 10: task_switch % 11: prev_switch % 12: prev_accuracy % 13: correct_response % 14: response_side (0,1 or2)
            % 15: rt % 16: rt_thresh_color % 17: rt_thresh_tilt % 18: accuracy % 19: position_color % 20: position_tilt % 21: position_target
            % 22: position_distractor % 23: sequence_position
            All=CondIdx2D(:,cond); NumTrialsPerBin(1,cond)=sum(All); 
            
            Switch=and(CondIdx2D(:,cond),Ti(:,10)==1); NumTrialsPerBin(2,cond)=sum(Switch);
            Switch=and(CondIdx2D(:,cond),Ti(:,10)==0); NumTrialsPerBin(3,cond)=sum(Switch);
            Switch=and(CondIdx2D(:,cond),Ti(:,10)==-1); NumTrialsPerBin(4,cond)=sum(Switch);
            
            PrevSwitch=and(CondIdx2D(:,cond),Ti(:,11)==1); NumTrialsPerBin(5,cond)=sum(PrevSwitch);
            PrevSwitch=and(CondIdx2D(:,cond),Ti(:,11)==0); NumTrialsPerBin(6,cond)=sum(PrevSwitch);
            PrevSwitch=and(CondIdx2D(:,cond),Ti(:,11)==-1); NumTrialsPerBin(7,cond)=sum(PrevSwitch);
            
            PrevCorr=and(CondIdx2D(:,cond),Ti(:,12)==-1); NumTrialsPerBin(8,cond)=sum(PrevCorr);
            PrevCorr=and(CondIdx2D(:,cond),Ti(:,12)==0); NumTrialsPerBin(9,cond)=sum(PrevCorr);
            PrevCorr=and(CondIdx2D(:,cond),Ti(:,12)==1); NumTrialsPerBin(10,cond)=sum(PrevCorr);
            PrevCorr=and(CondIdx2D(:,cond),Ti(:,12)==2); NumTrialsPerBin(11,cond)=sum(PrevCorr);
            
            TargetRedOrLeft=and(CondIdx2D(:,cond),Ti(:,7)==0); NumTrialsPerBin(12,cond)=sum(TargetRedOrLeft);
            TargetRedOrLeft=and(CondIdx2D(:,cond),Ti(:,7)==1); NumTrialsPerBin(13,cond)=sum(TargetRedOrLeft);
            
            DistractorRedOrLeft=and(CondIdx2D(:,cond),Ti(:,8)==0); NumTrialsPerBin(14,cond)=sum(DistractorRedOrLeft);
            DistractorRedOrLeft=and(CondIdx2D(:,cond),Ti(:,8)==1); NumTrialsPerBin(15,cond)=sum(DistractorRedOrLeft);
            
            CueAx=and(CondIdx2D(:,cond),Ti(:,6)==0); NumTrialsPerBin(16,cond)=sum(CueAx);
            CueAx=and(CondIdx2D(:,cond),Ti(:,6)==1); NumTrialsPerBin(17,cond)=sum(CueAx);

            Block=and(CondIdx2D(:,cond),Ti(:,2)==4); NumTrialsPerBin(18,cond)=sum(Block);
            Block=and(CondIdx2D(:,cond),Ti(:,2)==5); NumTrialsPerBin(19,cond)=sum(Block);
            Block=and(CondIdx2D(:,cond),Ti(:,2)==6); NumTrialsPerBin(20,cond)=sum(Block);
            Block=and(CondIdx2D(:,cond),Ti(:,2)==7); NumTrialsPerBin(21,cond)=sum(Block);
            Block=and(CondIdx2D(:,cond),Ti(:,2)==8); NumTrialsPerBin(22,cond)=sum(Block);
            Block=and(CondIdx2D(:,cond),Ti(:,2)==9); NumTrialsPerBin(23,cond)=sum(Block);
            Block=and(CondIdx2D(:,cond),Ti(:,2)==10); NumTrialsPerBin(24,cond)=sum(Block);

            ColorLoc=and(CondIdx2D(:,cond),Ti(:,19)==1); NumTrialsPerBin(25,cond)=sum(ColorLoc);
            ColorLoc=and(CondIdx2D(:,cond),Ti(:,19)==2); NumTrialsPerBin(26,cond)=sum(ColorLoc);
            ColorLoc=and(CondIdx2D(:,cond),Ti(:,19)==3); NumTrialsPerBin(27,cond)=sum(ColorLoc);
            ColorLoc=and(CondIdx2D(:,cond),Ti(:,19)==4); NumTrialsPerBin(28,cond)=sum(ColorLoc);
            ColorLoc=and(CondIdx2D(:,cond),Ti(:,19)==5); NumTrialsPerBin(29,cond)=sum(ColorLoc);
            ColorLoc=and(CondIdx2D(:,cond),Ti(:,19)==6); NumTrialsPerBin(30,cond)=sum(ColorLoc);
            ColorLoc=and(CondIdx2D(:,cond),Ti(:,19)==7); NumTrialsPerBin(31,cond)=sum(ColorLoc);
            ColorLoc=and(CondIdx2D(:,cond),Ti(:,19)==8); NumTrialsPerBin(32,cond)=sum(ColorLoc);

            TiltLoc=and(CondIdx2D(:,cond),Ti(:,20)==1); NumTrialsPerBin(33,cond)=sum(TiltLoc);
            TiltLoc=and(CondIdx2D(:,cond),Ti(:,20)==2); NumTrialsPerBin(34,cond)=sum(TiltLoc);
            TiltLoc=and(CondIdx2D(:,cond),Ti(:,20)==3); NumTrialsPerBin(35,cond)=sum(TiltLoc);
            TiltLoc=and(CondIdx2D(:,cond),Ti(:,20)==4); NumTrialsPerBin(36,cond)=sum(TiltLoc);
            TiltLoc=and(CondIdx2D(:,cond),Ti(:,20)==5); NumTrialsPerBin(37,cond)=sum(TiltLoc);
            TiltLoc=and(CondIdx2D(:,cond),Ti(:,20)==6); NumTrialsPerBin(38,cond)=sum(TiltLoc);
            TiltLoc=and(CondIdx2D(:,cond),Ti(:,20)==7); NumTrialsPerBin(39,cond)=sum(TiltLoc);
            TiltLoc=and(CondIdx2D(:,cond),Ti(:,20)==8); NumTrialsPerBin(40,cond)=sum(TiltLoc);

            TargetLoc=and(CondIdx2D(:,cond),Ti(:,21)==1); NumTrialsPerBin(41,cond)=sum(TargetLoc);
            TargetLoc=and(CondIdx2D(:,cond),Ti(:,21)==2); NumTrialsPerBin(42,cond)=sum(TargetLoc);
            TargetLoc=and(CondIdx2D(:,cond),Ti(:,21)==3); NumTrialsPerBin(43,cond)=sum(TargetLoc);
            TargetLoc=and(CondIdx2D(:,cond),Ti(:,21)==4); NumTrialsPerBin(44,cond)=sum(TargetLoc);
            TargetLoc=and(CondIdx2D(:,cond),Ti(:,21)==5); NumTrialsPerBin(45,cond)=sum(TargetLoc);
            TargetLoc=and(CondIdx2D(:,cond),Ti(:,21)==6); NumTrialsPerBin(46,cond)=sum(TargetLoc);
            TargetLoc=and(CondIdx2D(:,cond),Ti(:,21)==7); NumTrialsPerBin(47,cond)=sum(TargetLoc);
            TargetLoc=and(CondIdx2D(:,cond),Ti(:,21)==8); NumTrialsPerBin(48,cond)=sum(TargetLoc);
        
            DistractorLoc=and(CondIdx2D(:,cond),Ti(:,22)==1); NumTrialsPerBin(49,cond)=sum(DistractorLoc);
            DistractorLoc=and(CondIdx2D(:,cond),Ti(:,22)==2); NumTrialsPerBin(50,cond)=sum(DistractorLoc);
            DistractorLoc=and(CondIdx2D(:,cond),Ti(:,22)==3); NumTrialsPerBin(51,cond)=sum(DistractorLoc);
            DistractorLoc=and(CondIdx2D(:,cond),Ti(:,22)==4); NumTrialsPerBin(52,cond)=sum(DistractorLoc);
            DistractorLoc=and(CondIdx2D(:,cond),Ti(:,22)==5); NumTrialsPerBin(53,cond)=sum(DistractorLoc);
            DistractorLoc=and(CondIdx2D(:,cond),Ti(:,22)==6); NumTrialsPerBin(54,cond)=sum(DistractorLoc);
            DistractorLoc=and(CondIdx2D(:,cond),Ti(:,22)==7); NumTrialsPerBin(55,cond)=sum(DistractorLoc);
            DistractorLoc=and(CondIdx2D(:,cond),Ti(:,22)==8); NumTrialsPerBin(56,cond)=sum(DistractorLoc);
        end
end