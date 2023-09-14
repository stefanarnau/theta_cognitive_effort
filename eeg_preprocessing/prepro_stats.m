clear all;

% PATH VARS
PATH_EEGLAB        = 'insert_path_here';
PATH_AUTOCLEANED   = 'insert_path_here';

% Subject list
subject_list = {'VP09', 'VP17', 'VP25', 'VP10', 'VP11', 'VP13', 'VP14', 'VP15', 'VP16', 'VP18',...
                'VP19', 'VP20', 'VP21', 'VP22', 'VP23', 'VP08', 'VP24', 'VP26', 'VP27', 'VP28',...
                'VP29', 'VP30', 'VP31', 'VP32', 'VP33', 'VP34'};

% Init eeglab
addpath(PATH_EEGLAB);
eeglab;

% Preprostats matrix
prepro_stats = [];

% Iterate subjects
for s = 1 : length(subject_list)

    % Load info
    EEG = pop_loadset('filename', [subject_list{s} '_cleaned.set'], 'filepath', PATH_AUTOCLEANED, 'loadmode', 'info');

    prepro_stats(s, 1) = length(EEG.chans_rejected);
    prepro_stats(s, 2) = length(EEG.rejected_epochs);
    prepro_stats(s, 3) = length(EEG.nobrainer);

end

% Get descriptive statistics 
mean(prepro_stats, 1)
std(prepro_stats, [], 1)