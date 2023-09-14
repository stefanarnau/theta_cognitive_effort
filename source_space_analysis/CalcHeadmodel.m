%% MRI preprocessing - create headmodel 
% cp. tutorial: https://www.fieldtriptoolbox.org/tutorial/headmodel_eeg_bem/

% by NL (nathalie.liegel@gmail.com), last change: 10.09.22

clear all
close all
clc


%% Parameters
LoadSegementation=1; % if 1, segemnted mri is loaded , otherwise mesh is loaded
    % only neede if LoadSegmentation=1:
    NumVertices = [3000 2000 1000]; % brain skull scalp, tutorial: = [3000 2000 1000];
scalpsmooth=10; 
scalpthreshold=0.07; 
Method='dipoli';  % Create a volume conduction model using Method='dipoli', 'openmeeg', or 'bemcp'
AnalysisSteps='HM001BN001SG001AL002'; 

%ProjectName='Auvibo';  
    %SubjList={'005','006','007','008','009','010','011'}; 
    %SubjList={'011','010'};
ProjectName='Bocotilt'; % '010','011','013','014','015','016','017'
    SubjList={'010','011','013','014','015','016','017','018','019','020','021','023','024','025','026','027','028','029','030','031','032','034'}; % 022
    %SubjList={'022'};  % for single subject plotting and param determination

if isunix, PathMain='/media/sf_nathaExtern'; else, PathMain = 'H:'; end  % alternatively: ispc for windows
    PathFT = fullfile(PathMain, 'Toolboxes','fieldtrip-20220104'); addpath(PathFT);ft_defaults;
    PathMriSeg=fullfile(PathMain,[ProjectName 'DataBase'],AnalysisSteps(11:end),['thre' num2str(scalpthreshold)],['smooth' num2str(scalpsmooth)]);
    PathMriMesh=fullfile(PathMain,[ProjectName 'DataBase'],AnalysisSteps(6:end),['thre' num2str(scalpthreshold)],['smooth' num2str(scalpsmooth)]);
    PathHeadmodel=fullfile(PathMain,[ProjectName 'DataBase'],AnalysisSteps); mkdir(PathHeadmodel);

%% 
for subj=1:length(SubjList)
    
    if LoadSegementation
        load (fullfile(PathMriSeg,[SubjList{subj} '.mat'])); % loaded is segmentedmri
        % Create mesh:
        cfg=[];
        cfg.tissue={'brain','skull','scalp'};
        cfg.numvertices = NumVertices;
        bnd=ft_prepare_mesh(cfg,segmentedmri);
    else
        load (fullfile(PathMriMesh,SubjList{subj})); % loaded is bnd
    end


    cfg        = [];
    cfg.method =Method; 
    headmodel  = ft_prepare_headmodel(cfg, bnd);


    % Save
    save(fullfile(PathHeadmodel,SubjList{subj}),'headmodel'); 

end