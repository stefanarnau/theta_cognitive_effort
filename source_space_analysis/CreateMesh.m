%% MRI preprocessing - Create meshes and plot
% cp. tutorial: https://www.fieldtriptoolbox.org/tutorial/headmodel_eeg_bem/

% by NL (nathalie.liegel@gmail.com), last change: 10.09.22

clear all
close all
clc


%% Parameter
PlotSgSubj=0; 
SavePlots=1; 

%ProjectName='Auvibo';  
    %SubjList={'005','006','007','008','009','010','011'}; 
    %SubjList={'011','010'};
ProjectName='Bocotilt'; % '010','011','013','014','015','016','017'
    SubjList={'010','011','013','014','015','016','017','018','019','020','021','023','024','025','026','027','028','029','030','031','032','034'};  % '022',
    %SubjList={'022','026'};  

if isunix, PathMain='/media/sf_nathaExtern'; else, PathMain = 'H:'; end  % alternatively ispc for windows
    PathFT = fullfile(PathMain, 'Toolboxes','fieldtrip-20220104'); addpath(PathFT); ft_defaults;

NumVertices = [3000 2000 1000]; % brain skull scalp, tutorial: = [3000 2000 1000];

for scalpsmooth=10
for scalpthreshold=0.07  % to do 0.07, 0.1, 0.3

AnalysisSteps='BN001SG001AL002'; 
    PathMriAligned=fullfile(PathMain,[ProjectName 'DataBase'],AnalysisSteps(11:end)); % needed to get coordinates of single voxels
    PathMriSeg=fullfile(PathMain,[ProjectName 'DataBase'],AnalysisSteps(6:end),['thre' num2str(scalpthreshold)],['smooth' num2str(scalpsmooth)]);
    PathMriMesh=fullfile(PathMain,[ProjectName 'DataBase'],AnalysisSteps,['thre' num2str(scalpthreshold)],['smooth' num2str(scalpsmooth)]);mkdir(PathMriMesh);
    PathPlots=fullfile(PathMain,[ProjectName 'DataBase'],AnalysisSteps,['thre' num2str(scalpthreshold)]); mkdir(PathPlots);


%% 
for subj=1:length(SubjList)

    load(fullfile(PathMriSeg,SubjList{subj})); % loaded is segmentedmri

    % Create mesh
    cfg=[];
    cfg.tissue={'brain','skull','scalp'};
    cfg.numvertices = NumVertices;
    bnd=ft_prepare_mesh(cfg,segmentedmri);

    load(fullfile(PathMriAligned,SubjList{subj})); segmentedmri.anatomy=mri.anatomy; 

    if SavePlots
        % Segmented mri
        figure(1); clf
        cfg=[]; cfg.funparameter = 'brain'; ft_sourceplot(cfg,segmentedmri); savefig(fullfile(PathPlots,['BrainSeg' SubjList{subj} 'Scalpsmooth' num2str(scalpsmooth)]))
        figure(2); clf
        cfg=[]; cfg.funparameter = 'skull'; ft_sourceplot(cfg,segmentedmri); savefig(fullfile(PathPlots,['SkullSeg' SubjList{subj} 'Scalpsmooth' num2str(scalpsmooth)]))
        figure(3); clf
        cfg=[]; cfg.funparameter = 'scalp'; ft_sourceplot(cfg,segmentedmri); savefig(fullfile(PathPlots,['ScalpSeg' SubjList{subj} 'Scalpsmooth' num2str(scalpsmooth)]))

        % Mesh 
        figure(4);clf
        ft_plot_mesh(bnd(3),'facecolor','none'); 
        savefig(fullfile(PathPlots,['Scalpmesh' SubjList{subj} 'Scalpsmooth' num2str(scalpsmooth)])) %scalp
        figure(5);clf
        ft_plot_mesh(bnd(2),'facecolor','none'); 
        savefig(fullfile(PathPlots,['Skullmesh' SubjList{subj} 'Scalpsmooth' num2str(scalpsmooth)])) %skull
        figure(6);clf
        ft_plot_mesh(bnd(1),'facecolor','none'); 
        savefig(fullfile(PathPlots,['Brainmesh' SubjList{subj} 'Scalpsmooth' num2str(scalpsmooth)])) %brain
    end
    
    if PlotSgSubj  % plot without saving
        % All three meshes in one figure: 
        figure(); clf
        % brain: % facecolor=color von Fläche zwischen kanten, edgecolor= color von punkten und Kanten
        ft_plot_mesh(bnd(1), 'facecolor','r', 'facealpha', 0.1,...
            'edgecolor', [1 1 1], 'edgealpha', 0.5); % ,'maskstyle','colormix', 'facecolor',[0.2 0.2 0.2],'edgealpha', 0.05
        hold on;
        ft_plot_mesh(bnd(2),'facecolor','b','edgecolor','none','facealpha',0.2); % skull
        hold on;
        ft_plot_mesh(bnd(3),'edgecolor','none','facecolor',[0.4 0.6 0.4],'facealpha', 0.1);  % scalp

        % 
    end

    save(fullfile(PathMriMesh,SubjList{subj}),'bnd');
end
end
end






