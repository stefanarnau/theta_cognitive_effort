%% Align electrode positions to headmodel and sanity check plots

clear all
close all
clc


%% Parameters 
ProjectName='Bocotilt';  % = 'Bocotilt' or ='Auvibo';
%SubjStr='026';   %  ='005';='006';='010';etc
SubjList={'010','011','013','014','015','016','017','018','019','020','021','022','023','024','025','026','027','028','029','030','031','032','034'}; % '008','033'(nodipoli headmodel, bad scan)

AnalysisSteps='EL001HM001BN001SG001AL002'; scalpthreshold=0.07; scalpsmooth=10; 
PathMain='H:';  %='/media/sf_nathaExtern';  ='H:';
    PathFT = fullfile(PathMain, 'Toolboxes','fieldtrip-20220104'); addpath(PathFT);ft_defaults;
    PathScans=fullfile(PathMain,[ProjectName 'DataBase'],'2 scanDone');
    PathMriAligned=fullfile(PathMain,[ProjectName 'DataBase'],AnalysisSteps(end-4:end));
    PathMriSeg=fullfile(PathMain,[ProjectName 'DataBase'],AnalysisSteps(end-9:end),['thre' num2str(scalpthreshold)],['smooth' num2str(scalpsmooth)]);
    PathMriMesh=fullfile(PathMain,[ProjectName 'DataBase'],AnalysisSteps(end-14:end),['thre' num2str(scalpthreshold)],['smooth' num2str(scalpsmooth)]);
    PathHeadmodel=fullfile(PathMain,[ProjectName 'DataBase'],AnalysisSteps(end-19:end));
    PathElec=fullfile(PathMain,[ProjectName 'DataBase'],AnalysisSteps); mkdir(PathElec); 


if 0
%% Load
load(fullfile(PathHeadmodel,SubjStr)) % loaded is 'headmodel'
%load(fullfile(PathMriMesh,SubjStr)) % loaded is 'bnd'
load(fullfile(PathMriSeg,SubjStr)) % loaded is 'segmentedmri'
load (fullfile(PathScans,['elec' SubjStr])); % loaded is elec
load (fullfile(PathScans, SubjStr)); % loaded is head_surface_new
load (fullfile(PathMriAligned, SubjStr)); % loaded is mri


%% Optional: check MRI
% cfg=[]; ft_sourceplot(cfg,mri);


%% Sanity check on segmentation
% figure; cfg=[]; cfg.funparameter = 'brain'; ft_sourceplot(cfg,segmentedmri);
% figure; cfg=[]; cfg.funparameter = 'skull'; ft_sourceplot(cfg,segmentedmri);
% figure; cfg=[]; cfg.funparameter = 'scalp'; ft_sourceplot(cfg,segmentedmri);
% 
% 
% %https://www.fieldtriptoolbox.org/example/coregistration_quality_control/
% % [0 0 60] is good 
% figure; ft_plot_ortho(mri.anatomy, 'transform', mri.transform, 'location', [0 0 0], 'intersectmesh', headmodel.bnd(1))
% figure; ft_plot_montage(mri.anatomy, 'transform', mri.transform, 'intersectmesh', headmodel.bnd(1))
% figure; ft_plot_ortho(mri.anatomy, 'transform', mri.transform, 'location', [0 0 0], 'intersectmesh', headmodel.bnd(2))
% figure; ft_plot_montage(mri.anatomy, 'transform', mri.transform, 'intersectmesh', headmodel.bnd(2))
% figure; ft_plot_ortho(mri.anatomy, 'transform', mri.transform, 'location', [0 0 0], 'intersectmesh', headmodel.bnd(3))
% figure; ft_plot_montage(mri.anatomy, 'transform', mri.transform, 'intersectmesh', headmodel.bnd(3))


%% Sanity check on headmodel
figure;
ft_plot_mesh(headmodel.bnd(3),'facecolor','none'); %scalp
figure;
ft_plot_mesh(headmodel.bnd(2),'facecolor','none'); %skull
figure;
ft_plot_mesh(headmodel.bnd(1),'facecolor','none'); %scalp


figure(); clf
% scalp: % facecolor=color of surface between edges, edgecolor= color of points and edges
ft_plot_mesh(headmodel.bnd(1), 'facecolor', [0.2 0.2 0.2], 'facealpha', 1,...
    'edgecolor',  [0.2 0.2 0.2], 'edgealpha', 0.1); % ,'maskstyle','colormix', 'facecolor',[0.2 0.2 0.2],'edgealpha', 0.05, 'facecolor','r',
hold on;
ft_plot_mesh(headmodel.bnd(2),'edgecolor', [0.4 0.4 0.4],'edgealpha',0.4, 'facecolor','g', 'facealpha',0.4); % skull
hold on;
ft_plot_mesh(headmodel.bnd(3),'facecolor','b','facealpha', 0.1);  % scalp


%% Sanity check on electrodes
% figure(), clf
% ft_plot_mesh(head_surface_new)
% hold on
% ft_plot_sens(elec,'label','on')


%% Align electrode positions
% Plot electrodes
figure();clf
ft_plot_mesh(headmodel.bnd(3), 'edgecolor','none','facealpha',0.8,'facecolor',[0.6 0.6 0.8]);
hold on;
ft_plot_sens(elec)

% % Electrodes can all be moved inward if necessary, see tut: 
% cfg = [];
% cfg.method     = 'moveinward';
% cfg.moveinward = 5; % 12=12 mm
% cfg.elec       = elec;
% elec_aligned = ft_electroderealign(cfg);


% % electrodes can be aligned to MRI
% % see lower part of https://www.fieldtriptoolbox.org/tutorial/headmodel_eeg_bem/
% % option 1: automatic alignement (does not really make sense, fiducials have alredy been aligned earlier
% % option 2: manual alignment:
% cfg           = [];
% cfg.method    = 'interactive';
% cfg.elec      = elec;
% cfg.headshape = headmodel.bnd(3);
% elec_aligned  = ft_electroderealign(cfg);
% 
% % Plot electrodes
% figure();clf
% ft_plot_mesh(headmodel.bnd(1), 'edgecolor','none','facealpha',0.8,'facecolor',[0.6 0.6 0.8]);
% hold on;
% %ft_plot_sens(elec_aligned)
% ft_plot_sens(elec_aligned,'label','on')

elec=elec_aligned; 
%save(fullfile(PathHeadmodel, ['elec' SubjStr ]), 'elec');
end


if 1
    for subj=1:length(SubjList)
        load (fullfile(PathScans,['elec' SubjList{subj}])); % loaded is elec
        save(fullfile(PathElec, SubjList{subj}), 'elec');
    end
end