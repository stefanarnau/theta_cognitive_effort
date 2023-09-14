%% Stats on source level
% this is the same than BStats, only that script is less flexible / customized for Bocotilt paper

% Stat part based on tutorial: https://macquarie-meg-research.github.io/MQ_MEG_Scripts/docs/mq_source_statistics_tutorial.html


clear all 
close all
clc

PlotSgSubj=0; 
% ProjectName='Auvibo';  
%    SubjList={'005','006','007','008','009','010','011','012','013','014','015','016','019','020','021','022','023','024','025','027','028','029','030','031','032','033','034','042','049','050'};   %32 subjects without 26
    
ProjectName='Bocotilt';  
  %SubjList={'005','006','007','008','009','011'}; % '010' no dipoli headmodel
   SubjList={'010','011','013','014','015','016','017','018','019','020','021','022','023','024','025','026','027','028','029','030','031','032','034'}; % '008','033'(nodipoli headmodel, bad scan)
   %SubjList={'008','010','011','013','014','015','016','017','018','019','020','021','023','024','025','026','027','028','029','030','031','032','033','034'}; % '022',


if isunix, PathMain='/media/sf_nathaExtern'; else PathMain = 'I:'; end  % alternatively ispc for windows
    PathEeglab = fullfile(PathMain,'Toolboxes','eeglab2021.1'); addpath(PathEeglab); eeglab; 
    %PathFT = fullfile(PathMain,'Toolboxes','fieldtrip-20221126'); addpath(PathFT); ft_defaults;
    PathFT = fullfile(PathMain,'Toolboxes','fieldtrip-20210413'); addpath(PathFT); ft_defaults; % attention: plot on one surface only implemented for fieldtrip-20210413
    PathLib1 = fullfile(PathMain,'Toolboxes','FunctionLib'); addpath(genpath(PathLib1));
AnalysisName='DI500PR010'; % analysis code of outcome of this script, if ='C004P103'; then P001 is loaded
    PathPreProcMat = fullfile(PathMain,[ProjectName AnalysisName(6:end)],'OutputMats');
    PathPreProc = fullfile(PathMain,[ProjectName AnalysisName(6:end)],'3 preprocessed'); 
    PathTrialInfo = fullfile(PathMain,[ProjectName AnalysisName(6:end)], 'Trialinfo' );
    PathSpatFiltMat = fullfile(PathMain,[ProjectName AnalysisName], 'All4' ); mkdir(PathSpatFiltMat);

% Choose one of the 4:
%FolderName='BonStanFreq120FreqRange20Time1300to1300Time2300to1300AllTime-800to2200V1'; AnalysisNamePlot='BonMinStanAlpha';
%FolderName='BonStanFreq60FreqRange20Time10to600Time20to600AllTime-800to2200V1'; AnalysisNamePlot='BonMinStanTheta';
FolderName='SwitchRepeatFreq90FreqRange20Time1700to1200Time2700to1200AllTime-800to2200V1'; AnalysisNamePlot='RepeatMinSwitchAlpha';
%FolderName='SwitchRepeatFreq55FreqRange20Time11500to2200Time21500to2200AllTime-800to2200V1'; AnalysisNamePlot='RepeatMinSwitchTheta';



%% Create input and exchange pos field
% pos field can be exchanged as a mni-wrapped grid model was used 

% load Template to get correct pos field. Must be the same templae used for crteating your sourcemodel: 
Template=fullfile(PathFT,'template','sourcemodel','standard_sourcemodel3d5mm'); 
load(Template); % loaded is sourcemodel


CellSourceCond1=cell(length(SubjList),1);
CellSourceCond2=cell(length(SubjList),1);
for subj=1:length(SubjList)

    load(fullfile(PathSpatFiltMat, FolderName, SubjList{subj})); % loaded is sourceAll, sourceCond1, sourceCond2
%     load(fullfile(PathMriAligned,SubjList{subj})) % loaded is mri
%     mriInd=mri;     

    sourceCond1.pos=sourcemodel.pos; 
    CellSourceCond1{subj}=sourceCond1;

    sourceCond2.pos=sourcemodel.pos; 
    CellSourceCond2{subj}=sourceCond2;

end


%% Stats
cfg                     = [];
cfg.dim                 = CellSourceCond2{1}.dim;
cfg.method              = 'montecarlo';
cfg.statistic           = 'ft_statfun_depsamplesT';
cfg.parameter           = 'pow';
cfg.correctm            = 'cluster';
cfg.computecritval      = 'yes';
cfg.numrandomization    = 4000;   % tut 4000
%cfg.tail                = 1;    % 1=One sided testing
cfg.clusteralpha = 0.05;
cfg.alpha = 0.05; % significance level

% Design Matrix
nsubj                   = numel(CellSourceCond2);
cfg.design(1,:)         = [1:nsubj 1:nsubj];
cfg.design(2,:)         = [ones(1,nsubj) ones(1,nsubj)*2];
% row of design matrix that contains unit variable (in this case: subjects)
cfg.uvar                = 1;
% row of design matrix that contains independent variable (the conditions)
cfg.ivar                = 2; 

% Perform statistical analysis
[stat]                  = ft_sourcestatistics(cfg,CellSourceCond1{:},CellSourceCond2{:});
% save stat
save(fullfile('I:\BocotiltPosterPaper',['Stat' AnalysisNamePlot '.mat']),'stat');





% added NL as if there is no negative (or positive, but this is seldom) cluster, then FT sometimes does not add the fields. This produces
% errors in ft_interpolate at line 129 (fixpos(functional)). This is fixed here:
if or( or (~isfield(stat.negclusters,'prob'),~isfield(stat.negclusters,'clusterstat')),or(~isfield(stat.negclusters,'stddev'),~isfield(stat.negclusters,'cirange') ) )
    statOld=stat; 
    stat.negclusters=struct('prob',[], 'clusterstat',[],'stddev',[],'cirange',[]);
end
if or( or (~isfield(stat.posclusters,'prob'),~isfield(stat.posclusters,'clusterstat')),or(~isfield(stat.posclusters,'stddev'),~isfield(stat.posclusters,'cirange') ) )
    statOld=stat; 
    stat.posclusters=struct('prob',[], 'clusterstat',[],'stddev',[],'cirange',[]);
end



%% Interpolate to template MRI and plot
% from stat tut: Interpolate onto SPM T1 Brain
mri                 = ft_read_mri([PathFT '\template\anatomy\single_subj_T1.nii']);
cfg                 = [];
cfg.voxelcoord      = 'no';
cfg.parameter       = 'stat';
cfg.interpmethod    = 'nearest';
statint             = ft_sourceinterpolate(cfg, stat, mri);


cfg                = [];
cfg.method         = 'surface';
cfg.funcolormap = 'jet';
%cfg.maskparameter  = 'mask';
cfg.funparameter   = 'stat';
cfg.projmethod     = 'nearest';
cfg.surfinflated   = 'surface_white_left.mat'; % surface_inflated_both_caret, surface_pial_both, surface_white_both
cfg.camlight       = 'no';
ft_sourceplot(cfg, statint);
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
light ('Position',[0 0 50])
light ('Position',[0 -50 0])
material dull;
drawnow;
view([100 0]);  % from face: view([0 0]); from the side: view([100 1]) or view([100 0]) or view([1000 0]);
ax = gca; exportgraphics(ax,fullfile('I:\BocotiltPosterPaper',[AnalysisNamePlot '1'  '.tif']),'Resolution',600); % 3:= 100, 4 = 1000 (see view) 

cfg                = [];
cfg.method         = 'surface';
cfg.funcolormap = 'jet';
%cfg.maskparameter  = 'mask';
cfg.funparameter   = 'stat';
cfg.projmethod     = 'nearest';
cfg.surfinflated   = 'surface_white_left.mat'; % surface_inflated_both_caret, surface_pial_both, surface_white_both
cfg.camlight       = 'no';
ft_sourceplot(cfg, statint);
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
light ('Position',[0 0 50])
light ('Position',[0 -50 0])
material dull;
drawnow;
view([1000 0]);  % from face: view([0 0]); from the side: view([100 1]) or view([100 0]) or view([1000 0]);
ax = gca; exportgraphics(ax,fullfile('I:\BocotiltPosterPaper',[AnalysisNamePlot '2'  '.tif']),'Resolution',600); % 3:= 100, 4 = 1000 (see view) 

cfg                = [];
cfg.method         = 'surface';
cfg.funcolormap = 'jet';
%cfg.maskparameter  = 'mask';
cfg.funparameter   = 'stat';
cfg.projmethod     = 'nearest';
cfg.surfinflated   = 'surface_white_right.mat'; % surface_inflated_both_caret, surface_pial_both, surface_white_both
cfg.camlight       = 'no';
ft_sourceplot(cfg, statint);
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
light ('Position',[0 0 50])
light ('Position',[0 -50 0])
material dull;
drawnow;
view([100 0]);  % from face: view([0 0]); from the side: view([100 1]) or view([100 0]) or view([1000 0]);
ax = gca; exportgraphics(ax,fullfile('I:\BocotiltPosterPaper',[AnalysisNamePlot '3'  '.tif']),'Resolution',600); % 3:= 100, 4 = 1000 (see view) 


cfg                = [];
cfg.method         = 'surface';
cfg.funcolormap = 'jet';
%cfg.maskparameter  = 'mask';
cfg.funparameter   = 'stat';
cfg.projmethod     = 'nearest';
cfg.surfinflated   = 'surface_white_right.mat'; % surface_inflated_both_caret, surface_pial_both, surface_white_both
cfg.camlight       = 'no';
ft_sourceplot(cfg, statint);
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
light ('Position',[0 0 50])
light ('Position',[0 -50 0])
material dull;
drawnow;
view([1000 0]);  % from face: view([0 0]); from the side: view([100 1]) or view([100 0]) or view([1000 0]);
ax = gca; exportgraphics(ax,fullfile('I:\BocotiltPosterPaper',[AnalysisNamePlot '4'  '.tif']),'Resolution',600); % 3:= 100, 4 = 1000 (see view) 



%% Significant data - Interpolate and plot
%Max=0.5;
%Calculate adjusted partial eta squared and replace t values with it: 
stat2=stat; 
%stat.stat(stat.stat>0)=0; 
adjpetasq=nan(length(stat.stat),1);   % channel x freq x time
for chan=1:length(stat.stat)
    petasq = (squeeze(stat.stat(chan, :, :)) .^ 2) ./ ((squeeze(stat.stat(chan, :, :)) .^ 2) + (numel(SubjList) - 1));
    adj_petasq = petasq - (1 - petasq) .* (1 / (numel(SubjList) - 1));
    adjpetasq(chan, :, :) = adj_petasq;
end
stat.stat=adjpetasq;
%stat=stat2; 

% calculate max for colorbar: 
Max=nanmax( abs(adjpetasq(:)) );


% Interpolate to template MRI and plot
% from stat tut: Interpolate onto SPM T1 Brain
mri                 = ft_read_mri([PathFT '\template\anatomy\single_subj_T1.nii']);
cfg                 = [];
cfg.voxelcoord      = 'no';
cfg.parameter       = 'stat';
cfg.interpmethod    = 'nearest';
statint             = ft_sourceinterpolate(cfg, stat, mri);
%only decomment if you have a significant cluster:
cfg.parameter       = 'mask';
maskint             = ft_sourceinterpolate(cfg, stat,mri);
statint.mask        = maskint.mask;
%just for debugging: save(fullfile('H:\BocotiltSF101PR010\OutputMats',['StatFromFreqAnalysis']),'stat','statint'); 

%left hemisphere, view 100 0
cfg                = [];
cfg.method         = 'surface';
cfg.funcolormap = 'jet';
cfg.maskparameter  = 'mask';
cfg.funparameter   = 'stat';
cfg.funcolorlim    = [0 Max];
cfg.projmethod     = 'nearest';
cfg.surfinflated   = 'surface_white_left.mat';
cfg.camlight       = 'no';
ft_sourceplot(cfg, statint);
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
%colormap(flipud(brewermap(64,'RdGy'))) % change the colormap, reversed colormap
colormap(brewermap(64,'Greens')) % change the colormap
light ('Position',[0 0 50])
light ('Position',[0 -50 0])
material dull;
drawnow;
view([100 0]);
ax = gca; exportgraphics(ax,fullfile('I:\BocotiltPosterPaper',[AnalysisNamePlot 'Sign' '1'  '.tif']),'Resolution',600); % 3:= 100, 4 = 1000 (see view) 

% left hemisphere, view 1000 0
cfg                = [];
cfg.method         = 'surface';
cfg.funcolormap = 'jet';
cfg.maskparameter  = 'mask';
cfg.funparameter   = 'stat';
cfg.funcolorlim    = [0 Max];
cfg.projmethod     = 'nearest';
cfg.surfinflated   = 'surface_white_left.mat';
cfg.camlight       = 'no';
ft_sourceplot(cfg, statint);
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
%colormap(flipud(brewermap(64,'RdGy'))) % change the colormap, reversed colormap
colormap(brewermap(64,'Greens')) % change the colormap
light ('Position',[0 0 50])
light ('Position',[0 -50 0])
material dull;
drawnow;
view([1000 0]);
ax = gca; exportgraphics(ax,fullfile('I:\BocotiltPosterPaper',[AnalysisNamePlot 'Sign' '2'  '.tif']),'Resolution',600); % 3:= 100, 4 = 1000 (see view) 

% right hemisphere, view 100 0
cfg                = [];
cfg.method         = 'surface';
cfg.funcolormap = 'jet';
cfg.maskparameter  = 'mask';
cfg.funparameter   = 'stat';
cfg.funcolorlim    = [0 Max];
cfg.projmethod     = 'nearest';
cfg.surfinflated   = 'surface_white_right.mat';
cfg.camlight       = 'no';
ft_sourceplot(cfg, statint);
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
%colormap(flipud(brewermap(64,'RdGy'))) % change the colormap, reversed colormap
colormap(brewermap(64,'Greens')) % change the colormap
light ('Position',[0 0 50])
light ('Position',[0 -50 0])
material dull;
drawnow;
view([100 0]);
ax = gca; exportgraphics(ax,fullfile('I:\BocotiltPosterPaper',[AnalysisNamePlot 'Sign' '3'  '.tif']),'Resolution',600); % 3:= 100, 4 = 1000 (see view) 

% right hemisphere, view 1000 0
cfg                = [];
cfg.method         = 'surface';
cfg.funcolormap = 'jet';
cfg.maskparameter  = 'mask';
cfg.funparameter   = 'stat';
cfg.funcolorlim    = [0 Max];
cfg.projmethod     = 'nearest';
cfg.surfinflated   = 'surface_white_right.mat';
cfg.camlight       = 'no';
ft_sourceplot(cfg, statint);
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
%colormap(flipud(brewermap(64,'RdGy'))) % change the colormap, reversed colormap
colormap(brewermap(64,'Greens')) % change the colormap
light ('Position',[0 0 50])
light ('Position',[0 -50 0])
material dull;
drawnow;
view([1000 0]);
ax = gca; exportgraphics(ax,fullfile('I:\BocotiltPosterPaper',[AnalysisNamePlot 'Sign' '4'  '.tif']),'Resolution',600); % 3:= 100, 4 = 1000 (see view) 


%% Work with atlas
%atlas = ft_read_atlas([PathFT '\template\atlas\brainnetome\BNA_MPM_thr25_1.25mm.nii']); % brainnetome atlas
atlas = ft_read_atlas([PathFT '\template\atlas\aal\ROI_MNI_V4.nii']); % aal

% Interpolate atlas on anatomical data (if you skip that step, hemispheres of aal atlas are exchanged and other atlas do not work anyway 
cfg = [];
cfg.interpmethod = 'nearest';
cfg.parameter = 'tissue';
atlasNew = ft_sourceinterpolate(cfg, atlas, mri);
atlas=atlasNew;

TissueReshaped=reshape(atlas.tissue,[902629 1]);
mask=statint.mask; 
mask=mask==1; % convert double to logical
TissueReshapedMasked=TissueReshaped; 
TissueReshapedMasked(~mask)=nan; 

ParcelNumsFrequencies=nan(length(atlas.tissuelabel),1); 
for parcel=1:length(atlas.tissuelabel)
    ParcelNumsFrequencies(parcel)=length(find(TissueReshapedMasked==parcel)); 
    %SumOfStatValue=
end
ParcelNumsFrequenciesCropped=ParcelNumsFrequencies(ParcelNumsFrequencies>0) 
Labels=atlas.tissuelabel(ParcelNumsFrequencies>0)'


%% Plot data with atlas
statint.coordsys='mni'; % NL: must be mni as data has been warped to fit to template sourcemodel and template sourcemodel is mni
atlas = ft_read_atlas([PathFT '\template\atlas\aal\ROI_MNI_V4.nii']); % aal
cfg                = [];
cfg.method         = 'ortho';
cfg.funcolormap = 'jet';
cfg.atlas          = atlas; 
cfg.funparameter   = 'stat';
cfg.maskparameter   = 'mask';
cfg.projmethod     = 'nearest';
ft_sourceplot(cfg, statint); 


%% Loop up atlas parcels with ft_volumelookup
% atlas = ft_read_atlas([PathFT '\template\atlas\aal\ROI_MNI_V4.nii']); % aal
% statint.coordsys='mni';
% cfg=[]; 
% cfg.atlas  =atlas;             %= string, filename of atlas to use, see FT_READ_ATLAS
% cfg.maskparameter='mask';        %= string, field in volume to be looked up, data in field should be logical
% cfg.minqueryrange  =1;    % = number, should be odd and <= to maxqueryrange (default = 1)
% cfg.maxqueryrange  =9;    % = number, should be odd and >= to minqueryrange (default = 1)
% OutputVolLookup = ft_volumelookup(cfg, statint); 