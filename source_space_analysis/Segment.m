%% MRI preprocessing - segment T1 scan 
% cp. tutorial: https://www.fieldtriptoolbox.org/tutorial/headmodel_eeg_bem/

% Segment the brain in scalp, skin and brain. You need that for the BEM headmodel (widely used for EEG source analysis)
% Segmentation is a very crucial step to get a good headmodel
% Even if it runs smoothly, the segmentation might not be very good. 
% Plot it using the script "CreateMesh" to see how good your segmentaiton is!
% I had problems with the scalp segementation and therefore run the script for several different scalpsmooth and scalpthreshold values. 
% You might want to loop over other parameters or delete the loop

% by NL (nathalie.liegel@gmail.com), last change: 10.09.22


%% Params
% ProjectName='Auvibo';  
%     SubjList={'005','006','007','008','009','010','011','012','013','014','015','016','019','020','021','022','024','025','027','029','030','031','032','033','034','035','036','049','050'}; %'023','026','028',
%     SubjList={'023','028','042'};
ProjectName='Bocotilt'; % '010','011','013','014','015','016','017'
    SubjList={'008','010','011','013','014','015','016','017','018','019','020','021','022','023','024','025','026','027','028','029','030','031','032','033','034'};
   %SubjList={'031'};  % for single subject plotting and param determination

% Specify PathFT and PathMriAligned: You might want to delete the variable AnalysisSteps and PathMain and just copy-paste your pathes!
% (AnalysisSteps refers to my modular path system: first alignement of mri and scan (ALxxx) was done using parametrs of version xxx, then segmentation (SG) ) :
AnalysisSteps='SG001AL002';  
% different PathMain depending on which system script is run:
if isunix, PathMain='/media/sf_nathaExtern'; else, PathMain = 'I:'; end  % alternatively ispc for windows
    PathFT = fullfile(PathMain, 'Toolboxes','fieldtrip-20221126'); addpath(PathFT); ft_defaults; % fieldtrip-20221126 (no changes made by NL), standard: fieldtrip-20220104 (change sin ft_volumesegment, made by NL)
    PathMriAligned=fullfile(PathMain,[ProjectName 'DataBase'],AnalysisSteps(6:end)); 


% Loop through different parameters (and save each variant) to get an idea which scalpsmooth works:

for scalpsmooth=10  % 5:20  =[10 8 5]
for scalpthreshold=[0.07]  %=[0.05 0.07 0.1 0.3 0.5]

    % where to save data:
    PathMriSeg=fullfile(PathMain,[ProjectName 'DataBase'],AnalysisSteps,['thre' num2str(scalpthreshold)],['smooth' num2str(scalpsmooth)]); mkdir(PathMriSeg);


    %% Params end
        
    for subj=1:length(SubjList)

    % Load
    load(fullfile(PathMriAligned,SubjList{subj}))


    cfg           = [];
    cfg.output    = {'brain','skull','scalp'};
    cfg.scalpthreshold = scalpthreshold; % 'no', or scalar, relative threshold value which is used to threshold the anatomical data in order to create a volumetric scalpmask, default = 0.1)       
    cfg.scalpsmooth  = scalpsmooth;  % = 'no', or scalar, the FWHM of the gaussian kernel in voxels (default = 5)
        % only optional: 
        %cfg.spmmethod ='new';      %= string, 'old', 'new', 'mars'
        %cfg.skullsmooth=0.5; % default: 0.5
%     segmentedmri  = ft_volumesegmentNl(cfg, mri,ThresholdZ);  % return serror: Undefined function 'align_ijk2xyz' for input arguments of type 'struct'.
%                     Error in ft_volumesegmentNl (line 309)
%                         [mri, permutevec, flipflags] = align_ijk2xyz(mri); 
    segmentedmri  = ft_volumesegment(cfg, mri);

    % Save
    save(fullfile(PathMriSeg,SubjList{subj}),'segmentedmri'); 


    end
end
end
