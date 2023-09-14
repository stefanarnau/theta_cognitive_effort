%% MRI preprocessing - create sourcemodel
% calculates regular grid wrapped in mni template sourcemodel,
% cp. tutorial: https://www.fieldtriptoolbox.org/tutorial/sourcemodel/

% by NL (nathalie.liegel@gmail.com), last change: 10.09.22

clear all
close all
clc


%% Parameters
% ProjectName='Auvibo';  % ='Auvibo';
% SubjList={'005','006','007','008','009','012','013','014','015','016','019','020','021','022','024','025','027','029','030','031','032','033','034','035','036','049','050'}; %
% SubjList={'010','011','023','028','042'};
ProjectName='Bocotilt'; 
SubjList={'008','010','011','013','014','015','016','017','018','019','020','021','022','023','024','025','026','027','028','029','030','031','032','033','034'}; 
%SubjList={'022'};

% Path defintions. SourcemodelNum refers to my path system 
% (which is modularly organized: SMyyyALxxx=Alignment of scan and mri was done, using version xxx, the sm=sourcemodel with parameters referring to version yyy). 
% You might just add your pathes here: 
SourcemodelNum='SM500AL002'; %template resolution=5 mm --> 500, template reso.: 10mm--> 100, template reso.=7,5 mm--> 750
PathMain='I:';  
    PathFT = fullfile(PathMain, 'Toolboxes','fieldtrip-20221126'); addpath(PathFT);ft_defaults;
    PathMriAligned=fullfile(PathMain,[ProjectName 'DataBase'],SourcemodelNum(6:end));
    PathSourcemodel=fullfile(PathMain,[ProjectName 'DataBase'],SourcemodelNum); mkdir(PathSourcemodel);
% Fieldtrip template used to normalize the regular grid. The resolution (here 5mm) of this templates defines the resolution of your grid: 
% Find more templates (e.g. 4mm, 7point5 mm, 7mm, 10mm) in in your local FT.
% For virtual channel analysis, Robert O. recommends a resolution of 8 mm or better
Template=fullfile(PathFT,'template','sourcemodel','standard_sourcemodel3d5mm'); %standard_sourcemodel3d7point5mm
    

%% 
for subj=1:length(SubjList)
    

    % Create sourcemodel
    load(fullfile(PathMriAligned,SubjList{subj})) % loaded is mri
    load(Template);

    % cp. in sourcemodel tut: section=subject-specific-grids-that-are-equivalent-across-subjects-in-normalized-space:
    template_grid = sourcemodel;
    template_grid=ft_convert_units(template_grid,'mm');

    clear sourcemodel;
    cfg           = [];
    cfg.warpmni   = 'yes';
    cfg.template  = template_grid;
    cfg.nonlinear = 'yes';
    cfg.mri       = mri;
    cfg.unit      ='mm';
    sourcemodel        = ft_prepare_sourcemodel(cfg);

    
    % Save
    save (fullfile(PathSourcemodel,SubjList{subj}),'sourcemodel'); 
end

