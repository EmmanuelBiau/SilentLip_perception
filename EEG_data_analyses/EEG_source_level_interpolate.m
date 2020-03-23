%% source grandaverage across participants
clearvars;clc; addpath XXX\fieldtrip; cd XXX\;
MI_source = dir('allsubj_MI_source_peak_VISUAL_n*.mat');
load(MI_source.name);

%% Interpolate MI source data;
close all; clear MI_diff_mina allsource_MI_mina no_mask_int_T2vsT1;

%subjects ID;
subjects = [2:9 11:17 19:26];

%Create a structure with all the MI_source data together;
allsource_MI_mina = cell(1, length(subjects));
allsource_MI_mina{1, length(subjects)} = [];

scount = 0;
for s = subjects
    
    MI_diff_mina = allsubj_MI_source_mina_peak_post{s};
    MI_diff_mina.powspctrm = [];
    MI_diff_mina.powspctrm = (allsubj_MI_source_mina_peak_post{s}.powspctrm(:,3) - allsubj_MI_source_mina_peak_pre{s}.powspctrm(:,3))./allsubj_MI_source_mina_peak_pre{s}.powspctrm(:,3);
    MI_diff_mina.freq = 3;
    
    scount = scount+1;
    allsource_MI_mina{scount} = MI_diff_mina;
    
    clear MI_diff_mina_VISUAL
    
end

%Grand average the MI_source;
cfg = [];
cfg.parameter = 'powspctrm';
GA_source_mina = ft_freqgrandaverage(cfg,allsource_MI_mina{:});

% Load mri template from FT and the grid template from one participant (e.g. subj2);
cd XXX\fieldtrip\template\headmodel;
load standard_mri;
cd XXX\headmodel_mat\subj2\;
load subj2_grid;

% Using template grid;
sourceTmpl = [];
sourceTmpl.inside = grid.inside;
sourceTmpl.dim = grid.dim;
sourceTmpl.pos = grid.pos;
sourceTmpl.unit = grid.unit;

% Apply this source 'template' to sound and movie data and mask the data;
%No-MASK;
no_mask_T2vsT1 = sourceTmpl;
no_mask_T2vsT1.powspctrm = nan(size(grid.pos,1),1); 
no_mask_T2vsT1.powspctrm(sourceTmpl.inside)= GA_source_mina.powspctrm;

%Interpolate the parameter 'pow';
cfg = [];
cfg.downsample = 2; % downsample the MRI resolution by 2 (i.e. half the resolution)
cfg.parameter = 'powspctrm';
%sounds;
no_mask_int_T2vsT1 = ft_sourceinterpolate(cfg, no_mask_T2vsT1, mri); 

%Plot source localizations in movies and sounds conditions;
close all;
atlas = ft_read_atlas('XXX\fieldtrip\template\atlas\aal\ROI_MNI_V4.nii');
cfg = [];
cfg.method = 'ortho';                   
cfg.funparameter = 'powspctrm';         
cfg.funcolorlim = 'maxabs';             
cfg.maskparameter = 'powspctrm';        
cfg.opacitylim =  'zeromax'; 
cfg.funcolormap = 'jet'; 
cfg.atlas = atlas;
cfg.location = 'max';
cfg.crosshair = 'no';
ft_sourceplot(cfg, no_mask_int_T2vsT1);
set(gcf,'name',['no_mask_N', num2str(length(subjects))]);


%% Localize voxels of interest (exploratory method on non-masked data);

%Plot source localizations in movies and sounds conditions;
close all;
cfg = [];
cfg.method = 'ortho';                  
cfg.funparameter = 'powspctrm';       
cfg.funcolorlim = 'maxabs';            
cfg.maskparameter = 'powspctrm';      
cfg.opacitylim =  'zeromax'; 
cfg.funcolormap = 'jet';
cfg.location = 'max';
ft_sourceplot(cfg, no_mask_T2vsT1);
set(gcf,'name',['no_mask_N', num2str(length(subjects))]);  


value_of_interest = -0.072939; % <--- ENTER THE VALUE OF INTEREST;

[source_values, source_idx] = sort(no_mask_T2vsT1.powspctrm,'descend'); %value = MI; idx = give the #grid corresponding to that MI value; Then go in grid.pos and find the coordinates of the source of interest;  
[~,value_of_interest_idx] = min(pdist2(value_of_interest, source_values));
value_of_interest_idx_grid = source_idx(value_of_interest_idx);
%Get the grid coordinates corresponding to the the source of interest;
value_of_interest_grid = grid.pos(value_of_interest_idx_grid,1:3);

%Check the position of the grid on the head model;
close all;
cd XXX\headmodel_mat\subj2\;
load subj2_grid;
cd XXX\headmodel_mat\subj2\;
load subj2_vol

source_grids = [value_of_interest_idx_grid];

grid.pos = grid.pos(source_grids,:);
%plot only the grid positions within the brain.
figure;
ft_plot_mesh(grid.pos)
% plot the BEM
ft_plot_mesh(vol.bnd(1), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
hold on;
ft_plot_mesh(vol.bnd(2),'edgecolor','none','facealpha',0.4);
hold off
cd XXX\headmodel_mat\subj2\;
load subj2_grid;

%%