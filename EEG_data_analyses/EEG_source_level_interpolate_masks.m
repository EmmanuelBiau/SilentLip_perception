%% source grandaverage across participants;
clearvars;clc; close all;addpath XXX\fieldtrip; 
cd XXX\;
MI_source = dir('allsubj_MI_source_peak_VISUAL_n*.mat');
load(MI_source.name);

%% Interpolate source data with mask to inspect the left audio ROI, left visual ROI and right visual ROI;
close all;
close all; clear MI_diff_mina allsource_MI_mina no_mask_int_T2vsT1;

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

%Load mri template from FT and the grid template from one participant (e.g. subj2);
cd XXX\fieldtrip\template\headmodel;
load standard_mri;
cd XXX\headmodel_mat\subj2\;
load subj2_grid;

%Using template grid;
sourceTmpl = [];
sourceTmpl.inside = grid.inside;
sourceTmpl.dim = grid.dim;
sourceTmpl.pos = grid.pos;
sourceTmpl.unit = grid.unit;

%Apply this source 'template' to sound and movie data and mask the data;
%mask the left hemisphere (auditory source detection;
xa = 0;
yl = 7;
zl = 8; 
left_audio_T2vsT1 = sourceTmpl;
left_audio_T2vsT1.powspctrm = nan(size(grid.pos,1),1);
left_audio_T2vsT1.powspctrm(sourceTmpl.inside & grid.pos(:,1) <= xa & grid.pos(:,2) <= yl & grid.pos(:,3) <= zl) = GA_source_mina.powspctrm(grid.pos(grid.inside,1) <= xa & grid.pos(grid.inside,2) <= yl & grid.pos(grid.inside,3) <= zl);

%mask the left visual hemisphere;
xl = 0;
yl = -7;
zl = 8; 
left_visual_T2vsT1 = sourceTmpl;
left_visual_T2vsT1.powspctrm = nan(size(grid.pos,1),1);
left_visual_T2vsT1.powspctrm(sourceTmpl.inside & grid.pos(:,1) <= xl & grid.pos(:,2) <= yl & grid.pos(:,3) <= zl) = GA_source_mina.powspctrm(grid.pos(grid.inside,1) <= xl & grid.pos(grid.inside,2) <= yl & grid.pos(grid.inside,3) <= zl);

%mask the right visual hemisphere;
xr = 0;
yr = -7;
zr = 8; 
right_visual_T2vsT1 = sourceTmpl;
right_visual_T2vsT1.powspctrm = nan(size(grid.pos,1),1);
right_visual_T2vsT1.powspctrm(sourceTmpl.inside & grid.pos(:,1) >= xr & grid.pos(:,2) <= yr & grid.pos(:,3) <= zr) = GA_source_mina.powspctrm(grid.pos(grid.inside,1) >= xr & grid.pos(grid.inside,2) <= yr & grid.pos(grid.inside,3) <= zr);

%no-MASK(to determine also a control source with no activity without masking);
no_mask_T2vsT1 = sourceTmpl;
no_mask_T2vsT1.powspctrm = nan(size(grid.pos,1),1); 
no_mask_T2vsT1.powspctrm(sourceTmpl.inside)= GA_source_mina.powspctrm;

%Interpolate the parameter 'pow';
cfg = [];
cfg.downsample = 2; 
cfg.parameter = 'powspctrm';
left_audio_int_T2vsT1 = ft_sourceinterpolate(cfg, left_audio_T2vsT1, mri); 
left_visual_int_T2vsT1 = ft_sourceinterpolate(cfg, left_visual_T2vsT1, mri); 
right_visual_int_T2vsT1 = ft_sourceinterpolate(cfg, right_visual_T2vsT1, mri); 
no_mask_int_T2vsT1 = ft_sourceinterpolate(cfg, no_mask_T2vsT1, mri); 

%Plot source localizations in masked regions of interest and find the greatest voxel manually (without atlas);
close all;
atlas = ft_read_atlas('XXX\fieldtrip\template\atlas\aal\ROI_MNI_V4.nii');
cfg = [];
cfg.method = 'slice';                   
cfg.funparameter = 'powspctrm';         
cfg.funcolorlim = 'maxabs';             
cfg.maskparameter = 'powspctrm';        
cfg.opacitylim =  'zeromax'; 
cfg.funcolormap = 'jet'; 
cfg.atlas = atlas;
cfg.location = 'max';
ft_sourceplot(cfg, left_audio_int_T2vsT1);
set(gcf,'name',['left_audio MASK-N', num2str(length(subjects))]);
ft_sourceplot(cfg, left_visual_int_T2vsT1);
set(gcf,'name',['left_VISUAL MASK-N', num2str(length(subjects))]);
ft_sourceplot(cfg, right_visual_int_T2vsT1);
set(gcf,'name',['right_VISUAL MASK-N', num2str(length(subjects))]);   
ft_sourceplot(cfg, no_mask_int_T2vsT1);
set(gcf,'name',['No-MASK-N', num2str(length(subjects))]);      
       
%% Now look at the masked data with the atlas and detect max voxel automatically and corresponding rid coordinate in each region for later source reconstruction; 
    
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
ft_sourceplot(cfg, left_audio_T2vsT1);
set(gcf,'name',['left_AUDIO MASK-N', num2str(length(subjects))]);
ft_sourceplot(cfg, left_visual_T2vsT1);
set(gcf,'name',['left_VISUAL MASK-N', num2str(length(subjects))]);
ft_sourceplot(cfg, right_visual_T2vsT1);
set(gcf,'name',['right_VISUAL MASK-N', num2str(length(subjects))]);
ft_sourceplot(cfg, no_mask_T2vsT1);
set(gcf,'name',['right_VISUAL MASK-N', num2str(length(subjects))]);    

%Get the left audio voxels activity and their indexes, ordered;
[left_audio_source_values, left_audio_source_idx] = sort(left_audio_T2vsT1.powspctrm,'descend'); %value = MI; idx = give the #grid corresponding to that MI value; Then go in grid.pos and find the coordinates of the source of interest;  
%Get the greatest source activity and its index;
left_audio_source_max_value = max(left_audio_source_values);
left_audio_source_max_idx = left_audio_source_idx(left_audio_source_values == max(left_audio_source_values));
%Get the grid coordinates corresponding to the the source of interest;
left_audio_source_grid = grid.pos(left_audio_source_max_idx,1:3);

%Get the left visual voxels activity and their indexes, ordered;
[left_visual_source_values, left_visual_source_idx] = sort(left_visual_T2vsT1.powspctrm,'descend');  
%Get the greatest source activity and its index;
left_visual_source_max_value = max(left_visual_source_values);
left_visual_source_max_idx = left_visual_source_idx(left_visual_source_values == max(left_visual_source_values));
%Get the grid coordinates corresponding to the the source of interest;
left_visual_source_grid = grid.pos(left_visual_source_max_idx,1:3);

%Get the right visual voxels activity and their indexes, ordered;
[right_visual_source_values, right_visual_source_idx] = sort(right_visual_T2vsT1.powspctrm,'descend'); 
%Get the greatest source activity and its index;
right_visual_source_max_value = max(right_visual_source_values);
right_visual_source_max_idx = right_visual_source_idx(right_visual_source_values == max(right_visual_source_values));
%Get the grid coordinates corresponding to the the source of interest;
right_visual_source_grid = grid.pos(right_visual_source_max_idx,1:3);

%-------------------------------------------------------------------------------------------------------------------------%
%Get the control source voxel and index (in no-mask data);
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
set(gcf,'name',['right_VISUAL MASK-N', num2str(length(subjects))]);

control_value = -0.241659; % <--- ENTER THE VALUE OF INTEREST;

[source_values, source_idx] = sort(no_mask_T2vsT1.powspctrm,'descend'); 
[~,control_source_value] = min(pdist2(control_value, source_values));
control_source_idx = source_idx(control_source_value);
%Get the grid coordinates corresponding to the the source of interest;
control_source_grid = grid.pos(control_source_idx,1:3);
%-------------------------------------------------------------------------------------------------------------------------%

%Now check the position of the selected ROI grids on the head model;
close all;
cd XXX\headmodel_mat\subj2\;
load subj2_grid;
cd XXX\headmodel_mat\subj2\;
load subj2_vol;

source_grids = [left_audio_source_max_idx; left_visual_source_max_idx; right_visual_source_max_idx; control_source_idx];
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

%design final matrix and export thier coordinates;
ROIs_coordinates = cell(length(source_grids), 4);
ROIs_coordinates = array2table(ROIs_coordinates, 'Variable',{'ROI','VOXEL_VALUE','GRID_NUMBER','GRID_POS'});
ROIs_coordinates.ROI{1} = 'left_audio_source';
ROIs_coordinates.ROI{2} = 'left_visual_source';
ROIs_coordinates.ROI{3} = 'right_visual_source';
ROIs_coordinates.ROI{4} = 'control_source';    
ROIs_coordinates.VOXEL_VALUE{1} = left_audio_source_max_value;
ROIs_coordinates.VOXEL_VALUE{2} = left_visual_source_max_value;
ROIs_coordinates.VOXEL_VALUE{3} = right_visual_source_max_value;  
ROIs_coordinates.VOXEL_VALUE{4} = control_value;   
ROIs_coordinates.GRID_NUMBER{1} = left_audio_source_max_idx;
ROIs_coordinates.GRID_NUMBER{2} = left_visual_source_max_idx;
ROIs_coordinates.GRID_NUMBER{3} = right_visual_source_max_idx;  
ROIs_coordinates.GRID_NUMBER{4} = control_source_value;     
ROIs_coordinates.GRID_POS{1} = left_audio_source_grid;
ROIs_coordinates.GRID_POS{2} = left_visual_source_grid;
ROIs_coordinates.GRID_POS{3} = right_visual_source_grid;   
ROIs_coordinates.GRID_POS{4} = control_source_grid;
cd XXX\;
save ROIs_coordinates ROIs_coordinates;
   
%% 