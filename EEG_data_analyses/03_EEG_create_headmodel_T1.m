% This first part is where we will check the correct orientation of the participant's T1 and align the electrodes to the mri; 
close all;clearvars;clc;
%prepare the paths to the different toolboxes of interest; 
restoredefaultpath
addpath 'XXX\fieldtrip';ft_defaults
addpath(genpath('XXX\headmodel_tools'));addpath(genpath('XXX\fieldtrip\external\freesurfer'));
addpath(genpath('XXX\fieldtrip\external\iso2mesh'));addpath(genpath('XXX\fieldtrip\external\openmeeg'));
% addpath(genpath('C:\toolbox\spm8'));
%-------------------------------------------------------------------------------------------------------------------------------

% 1-Reorient the mri in the conventional RAS orientation (positive values at the Right, Anterior and Superior poles of respective X, Y and Z axis);

%subject ID;
subjects = [3];

%Load the subject's T1 and create the mri for fieldtrip;
s = subjects;
cd (['XXX\subj_',num2str(s)]);
mri = ft_read_mri('xxx.nii');
%save subject's individual mri;
cd XXX\headmodel_mat\
mkdir(['subj',num2str(s)]);
cd(['XXX\headmodel_mat\subj',num2str(s)]);
save(['subj',num2str(s),'_mri'],'mri');

%Check the handedness of the transformation matrix(i.e. check whether the Right on the X-axis is positive or not);
det(mri.transform); % 1:right-handed; -1:left-handed(so needs to be realigned to be right-handed);

%Determine the coordinate-system (check visually the +/- polarity of each axis);
mri = ft_determine_coordsys(mri,'interactive','yes');

%Realign coordsys to spm coordinates (find a,p,s); asks first to enter the polarity of each axis in the original mri;
%Then reorient them manually with the anterior/posterior commisures and right side; 
cfg = [];
cfg.coordsys = 'spm';
cfg.method = 'interactive';
mri_spm = ft_volumerealign(cfg, mri);

%Reslice the anatomical volume in a way that each slice will be equalythick;
cfg = [];
mri_rs = ft_volumereslice(cfg, mri_spm);

%Check the handedness of the tranformation matrix again: Here, this would be always +1. If not, start mri reorientiation again;  
det(mri_rs.transform);

%Only if handedness is correct, then save the new well-oriented mri;
cd(['XXX\headmodel_mat\subj',num2str(s)]);
ft_write_mri(['subj',num2str(s),'_spm.nii'], mri_rs.anatomy, 'dataformat', 'nifti', 'transform', mri_rs.transform);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2- Realign the electrodes on individual mri by using fiducial markers;

%subject ID;
subjects = [3];

%Read the polhemus file of the participant; 
s = subjects;
cd(['XXX\subj_',num2str(s)]);
pos_file = (dir(['subj',num2str(s),'*.pos']));
elec = ft_read_sens(pos_file.name); 
%convert from cm to mm units;
elec = ft_convert_units(elec, 'mm');

%Read the reoriented mri of the participant;
cd(['XXX\headmodel_mat\subj',num2str(s)]);
mri_rs = ft_read_mri(['subj',num2str(s),'_spm.nii']);

%Bring into ctf system (not really we just want the original fid_coordinates);
%Here you have to mark manually the coordinates of the nasion, right and left ear lobes (as in polhemus marking); needs a bit of practice;
cfg = [];
cfg.coordsys = 'ctf';
cfg.method = 'interactive';
mri_ctf = ft_volumerealign(cfg, mri_rs);

%Get the fiducial coordinates in voxelspace of the original MRI;
%fiducials saved in mri structure;
vox_Nas = mri_ctf.cfg.fiducial.nas;  
vox_Lpa = mri_ctf.cfg.fiducial.lpa;     
vox_Rpa = mri_ctf.cfg.fiducial.rpa;

%Use the original transformation Matrix;
Transform = mri_rs.transform;

%Transform the voxel indices of the mri to head coordinates in mm (SPM)
head_Nas = ft_warp_apply(Transform, vox_Nas, 'homogenous'); % nasion 
head_Lpa = ft_warp_apply(Transform, vox_Lpa, 'homogenous'); % Left preauricular
head_Rpa = ft_warp_apply(Transform, vox_Rpa, 'homogenous'); % Right preauricular

% Now, align the electrodes to the mri head positions; 
cfg = [];
cfg.method = 'template';
cfg.elec = elec;
cfg.fiducial = {'Nasion', 'LPA', 'RPA'};
cfg.template.elecpos(1,:) = head_Nas;
cfg.template.elecpos(2,:) = head_Lpa;
cfg.template.elecpos(3,:) = head_Rpa;
cfg.template.label = {'Nasion', 'LPA', 'RPA'};
elec_aligned = ft_electroderealign(cfg);

cd(['XXX\headmodel_mat\subj',num2str(s)]);
save(['subj',num2str(s),'_elec_aligned'],'elec_aligned');

%% Segmentation: 
%This step uses a segmentation of the anatomy, where gray, white and the cerebro-spinal fluid compartments are differentiated, to create a skull-stripped anatomy and a brainmask. 
%The brainmask is a binary mask of the inner skull.

%subject ID;
subjects = [3];

%Load the mri of the participant;
s = subjects;
folder = 'XXX\headmodel_mat\';
mri_file = [folder '\subj' num2str(s) '\subj',num2str(s) '_spm.nii'];
mri_rs = ft_read_mri(mri_file);

addpath(genpath('XXX\spm8')); % We need to call a spm function for this step so we add the path here but need to remove the path after because of path issues FT/spm;
seg = create_bem_segmentation('t1filename',mri_file,'writeseg','yes');
rmpath(genpath('XXX\spm8')); % remove spm path here now; 

%fake knowledge;
seg.transform = mri_rs.transform;
seg.dim = mri_rs.dim;
seg.unit = mri_rs.unit;
seg.coordsys = 'spm';

%Save segmented subject's MRI with the four compartments;
cd(['XXX\headmodel_mat\subj',num2str(s)]);
save(['subj',num2str(s),'_seg'],'seg');

%% Plot masks to check for correctness of compartments segmentation (i.e. scalp, skull,csf and brain);

%subject ID;
subjects = [3];

s = subjects;
cd(['XXX\headmodel_mat\subj',num2str(s)]);
load(['subj',num2str(s),'_seg']);
figure, ft_plot_slice(seg.scalp,'location',[88 128 128],'orientation',[0 0 1]), hold on
ft_plot_slice(seg.skull,'location',[88 128 128],'orientation',[0 1 0]);
ft_plot_slice(seg.csf+seg.brain,'location',[88 128 128],'orientation',[1 0 0]);

%% Generate compartmented meshes;

%Find first instance of skull
for i = 1:size(seg.skull,3)
    
    if any(any(seg.skull(:,:,i) > 0))
        bottom_idx = i;
        break
    end
    
end

%Fill the holes in the bottom slices of the different elements. If there
%are holes in the bottom slices, there will be errors in the mesh constructions.
seg.scalp(:,:,1)=imfill(seg.scalp(:,:,1),'holes'); % fill bottom slice of scalp

seg.skull(:,:,1)=false(size(seg.skull,1),size(seg.skull,2)); % empty bottom slice of skull
seg.skull(:,:,2)=imfill(seg.skull(:,:,2),'holes'); % fill 2nd from bottom slice of skull

seg.csf(:,:,1:2)=false(size(seg.csf,1),size(seg.csf,2),2); % empty bottom slice of skull
seg.csf(:,:,3)=imfill(seg.csf(:,:,3),'holes'); % fill 2nd from bottom slice of CSF;

clear cfg
cfg.tissue = {'scalp', 'skull', 'csf', 'brain'};
cfg.method = 'iso2mesh'; % 'projectmesh';
cfg.numvertices = 10000; % We'll decimate later;
bnd = ft_prepare_mesh(cfg, seg);

%Decimate;
[bnd(1).pos, bnd(1).tri] = meshresample(bnd(1).pos, bnd(1).tri, 1000/size(bnd(1).pos,1));
[bnd(2).pos, bnd(2).tri] = meshresample(bnd(2).pos, bnd(2).tri, 2000/size(bnd(2).pos,1));
[bnd(3).pos, bnd(3).tri] = meshresample(bnd(3).pos, bnd(3).tri, 3000/size(bnd(3).pos,1));
[bnd(4).pos, bnd(4).tri] = meshresample(bnd(4).pos, bnd(4).tri, 3000/size(bnd(4).pos,1));

%Check and repair individual meshes using iso2mesh;
for ii = 1:length(bnd)
    
    [bnd(ii).pos, bnd(ii).tri] = meshcheckrepair(bnd(ii).pos, bnd(ii).tri, 'dup');
    [bnd(ii).pos, bnd(ii).tri] = meshcheckrepair(bnd(ii).pos, bnd(ii).tri, 'isolated');
    [bnd(ii).pos, bnd(ii).tri] = meshcheckrepair(bnd(ii).pos, bnd(ii).tri, 'deep');
    [bnd(ii).pos, bnd(ii).tri] = meshcheckrepair(bnd(ii).pos, bnd(ii).tri, 'meshfix');    
    bnd(ii).pnt = bnd(ii).pos;
 
end

%Ensure there are no overlaps;
bnd = decouplesurf(bnd); %decouplesurf is an unimplemented subfunction temporarily stashed in prepare_mesh_segmentation;

%save the individual BEM model;
cd(['XXX\headmodel_mat\subj',num2str(s)]);
save(['subj',num2str(s),'_bnd'],'bnd');

%% Plot headmodel's compartments together with the electrodes (not perfectely aligned yet);

%subject ID;
subjects = [3];

s = subjects;
cd(['XXX\headmodel_mat\subj',num2str(s)]);
load(['subj',num2str(s),'_bnd']); load(['subj',num2str(s),'_elec_ft']);
figure;
ft_plot_mesh(bnd(1),'facecolor','none'); %scalp
hold;
ft_plot_mesh(bnd(2), 'facecolor',[0.1 0.1 0.1], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05); %skull
hold ;
ft_plot_mesh(bnd(3),'edgecolor','none','facealpha',0.4); %csf
hold;
ft_plot_mesh(bnd(4),'edgecolor','none','facecolor',[0.6 0.4 0.4], 'facealpha', 0.7); %brain
hold;set(gcf, 'color', 'w')
hold;
ft_plot_sens(elec_ft ,'style', 'b*');

%% Realign the electrodes on scap with interactive align the electrodes to the head: Here you manually rotate,translate and scale in X,Y and Z axis;

elec_ft = ft_convert_units(elec_ft,'mm');
cfg = [];
cfg.method = 'interactive';
cfg.elec = elec_ft;
cfg.headshape = vol.bnd(1);
elec_ft = ft_electroderealign(cfg);

%Once correct, save the realigned electrode coordinates; 
cd(['XXX\headmodel_mat\subj',num2str(s)]);
save(['subj',num2str(s),'_elec_ft'],'elec_ft');


%% Plot the meshes and aligned electrode toegether for sanity check that everything is correct before going further;

%subject ID;
subjects = [3];

s = subjects;
cd(['XXX\headmodel_mat\subj',num2str(s)]);
load(['subj',num2str(s),'_bnd']); load(['subj',num2str(s),'_elec_ft']);load(['subj',num2str(s),'_vol']);


elec_ft = ft_convert_units(elec_ft,'cm');
figure;
ft_plot_mesh(elec_ft.chanpos,'vertexcolor',[1 0 0],'vertexmarker','o') % plot all electrodes as red circles
ft_plot_mesh(vol.bnd(1), 'facealpha', 0.25), hold on
ft_plot_mesh(vol.bnd(2), 'facealpha', 0.25, 'facecolor', 'blue')
ft_plot_mesh(vol.bnd(3), 'facealpha', 0.25, 'facecolor', 'yellow')
hold;
ft_plot_sens(elec_ft ,'style', 'g*');

%% Construct a volume conduction model (dipoli doesn't work for windows - needs to run this step on linux machine from the server);

%subject ID;
subjects = [3];

%Adpat path according to linux machine folders;
s = subjects;
folder = ['XXX\headmodel_mat\subj',num2str(s),'\'];
load([folder,'subj', num2str(s),'_bnd']); load([folder, 'subj', num2str(s),'_elec_ft']);load([folder,'subj', num2str(s),'_vol']);

%Remove the CSF element to allow ft_prepare to work properly (and save the original bnd somewhere);
bnd_original = bnd; 
save(['subj',num2str(s),'_bnd_original'],'bnd_original');
bnd = bnd([1:2 4]);

cfg = [];
cfg.method = 'dipoli'; 
cfg.conductivity = [0.3300 0.0042 0.3300]; %cfg.conductivity = [0.3300 0.0042 1.0000 0.3300];
cfg.tissue = {'scalp', 'skull', 'brain'}; % cfg.tissue = {'scalp', 'skull', 'csf', 'brain'};
vol = ft_prepare_headmodel(cfg, bnd);

%Plot the 3 elements separately;
figure; ft_plot_mesh(vol.bnd(3),'facecolor','none'); %brain
figure; ft_plot_mesh(vol.bnd(2),'facecolor','none'); %skull
figure; ft_plot_mesh(vol.bnd(1),'facecolor','none'); %scalp

%Plot the 3 elements together;
figure;
ft_plot_mesh(vol.bnd(1),'facecolor','none'); %scalp
hold;
ft_plot_mesh(vol.bnd(2), 'facecolor',[0.1 0.1 0.1], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05); %skull
hold ;
% ft_plot_mesh(vol.bnd(3),'edgecolor','none','facealpha',0.4); %csf
% hold;
ft_plot_mesh(vol.bnd(3),'edgecolor','none','facecolor',[0.6 0.4 0.4], 'facealpha', 0.7); %brain
hold;set(gcf, 'color', 'w')

%plot the fine tuned electrodes to the headmodel
ft_plot_mesh(elec_ft.chanpos,'vertexcolor',[1 0 0],'vertexmarker','o') % plot all electrodes as red circles
ft_plot_mesh(vol.bnd(1), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);

%Prepare source model using pre-specifed [x,y,z] grid; 
vol = ft_convert_units(vol,'cm');
elec_ft = ft_convert_units(elec_ft,'cm');

cfg = [];  
cfg.grid.resolution = 1; %cm
cfg.headmodel = vol;
cfg.elec = elec_ft;
grid = ft_prepare_sourcemodel(cfg);

%Plot the head model of the participant for a very last check;
%Plot only the grid positions within the brain;
figure;
ft_plot_mesh(grid.pos(grid.inside,:))
%Plot the BEM;
ft_plot_mesh(vol.bnd(1), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
hold on;
ft_plot_mesh(vol.bnd(2),'edgecolor','none','facealpha',0.4);
hold on;
ft_plot_mesh(vol.bnd(3),'edgecolor','none','facecolor',[0.4 0.6 0.4],'facealpha',0.6);

%Save the individual headmodel templates of the subject (elec position, grid and vol templates);
save([folder,'subj',num2str(s),'_elec_ft'], 'elec_ft');
save([folder,'subj',num2str(s),'_grid'], 'grid');
save([folder,'subj',num2str(s),'_vol'], 'vol');

%% Final check realigned electrodes on the scalp;

%subject ID;
subjects = [3];

s = subjects;

cd(['XXX\headmodel_mat\subj',num2str(s)]);
load(['subj',num2str(s),'_elec_ft']);
load(['subj',num2str(s),'_grid']);
load(['subj',num2str(s),'_vol']);
%Display interactive electrode positions on scalp model;
cfg = [];
cfg.method = 'interactive';
cfg.elec = elec_ft;
cfg.headshape = vol.bnd(1);
ft_electroderealign(cfg);

%% Using a template grid to warp individuals' grid points to the same location in normalized MNI space so that source-constructed activity can
%be directly averaged across subjects.

%subject ID;
subjects = [3];

%Normalise individual mri to MNI space;
s = subjects;
cd(['XXX\headmodel_mat\subj',num2str(s)]); 
mri_rs = ft_read_mri(['subj' num2str(s) '_spm.nii']);
mri_rs.coordsys = 'spm';

%Load the grid template from one participant (e.g. subj2);
cd XXX\fieldtrip\template\headmodel;
load template_grid;

cfg = [];
cfg.grid.warpmni = 'yes';
cfg.grid.template = template_grid;
cfg.grid.nonlinear = 'yes'; % use non-linear normalization
cfg.mri = mri_rs;
grid_n = ft_prepare_sourcemodel(cfg,grid);
grid_n.label = preproc_data.label(1:128,1);

%Plot the head model with the normalized grid positions within the brain;
figure;
ft_plot_mesh(grid.pos(grid.inside,:))
%Plot the BEM;
ft_plot_mesh(vol.bnd(1), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
hold on;
ft_plot_mesh(vol.bnd(2),'edgecolor','none','facealpha',0.4);
hold on;
ft_plot_mesh(vol.bnd(3),'edgecolor','none','facecolor',[0.4 0.6 0.4],'facealpha',0.6);

%save the normalized grid;
cd(['XXX\headmodel_mat\subj',num2str(s)]);
save(['subj',num2str(s),'_grid_n'], 'grid_n');

%%