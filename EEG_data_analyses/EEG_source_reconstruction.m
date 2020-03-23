%% Source reconstruction at determined left audio, left visual and right visual ROIs;
clearvars;clc;

%subjects ID;
subjects = [2 4:9 11:17 19:26];
% subjects = [3]; %FOR SUBJECT 3, CHANGE grid by grid_n (to use the MNI normalized grid);

for s = subjects
    
%load preproc data;
cd (['XXX\subj_' num2str(s)]);
load(['subj_',num2str(s), '_after_preproc']);    

%Load Template 128-Easycap headmodel;
cd (['XXX\headmodel_mat\subj' num2str(s)]);
load(['subj',num2str(s),'_elec_ft']);
load(['subj',num2str(s),'_grid']);
%load(['subj',num2str(s),'_grid_n']); % <--- if subj3, LOAD grid_n; 
load(['subj',num2str(s),'_vol']);

%Time-lock the data get common source filters
cfg = [];
cfg.keeptrials = 'yes';
erp = ft_timelockanalysis(cfg, preproc_data);

%Keep left/right hemispheres seperately;
cfg = [];
cfg.layout = 'biosemi128';
lay = ft_prepare_layout(cfg, []);

cfg = []; 
cfg.method = 'lcmv';
cfg.grid = grid;
%cfg.grid = grid_n; %<-- if normalized in MNI (subj3);
cfg.lcmv.lambda = 5;
cfg.headmodel = vol;
cfg.elec = elec_ft;
cfg.channel = 1:128;
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.fixedori = 'yes';
source_tmp = ft_sourceanalysis(cfg, erp);
filters = cell2mat(source_tmp.avg.filter(source_tmp.inside)); % keep filters

%Create a data structure for the source data visual;
source = [];
source.sampleinfo = preproc_data.sampleinfo; % transfer sample information
source.time = preproc_data.time;       % transfer time information
source.trialinfo = preproc_data.trialinfo;  % transfer trial information
source.fsample = preproc_data.fsample;     

% create labels for each virtual electrode
for c = 1 : sum(grid.inside)
    label{c,1} = ['S' num2str(c)]; 
end

source.label = label;

% for each trial, apply filters to the recorded data
for j = 1 : numel(preproc_data.trial) 
    source.trial{1,j} = single(filters*preproc_data.trial{1,j}(1:128,:)); 
end

%WE NEED TO DEFINE the ROIS: AR, AL and VIS; 
audio_left = 1743;
visual_left = 1639;
visual_right = 1947;
% control = 3668;

cfg = [];
cfg.channel = sum(grid.inside(1:audio_left)); % <--- if subj3, use grid_n; 
source_audio_left = ft_preprocessing(cfg,source);

cfg = [];
cfg.channel = sum(grid.inside(1:visual_left)); % <--- if subj3, use grid_n; 
source_visual_left = ft_preprocessing(cfg,source);

cfg = [];
cfg.channel = sum(grid.inside(1:visual_right)); % <--- if subj3, use grid_n; 
source_visual_right = ft_preprocessing(cfg,source);

% cfg = [];
% cfg.channel = sum(grid.inside(1:control)); % <--- if subj3, use grid_n; 
% source_control = ft_preprocessing(cfg,source);


%save Data_source for each participant;
cd XXX;
mkdir(['subj' num2str(s)]);
cd (['XXX\subj',num2str(s)]);
save (['source_reconstruction_subj' num2str(s)], 'source_audio_left','source_visual_left','source_visual_right','filters');

end

%% Compute the visual ERP at the 3 reconstructed sources level in the movie condition for visual inspection and flipping correction;
clearvars;clc;

%subjects ID;
subjects = [2:9 11:17 19:26];


for s = subjects
   
    cd (['XXX\subj' num2str(s)]);
    load(['source_reconstruction_subj' num2str(s)]); 

    %Compute the visual ERP at the 3 reconstructed sources level in the movie condition;
    %lowpass filter the eeg data at the 3 sources of interest;
    cfg = [];
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 30;
    source_audio_left = ft_preprocessing(cfg, source_audio_left);
    source_visual_left = ft_preprocessing(cfg, source_visual_left);
    source_visual_right = ft_preprocessing(cfg, source_visual_right);
    %select silent movie trials only;
    cfg = [];
    cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.condition,2), source_audio_left.trialinfo, 'UniformOutput', false))); 
    erp_audio_left = ft_timelockanalysis(cfg, source_audio_left);
    erp_visual_left = ft_timelockanalysis(cfg, source_visual_left);
    erp_visual_right = ft_timelockanalysis(cfg, source_visual_right);

    %Compute the visual ERP at scalp level with the movieonly condition;
    cd (['XXX\subj_' num2str(s)]);
    load(['subj_',num2str(s), '_after_preproc']);  

    %lowpass filter the eeg data;
    cfg = [];
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 30;
    preproc_data = ft_preprocessing(cfg, preproc_data);
    %Movie only ERP;   
    cfg = [];
    cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.condition,2), preproc_data.trialinfo, 'UniformOutput', false))); 
    erp_movie = ft_timelockanalysis(cfg, preproc_data);

    cd (['XXX\subj' num2str(s)]);
    save (['ERPs_source_subj' num2str(s)], 'erp_audio_left', 'erp_visual_left', 'erp_visual_right', 'erp_movie');
   
end

%% Find the best channels to detect the visual component;
clearvars;clc;

%subjects ID;
subjects = [2:9 11:17 19:26];

%creat your matrix to store the channels capturing best the visual component;
FLIP_INSPECT = cell(length(subjects),6);
FLIP_INSPECT = array2table(FLIP_INSPECT, 'Variable',{'PARTICIPANT','uni_vis_ch','ind_lat_P2','flip_left_audio','flip_left_visual','flip_right_visual'}); 

%Here we make a bit of automatic detection of the visual component; It is supposed to detect the channel with the bigger component but in the following section you can control it visually for each participant and
%correct it to pick the best channel with the best visual component;

ind_lat_p2 = [0.08 0.12]; %choose the TW containing the visual component;
channel_subject = cell(length(subjects), 1); 
channel_subject{length(subjects), 1} = [];

for s = subjects
    
    cd (['XXX\subj' num2str(s)]);
    load (['ERPs_source_subj' num2str(s)]);
       
    %Select time-window of P2 visual component; 
    cfg = [];
    cfg.latency = ind_lat_p2;
    cfg.avgovertime = 'yes';
    cfg.channel = 1:128;
    uni_p1 = ft_selectdata(cfg, erp_movie);
    
    %Find the channel with the bigger P1 visual component;   
    [~, uni_vis_ch] = max(uni_p1.avg);
    channel_subject{s,1} = erp_movie.label{uni_vis_ch};
         
end

cd XXX\;
save FLIP_INSPECTION FLIP_INSPECT;  

%% Visual inspection of the ERPs of interest for every subject individually to localise the best electrode for viusal component and its time-window;
clearvars;clc;
cd XXX\;
load FLIP_INSPECTION; 

%subject ID;
subjects = 2;

%Now inspect each participant individually and find the best channel for visual component;
%Change it in the FLIP_INSPECT corresponding cell;

for s = subjects
    
    cd (['XXX\subj' num2str(s)]);
    load (['ERPs_source_subj' num2str(s)]);
    close all;
    cfg = [];
    cfg.parameter = 'avg';
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 8;
    cfg.xlim = [0 1];
    cfg.layout = 'biosemi128';
    cfg.marker = 'labels';
    figure; ft_topoplotER(cfg,erp_movie);

end

%save the updated FLIP INSPECTION after individual correction;
cd XXX;
save FLIP_INSPECTION FLIP_INSPECT; 

%% This is the section were we determine the direction of the dipole which may have been randomized after the beamforming reconstruction;
%for each participant, we will compare the visual ERP direction at scalp level and source level (at the 3 ROIs) to match the original direction and adjust by -1
%if the polarity has been inverted by the beamformer method;
clearvars;clc;
cd XXX\;
load FLIP_INSPECTION; 


%subject ID;
subjects = 2; 

%load subject ERPs at source level;  
s = subjects;
cd (['XXX\subj' num2str(s)]);
load (['ERPs_source_subj' num2str(s)]);
%select movie scalp erps;
cfg = [];
cfg.channel = FLIP_INSPECT_visual_entrainment.uni_vis_ch{s};
erp_uni_vis = ft_selectdata(cfg, erp_movie);

%Compute the ERPs at sources, in multimodal e unimodal conditions to compare the polarity;
cfg = [];
cfg.lpfilter = 'yes';
cfg.lpfreq = 15;
%scalp ERPs;
erp_uni_vis = ft_preprocessing(cfg, erp_uni_vis);
%source ERPs;
erp_audio_left = ft_preprocessing(cfg,erp_audio_left);
erp_visual_left = ft_preprocessing(cfg, erp_visual_left);
erp_visual_right = ft_preprocessing(cfg, erp_visual_right);
%normalize the signals;
erp_uni_vis.avg = (erp_uni_vis.avg - mean(erp_uni_vis.avg))./std(erp_uni_vis.avg);
erp_audio_left.avg = (erp_audio_left.avg - mean(erp_audio_left.avg))./std(erp_audio_left.avg);
erp_visual_left.avg = (erp_visual_left.avg - mean(erp_visual_left.avg))./std(erp_visual_left.avg);
erp_visual_right.avg = (erp_visual_right.avg - mean(erp_visual_right.avg))./std(erp_visual_right.avg);
   
%Always start inspections with +1 indice to make it easy;
%if polarity of the ERP at the sourc eis the same as the scalp ERP, flip_xxx = 1; if polarity of the ERP at the source is inversed to the scalp ERP, flip_xxx = -1; 
flip_la = 1;
flip_lv = 1;
flip_rv = 1;
    
%plot the signal together with the average scalp ERPs: Play around the parameters to make sure you can determine the right polarity respect to the scalp ERPs for each source ERPs; 
close all; figure; 
plot(time,flip_la*erp_audio_left.avg(:,512:2048), 'LineWidth',1, 'Color', [0 0 1]);
xlim([-0.1 1]);
hold on 
plot(time,flip_lv*erp_visual_left.avg(:,512:2048), 'LineWidth',1, 'Color', [1 0.6 0]);
hold on 
plot(time,flip_rv*erp_visual_right.avg(:,512:2048), 'LineWidth',1, 'Color', [1 0 0]);
hold on
plot(time,erp_uni_vis.avg(:,512:2048), 'LineWidth',3, 'Color', [0 0 0]);
legend('Audio left','Visual left','Visual right','ERP scalp');

    
%Once you have determined the polarity of each source ERP, put 1 or -1 in the corresponding cell of the FLIP_INSPECTION table;
%Then, save (my advice is to save after each subject); 
cd XXX\;
FLIP_INSPECT_visual_entrainment.flip_left_audio{s} = flip_la;
FLIP_INSPECT_visual_entrainment.flip_left_visual{s} = flip_lv;
FLIP_INSPECT_visual_entrainment.flip_right_visual{s} =flip_rv;
save FLIP_INSPECTION FLIP_INSPECT;

    
%% Now Flip the source data by 1 or -1 to correct the dipole directionality;
clearvars;clc;
cd XXX\;
load FLIP_INSPECTION; 

%subjects ID;
subjects = [2:9 11:17 19:26];


for s = subjects

    cd (['XXX\subj' num2str(s)]);
    load(['source_reconstruction_subj' num2str(s)])
    source_left_audio_realigned = source_audio_left;
    source_left_visual_realigned = source_visual_left;
    source_right_visual_realigned = source_visual_left;
    
    for i =1:length(source_left_audio_realigned.trial)
        source_left_audio_realigned.trial{i} = source_left_audio_realigned.trial{i}*FLIP_INSPECT.flip_left_audio{s}; 
        source_left_visual_realigned.trial{i} = source_left_visual_realigned.trial{i}*FLIP_INSPECT.flip_left_visual{s}; 
        source_right_visual_realigned.trial{i} = source_right_visual_realigned.trial{i}*FLIP_INSPECT.flip_right_visual{s};

    end
        
    save (['source_reconstruction_realigned_subj' num2str(s)], 'source_left_audio_realigned','source_left_visual_realigned','source_right_visual_realigned');      
        
end

%% Finally, we look at the reconstructed ERPs at the ROI sources before and after the flipping realignment to see any improvement;
%This is not exact science so maybe sometimes there is not obvious improvement - patience;
clearvars;clc;
cd XXX\;
load FLIP_INSPECTION; 

%subjects ID;
subjects = [2:9 11:17 19:26];

sc = 0;
for s = subjects

    cd (['XXX\subj' num2str(s)]);
    load (['ERPs_source_subj' num2str(s)]);

    sc = sc + 1;
    erp_before.erp_left_audio{sc,1} = erp_audio_left;
    erp_before.erp_left_visual{sc,1} = erp_visual_left;
    erp_before.erp_right_visual{sc,1} = erp_visual_right;
    
    load (['source_reconstruction_realigned_subj' num2str(s)]);

    cfg = [];
    cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.condition,2), source_left_audio_realigned.trialinfo, 'UniformOutput', false))); 
    erp_left_audio_realigned = ft_timelockanalysis(cfg, source_left_audio_realigned);
    erp_left_visual_realigned = ft_timelockanalysis(cfg, source_left_visual_realigned);
    erp_right_visual_realigned = ft_timelockanalysis(cfg, source_right_visual_realigned);    
    
    erp_after.erp_left_audio{sc,1} = erp_left_audio_realigned;
    erp_after.erp_left_visual{sc,1} = erp_left_visual_realigned;
    erp_after.erp_right_visual{sc,1} = erp_right_visual_realigned;


end    

%Grande average of ERps before re-aligning (reflip); 
cfg = [];
GA_before_left_audio = ft_timelockgrandaverage(cfg,erp_before.erp_left_audio{:});
GA_before_left_visual = ft_timelockgrandaverage(cfg,erp_before.erp_left_visual{:});
GA_before_right_visual = ft_timelockgrandaverage(cfg,erp_before.erp_right_visual{:});
%Grande average of ERps after re-aligning (reflip); 
cfg = [];
GA_after_left_audio = ft_timelockgrandaverage(cfg,erp_after.erp_left_audio{:});
GA_after_left_visual = ft_timelockgrandaverage(cfg,erp_after.erp_left_visual{:});
GA_after_right_visual = ft_timelockgrandaverage(cfg,erp_after.erp_right_visual{:});
    
cfg = [];
close all;
cfg.preproc.lpfilter = 'yes';
cfg.lpfreq = 30;
cfg.linewidth = 1;
figure; ft_singleplotER(cfg, GA_before_left_audio, GA_after_left_audio);
legend('before realigned','after');
title('ERP Left Audio - S1743');

figure; ft_singleplotER(cfg, GA_before_left_visual, GA_after_left_visual);
legend('before realigned','after');
title('ERP Left Visual  - 1639');

figure; ft_singleplotER(cfg, GA_before_right_visual);
legend('before realigned','after');
title('ERP Right Visual - 1947');
   
%%