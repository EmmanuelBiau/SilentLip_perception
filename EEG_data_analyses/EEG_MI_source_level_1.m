%% Calculate filters and source activity in Movie Trials;
clearvars; close all;

%subject ID;
subjects = [2 4:9 11:17 19:26]; 
% subjects = [3]; 

for s = subjects

%Load preprocessed data;
cd(['XXX\subj_',num2str(s)]);
load(['subj_',num2str(s),'_after_preproc']);

%Load participant's headmodel;
cd(['XXX\headmodel_mat\subj',num2str(s)]);
load(['subj',num2str(s),'_elec_ft']);
load(['subj',num2str(s),'_grid']); 
% load(['subj',num2str(s),'_grid_n']); % <--- if subject 3 with MRI scans, LOAD grid_n; 
load(['subj',num2str(s),'_vol']);
%convert elec coordinates in mm;
elec_ft = ft_convert_units(elec_ft,'mm');

%downsample eeg signal to 250Hz (as Lip signal from stim);
cfg = [];
cfg.resamplefs = 250;
preproc_data = ft_resampledata(cfg, preproc_data);

%Select movie trials;
cfg = [];
cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.condition,2), preproc_data.trialinfo, 'UniformOutput', false)));
movie_trials = ft_preprocessing(cfg,preproc_data);

%Compute first the normal leadfield (for movie condition later);
cfg = [];
cfg.grid = grid; 
% cfg.grid = grid_n; %<-- if normalized in MNI (subject #3 with T1 scans);
cfg.vol = vol;
cfg.elec = elec_ft;
grid = ft_prepare_leadfield(cfg,preproc_data);

%Timelock movie trials;
cfg = [];
cfg.keeptrials = 'yes';
cfg.covariance = 'yes';
cfg.covariancewindow = 'all';
tlmovie_data = ft_timelockanalysis(cfg, movie_trials);

%Move the timelocked movie data into source space and create the filter -MOVIES;
cfg = []; 
cfg.method = 'lcmv';
cfg.grid = grid;
cfg.headmodel = vol;
cfg.elec = elec_ft;
cfg.channel = 1:128;
cfg.lcmv.keepfilter = 'yes';
cfg.lcmv.fixedori = 'yes';
cfg.lcmv.lambda = 5;
source_tmp_movie = ft_sourceanalysis(cfg, tlmovie_data);
filters_movie = cell2mat(source_tmp_movie.avg.filter(source_tmp_movie.inside)); % keep filters

%save the movie filters;
cd XXX; 
mkdir(['subj',num2str(s)]); 
cd(['XXX\subj',num2str(s)]);
save(['subj',num2str(s),'_filters_movie'],'filters_movie');

%Create a source data structure for movie data;
source_movie = [];
source_movie.sampleinfo = movie_trials.sampleinfo;
source_movie.time = movie_trials.time; 
source_movie.trialinfo = movie_trials.trialinfo;
source_movie.fsample = movie_trials.fsample;  

%Create labels for each virtual electrode in source movie data;
for c = 1 : sum(grid.inside) 
    label{c,1} = ['S' num2str(c)]; 
end
source_movie.label = label;

for j = 1 : numel(movie_trials.trial) 
    source_movie.trial{1,j} = single(filters_movie*movie_trials.trial{1,j}(1:128,:)); 
end

%Save the final sound and movie source data, ready for MI at source level analyses;
cd(['XXX\subj',num2str(s)]) 
save(['subj',num2str(s),'_source_movie'],'source_movie');

clearvars -except subjects

end

%% Sort/prepare the epochs according to the corresponding Lips movement signal;
clearvars; clc; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------PREPARE SUBJECT'S TRIALS ---------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%subjects ID;
subjects = [2:9 11:17 19:26];


for s = subjects 

%Load headmodels from source localization in Movie condition;
cd(['XXX\subj',num2str(s)]);    
load(['subj',num2str(s),'_source_movie'],'source_movie'); % Head model from Movie condition; 
    
%Cut epoch length to stimuli length(0-5s);
cfg = [];
cfg.trials = 'all';
cfg.begsample = 501;  % t = 0s; 
cfg.endsample = 1750; % t = 5s;
source_movie_data_cut = ft_redefinetrial(cfg, source_movie);

%Reorganize the EEG data of Movie only block;
for i = 1:numel(source_movie_data_cut.trial)
    for j = 1:size(source_movie_data_cut.trial{1},1)
        source_movie_data_cut2{j}(i,:) = source_movie_data_cut.trial{i}(j,:);
    end
end

%Find missing trials;
for i = 1:numel(source_movie_data_cut.trialinfo)
    videonum(i,:) = source_movie_data_cut.trialinfo{i}.trialstim_num;
end

allvideonum = 1:60;
miss_movies = find(~ismember(allvideonum , videonum));

%Add NaN values for each missing trial;
for i = 1:numel(source_movie_data_cut2)
    source_movie_data_cut2{i} = [source_movie_data_cut2{i}; NaN(length(miss_movies),size(source_movie_data_cut2{i},2))];
end
  
newvideonum = [videonum; miss_movies'];
[~,sortvideoidx] = sort(newvideonum);

%Sort eeg_movies according to sortvideoidx;
for i = 1:numel(source_movie_data_cut2)
    eeg_source_movies{i} = source_movie_data_cut2{i}(sortvideoidx,:);
end

cd(['XXX\subj',num2str(s)]) 
save(['subj',num2str(s),'_source_movie_cut'],'eeg_source_movies','miss_movies','source_movie_data_cut');

clearvars -except subjects

end


%% NOW Calculate Mutual Information in the separate 0-2.5s (i.e. equivalent T1) and 2.5-5s (i.e. equivalent T2) Time-windows during silent movie presentation at source level;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------- NOW WE CAN PERFORM MI CALCULATIONS AT SOURCE LEVEL - T1 Time-window --------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars; clc; addpath XXX\AV_info_stimuli\gauss_info\; addpath XXX\AV_info_stimuli\gauss_info\mex\;


%subjects ID;
subjects = [2:9 11:17 19:26];

for s = subjects 

cd(['XXX\subj',num2str(s)])    
load(['subj',num2str(s),'_source_movie_cut']);  

%Now we calculate the MI between Movies and source epochs; Load your lips signal here;
stim_info_dir = 'XXX\AV_info_stimuli\';
cd(stim_info_dir);

list_v = dir('Info_lips_*.mat');

for a = 1:length(list_v)
    list_v_num(a) = str2num(list_v(a).name(16:end-4));
end

[newnum, inx] = sort(list_v_num);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- MUTUAL INFORMATION Movies-source_EEGs ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Prepare the list of remaining movie trials after processing;  
newnum_valid_movies = newnum;
newnum_valid_movies(miss_movies) = [];

%Create a source data structure
source_template = [];
source_template.sampleinfo = source_movie_data_cut.sampleinfo(1,:);
source_template.time{1} = source_movie_data_cut.time{1};
source_template.trialinfo = 1;
source_template.fsample = source_movie_data_cut.fsample;  

    for chnl = 1:length(eeg_source_movies)

        for a = 1:length(newnum)

            load(list_v(newnum(a)).name,'ds_mina');
            dat_mina(a,:) = ds_mina; clear ds_mina 

        end

        %Restructure the data in a FT template: this data structure has the same sampling rate (250 Hz);
        dat3 = []; dat3 = source_template; dat3.trial{1,1} = dat_mina;
        dat4 = []; dat4 = source_template; dat4.trial{1,1} = eeg_source_movies{1,chnl};

        dat3.time{1,1} = source_template.time{1,1}(1:1250);
        dat4.time{1,1} = source_template.time{1,1}(1:1250);

        %This label stuff is just for data structure. It create labels for each virtual electrode;
        for c = 1:size(dat3.trial{1},1)

            dat3.label{c,1} = ['T' num2str(c)];
            dat4.label{c,1} = ['T' num2str(c)];

        end

        %Calculate the instantaneous phase for each freq.bin in lips/eeg time-series;
        cfg = [];
        cfg.method = 'wavelet';
        cfg.output = 'fourier';
        cfg.foi= 2:1:20;   
        cfg.toi = 0:0.004:2.496; 
        cfg.width = 5;
        
        freq3 = ft_freqanalysis(cfg,dat3);
        freq4 = ft_freqanalysis(cfg,dat4);

        %Phase information for movies and source_eeg;   
        tfdata3 = freq3.fourierspctrm;
        phase3 = tfdata3./abs(tfdata3);
        phase3 = squeeze(phase3);% mina phase info;

        tfdata4 = freq4.fourierspctrm;
        phase4 = tfdata4./abs(tfdata4);
        phase4 = squeeze(phase4);% eeg phase info;

        %Calculate MI speech envelope information with EEG;

        for f = 1:length(freq3.freq) 

            phase3_temp = squeeze(phase3(:,f,:));
            phase4_temp = squeeze(phase4(:,f,:));

            for c = newnum_valid_movies 

                dat_mina = []; dat_mina = phase3_temp(c,:);
                dat_movie_eeg = []; dat_movie_eeg = phase4_temp(c,:);

                cdat_mina = []; cdat_mina = copnorm([real(dat_mina)' imag(dat_mina)']);
                cdat_movie_eeg = []; cdat_movie_eeg = copnorm([real(dat_movie_eeg)' imag(dat_movie_eeg)']);

                MI_source_mina_eeg(c,f) = info_gg(cdat_movie_eeg, cdat_mina, true, false, false);

            end
        end

        %Feed matrix with MI lips-eeg (and add zero values for the missing trials);
        all_MI_source_mina_pre{chnl} = MI_source_mina_eeg;

        clear dat_mina dat_movie_eeg 

    end

%save MI between Lips and EEG for each participant (pre = T1 window = 0 - 2.5s);
cd(['XXX\subj',num2str(s)])
save (['MI_source_pre_subj_', num2str(s)],'all_MI_source_mina_pre','miss_movies');

clearvars -except subjects addpath
                                         
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------- NOW WE CAN PERFORM MI CALCULATIONS AT SOURCE LEVEL - T2 Time-window --------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; addpath XXX\AV_info_stimuli\gauss_info\; addpath XXX\AV_info_stimuli\gauss_info\mex\;

%subjects ID;
subjects = [2:9 11:17 19:26];   

for s = subjects 

cd(['XXX\subj',num2str(s)])    
load(['subj',num2str(s),'_source_movie_cut']);  

%Now we calculate the MI between Movies and source epochs; Load your lips signal here;
stim_info_dir = 'XXX\AV_info_stimuli\';
cd(stim_info_dir);

list_v = dir('Info_lips_stim_*.mat');

for a = 1:length(list_v)
    list_v_num(a) = str2num(list_v(a).name(16:end-4));
end

[newnum, inx] = sort(list_v_num);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- MUTUAL INFORMATION Movies-source_EEGs ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Prepare the list of remaining movie trials after processing;  
newnum_valid_movies = newnum;
newnum_valid_movies(miss_movies) = [];

%create a source data structure
source_template = [];
source_template.sampleinfo = source_movie_data_cut.sampleinfo(1,:); 
source_template.time{1} = source_movie_data_cut.time{1}; 
source_template.trialinfo = 1; 
source_template.fsample = source_movie_data_cut.fsample;  

    for chnl = 1:length(eeg_source_movies)

        for a = 1:length(newnum)

            load(list_v(newnum(a)).name,'ds_mina');
            dat_mina(a,:) = ds_mina; clear ds_mina 

        end

        %Restructure the data in a FT template: this data structure has the same sampling rate (250 Hz)
        dat3 = []; dat3 = source_template; dat3.trial{1,1} = dat_mina;
        dat4 = []; dat4 = source_template; dat4.trial{1,1} = eeg_source_movies{1,chnl};

        dat3.time{1,1} = source_template.time{1,1}(1:1250);
        dat4.time{1,1} = source_template.time{1,1}(1:1250);

        %This label stuff is just for data structure. It create labels for each virtual electrode
        for c = 1:size(dat3.trial{1},1)

            dat3.label{c,1} = ['T' num2str(c)];
            dat4.label{c,1} = ['T' num2str(c)];

        end

        %Calculate the instantaneous phase for each freq.bin in audio/eeg time-series;
        cfg = [];
        cfg.method = 'wavelet';
        cfg.output = 'fourier';
        cfg.foi= 2:1:20;  
        cfg.toi = 2.50:0.004:4.996;   
        cfg.width = 5;

        freq3 = ft_freqanalysis(cfg,dat3);
        freq4 = ft_freqanalysis(cfg,dat4);

        %Phase information for movies and source_eeg;   
        tfdata3 = freq3.fourierspctrm;
        phase3 = tfdata3./abs(tfdata3);
        phase3 = squeeze(phase3); 

        tfdata4 = freq4.fourierspctrm;
        phase4 = tfdata4./abs(tfdata4);
        phase4 = squeeze(phase4);

        %Calculate MI speech envelope information with EEG;

        for f = 1:length(freq3.freq) % Here from 2-20Hz;

            phase3_temp = squeeze(phase3(:,f,:));
            phase4_temp = squeeze(phase4(:,f,:));

            for c = newnum_valid_movies 

                dat_mina = []; dat_mina = phase3_temp(c,:);
                dat_movie_eeg = []; dat_movie_eeg = phase4_temp(c,:);

                cdat_mina = []; cdat_mina = copnorm([real(dat_mina)' imag(dat_mina)']);
                cdat_movie_eeg = []; cdat_movie_eeg = copnorm([real(dat_movie_eeg)' imag(dat_movie_eeg)']);

                MI_source_mina_eeg(c,f) = info_gg(cdat_movie_eeg, cdat_mina, true, false, false);

            end
        end

        %Feed matrix with MI lips-eeg (and add zero values for the missing trials);
        all_MI_source_mina_post{chnl} = MI_source_mina_eeg;

        clear dat_mina dat_movie_eeg 

    end

%save MI between Lips and EEG for each participant (post = T2 window = 2.5-5s);
cd(['XXX\subj',num2str(s)])
save (['MI_source_post_subj_', num2str(s)],'all_MI_source_mina_post','miss_movies');

clearvars -except subjects addpath
                                         
end

%%