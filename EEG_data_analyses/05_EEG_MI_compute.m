%% Calculate filters and source activity in Movie Trials;
clearvars; close all;

%subject ID;
subjects = [2 4:9 11:17 19:26]; 

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

%% NOW Calculate Mutual Information in the separate 0.5-2.5s (i.e. equivalent T1) and 2.5-4.5s (i.e. equivalent T2) Time-windows during silent movie presentation at source level;
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
            cfg.toi = 0.5:0.004:2.496; 
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
            cfg.toi = 2.50:0.004:4.5;   
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

%% Here we clean and rearrange the MI data at source level before interpolation; 
%If the last trial (#60) is missing, it fixes the bug by adding an extra zeros line before going further;
clearvars;clc;cd XXX\; 

%subject ID;
subjects = [2:9 11:17 19:26];  

for s = subjects
    
    cd(['XXX\subj',num2str(s)])
    load (['MI_source_pre_subj_', num2str(s)]);
    load (['MI_source_post_subj_', num2str(s)]);
    %save una copia, non si sa mai;
    save(['MI_source_pre_subj_', num2str(s), '_old'], 'all_MI_source_mina_pre','miss_movies'); 
    save(['MI_source_post_subj_', num2str(s), '_old'], 'all_MI_source_mina_post','miss_movies');                                               

    for chnl = 1:length(all_MI_source_mina_pre)

        %Add line of zeros if trial 60 is missing;
        if 60 - length(all_MI_source_mina_pre{chnl})> 0
            all_MI_source_mina_pre{1,chnl} = [all_MI_source_mina_pre{chnl}; zeros((60 - length(all_MI_source_mina_pre{chnl})),length(all_MI_source_mina_pre{chnl}((60 - length(all_MI_source_mina_pre{chnl})),:)))];    
        end

        if 60 - length(all_MI_source_mina_post{chnl})> 0
            all_MI_source_mina_post{1,chnl} = [all_MI_source_mina_post{chnl}; zeros((60 - length(all_MI_source_mina_post{chnl})),length(all_MI_source_mina_post{chnl}((60 - length(all_MI_source_mina_post{chnl})),:)))];    
        end

    end

    %save the new ones;
    cd(['XXX\subj',num2str(s)]);
    save(['MI_source_pre_subj_', num2str(s)], 'all_MI_source_mina_pre','miss_movies'); 
    save(['MI_source_post_subj_', num2str(s)], 'all_MI_source_mina_post','miss_movies'); 

    clearvars -except subjects
    
end

%% Realign the MI_source on FreqPeak for each participant (Freqpeak = frequency in which the MI peaks between lips and audio envelope signals, see stimuli analyses);
clearvars; clc; cd XXX\; 
load ('stim_freqpeak.mat'); %this matrix contains the list of the 60 videos and their frequency of MI peak between lips and sound (in AV_info_stimuli folder);

%subject ID;
subjects = [2:9 11:17 19:26];  

for s = subjects 
    
    %Load the MI lips/speech eeg data of the subject;
    cd(['XXX\subj',num2str(s)])
    load (['MI_source_pre_subj_', num2str(s)]);
    load (['MI_source_post_subj_', num2str(s)]);

    %Keep the MI only for the frequency peak of the video/sound; 
    for chnl = 1:length(all_MI_source_mina_pre)

        for i = 1:length(all_MI_source_mina_pre{1,chnl})   

            MI_source_mina_peak_pre{1,chnl}(i,:) = all_MI_source_mina_pre{chnl}(i, stim_freqpeak(i,2)-3:stim_freqpeak(i,2)+1);
            MI_source_mina_peak_post{1,chnl}(i,:) = all_MI_source_mina_post{chnl}(i, stim_freqpeak(i,2)-3:stim_freqpeak(i,2)+1);

        end
    end

    %Save realigned MI_source info on MI freqpeak;
    cd(['XXX\subj',num2str(s)]);
    save (['MI_source_pre_subj_', num2str(s), '_peak'],'MI_source_mina_peak_pre','miss_movies');
    save (['MI_source_post_subj_', num2str(s), '_peak'],'MI_source_mina_peak_post','miss_movies');

    clearvars -except stim_freqpeak subjects

end

%% Give MI_source a Fieldtrip structure for source localization;
clearvars; clc;

%subject ID;
subjects = [2:9 11:17 19:26];      

allsubj_MI_source_mina_peak_pre = cell(1,max(subjects)); allsubj_MI_source_mina_peak_pre{1,length(subjects)} = [];
allsubj_MI_source_mina_peak_post = cell(1,max(subjects)); allsubj_MI_source_mina_peak_post{1,length(subjects)} = [];

for s = subjects

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %--------- 1 LOAD ALL WHAT WE NEED -------
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Prepare Movie logfiles of the subject;
    preproc_dir = ['XXX\subj_', num2str(s)];
    cd(preproc_dir);
    % Load the Sound and Movie logfiles from the recording session;
    logfile_dir = ['XXX\subj_', num2str(s)];
    cd(logfile_dir);
    %Movies Logfile;
    moviefile_temp = dataset('File', ['subj',num2str(s),'-Movie.csv'], 'Delimiter' ,',');
    moviefile = moviefile_temp.TRIAL;
    moviefile(:,2) = moviefile_temp.STIM_NUMB;
    moviefile(:,3) = moviefile_temp.FREQPEAK;
    clear moviefile_temp;

    % Load source data and MI_source of the subject;
    cd(['XXX\subj',num2str(s)]);
    load (['MI_source_pre_subj_', num2str(s), '_peak']); load (['MI_source_post_subj_', num2str(s), '_peak']); load (['subj', num2str(s), '_source_movie']);
    cd(['XXX\headmodel_mat\subj',num2str(s)]);
    load(['subj', num2str(s), '_elec_ft']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- 1 MUTUAL INFORMATION Movies-EEGs ---
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% MINA INFOS %%%
    mi_mina_pre_temp = [];
    mi_mina_pre_temp.label = source_movie.label;
    mi_mina_pre_temp.freq = 1:5;
    mi_mina_pre_temp.dimord = 'rpt_chan_freq';
    mi_mina_post_temp = [];
    mi_mina_post_temp.label = source_movie.label;
    mi_mina_post_temp.freq = 1:5;
    mi_mina_post_temp.dimord = 'rpt_chan_freq';

    %MI Mina pre (T1 window) and post (T2 window);
    mi_mina_data_pre = zeros(size(MI_source_mina_peak_pre{1},1),numel(mi_mina_pre_temp.label),length(mi_mina_pre_temp.freq));
    mi_mina_data_post = zeros(size(MI_source_mina_peak_post{1},1),numel(mi_mina_post_temp.label),length(mi_mina_post_temp.freq));

    for i = 1:numel(MI_source_mina_peak_pre)    
        for j = 1:size(MI_source_mina_peak_pre{1},1)
            mi_mina_data_pre(j,i,:) = MI_source_mina_peak_pre{i}(j,:);
        end    
    end
        mi_mina_pre_temp.powspctrm = mi_mina_data_pre;
        mi_mina_pre_temp.cumsumcnt(:,1) = length(MI_source_mina_peak_pre{1,1}).*ones(size(MI_source_mina_peak_pre{1},1),1);
        mi_mina_pre_temp.cumtapcnt(:,1) = ones(size(MI_source_mina_peak_pre{1},1),1);
        mi_mina_pre_temp.elec = elec_ft;

    for i = 1:size(mi_mina_pre_temp.powspctrm,1)
        mi_mina_pre_temp.trialinfo{i,1}.condition = 2;
        mi_mina_pre_temp.trialinfo{i,1}.trialnum = moviefile(moviefile(:,2) == i,1);
        mi_mina_pre_temp.trialinfo{i,1}.trialstim_num = moviefile(moviefile(:,2) == i,2);
        mi_mina_pre_temp.trialinfo{i,1}.freq_peak = moviefile(moviefile(:,2) == i,3);
    end


    for i = 1:numel(MI_source_mina_peak_post)    
        for j = 1:size(MI_source_mina_peak_post{1},1)
            mi_mina_data_post(j,i,:) = MI_source_mina_peak_post{i}(j,:);
        end    
    end
        mi_mina_post_temp.powspctrm = mi_mina_data_post;
        mi_mina_post_temp.cumsumcnt(:,1) = length(MI_source_mina_peak_post{1,1}).*ones(size(MI_source_mina_peak_post{1},1),1);
        mi_mina_post_temp.cumtapcnt(:,1) = ones(size(MI_source_mina_peak_post{1},1),1);
        mi_mina_post_temp.elec = elec_ft;

    for i = 1:size(mi_mina_post_temp.powspctrm,1)
        mi_mina_post_temp.trialinfo{i,1}.condition = 2;
        mi_mina_post_temp.trialinfo{i,1}.trialnum = moviefile(moviefile(:,2) == i,1);
        mi_mina_post_temp.trialinfo{i,1}.trialstim_num = moviefile(moviefile(:,2) == i,2);
        mi_mina_post_temp.trialinfo{i,1}.freq_peak = moviefile(moviefile(:,2) == i,3);
    end

    %Remove missing movie only trials;
    cfg = [];
    cfg.trials = find(cell2mat(cellfun(@(x) ~ismember(x.trialstim_num,miss_movies), mi_mina_pre_temp.trialinfo, 'UniformOutput', false)));
    mi_mina_pre = ft_selectdata(cfg,mi_mina_pre_temp);
    cfg = [];
    cfg.trials = find(cell2mat(cellfun(@(x) ~ismember(x.trialstim_num,miss_movies), mi_mina_post_temp.trialinfo, 'UniformOutput', false)));
    mi_mina_post = ft_selectdata(cfg,mi_mina_post_temp);

    %Average MI_mina_eeg data across freqpeak (i.e. 4-8Hz theta) and per frequency peak; theta(4-8Hz);
    cfg = [];
    cfg.keeptrials = 'yes';
    cfg.parameter = 'powspctrm';
    MI_mina_peak_pre = ft_freqdescriptives(cfg,mi_mina_pre);
    cfg = [];
    cfg.parameter = 'powspctrm';
    MI_mina_peak_post = ft_freqdescriptives(cfg,mi_mina_post);

    %Save together for further grand averages;
    allsubj_MI_source_mina_peak_pre{1,s} = MI_mina_peak_pre;
    allsubj_MI_source_mina_peak_post{1,s} = MI_mina_peak_post;

    clear MI_mina_peak_pre MI_mina_peak_post
    
end

%Save all the MI info for all the subjects together;
MI_source_dir = 'XXX';
save([MI_source_dir, 'allsubj_MI_source_peak_VISUAL_n', num2str(length(subjects))], 'allsubj_MI_source_mina_peak_pre', 'allsubj_MI_source_mina_peak_post');
