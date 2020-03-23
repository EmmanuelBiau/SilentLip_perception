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
save (['MI_source_pre_subj_', num2str(s), '_old'], 'all_MI_source_mina_pre','miss_movies'); 
save (['MI_source_post_subj_', num2str(s), '_old'], 'all_MI_source_mina_post','miss_movies');                                               

for chnl = 1:length(all_MI_source_mina_pre)

    %Add line of zeros if trial 60 is missing;
    if 60 - length(all_MI_source_mina_pre{chnl})> 0
        all_MI_source_mina_pre{1,chnl} = [all_MI_source_mina_pre{chnl}; zeros((60 - length(all_MI_source_mina_pre{chnl})),length(all_MI_source_mina_pre{chnl}((60 - length(all_MI_source_mina_pre{chnl})),:)))];    
    end

    if 60 - length(all_MI_source_mina_post{chnl})> 0
        all_MI_source_mina_post{1,chnl} = [all_MI_source_mina_post{chnl}; zeros((60 - length(all_MI_source_mina_post{chnl})),length(all_MI_source_mina_post{chnl}((60 - length(all_MI_source_mina_post{chnl})),:)))];    
    end

end

cd(['XXX\subj',num2str(s)])
save (['MI_source_pre_subj_', num2str(s)], 'all_MI_source_mina_pre','miss_movies'); 
save (['MI_source_post_subj_', num2str(s)], 'all_MI_source_mina_post','miss_movies'); 

clearvars -except subjects
    
end


%% Realign the MI_source on FreqPeak for each participant (Freqpeak = frequency in which the MI peaks between lips and audio envelope signals, see stimuli analyses);
clearvars; clc; cd XXX\; 
load ('stim_freqpeak.mat'); %this matrix contains the list of the 60 videos and their frequency of MI peak between lips and sound;

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

%Save realigned MI_source info on MI freqpeak (in Movie condition);
cd(['XXX\subj',num2str(s)])
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
cd(['XXX\subj',num2str(s)])
load (['MI_source_pre_subj_', num2str(s), '_peak']); load (['MI_source_post_subj_', num2str(s), '_peak']); load (['subj', num2str(s), '_source_movie']);
cd (['XXX\headmodel_mat\subj',num2str(s)])
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

%%