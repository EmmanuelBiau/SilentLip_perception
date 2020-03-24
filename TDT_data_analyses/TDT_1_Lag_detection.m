%%Calculate the lag between the first frame of the movie and the sound onset for every trial and store it for later; 
clearvars;clc;

%subjects ID;
subjects = 1:29;



for s = subjects
    
    %load the subject logfile;    
    cd (['XXX\subj' num2str(s)]);
    load(['subj',num2str(s),'-Tone_detectionTask']);

    %Import data;
    cfg = [];
    cfg.dataset = ['subj',num2str(s),'.bdf']; 
    cfg.channel = {'A1','Ana1','Ana2','Ana3'}; %only select the scalp channel plus 3 external eye channels;
    raw_data = ft_preprocessing(cfg);

    %define trials depending on eventtype and eventvalue
    cfg = [];
    cfg.dataset = ['subj',num2str(s),'.bdf']; 
    cfg.trialdef.eventtype = 'STATUS'; 
    cfg.trialdef.eventvalue = 65282; 
    cfg.trialdef.prestim = 1; 
    cfg.trialdef.poststim = 8; 
    cfg = ft_definetrial(cfg); 
    TdT_triggers = ft_redefinetrial(cfg,raw_data);
    TdT_triggers.label = {'A1','Visual','Audio1','Audio2'};

    %stick condition code to each epoch;
for i = 1 : length(TdT_triggers.trialinfo)
  
    if respMat.FREQPEAK{i} == 0
        res_data_temp.trialinfo{i,1}.type = 'trial';
        res_data_temp.trialinfo{i,1}.trigger = TdT_triggers.trialinfo(i,1);
        res_data_temp.trialinfo{i,1}.video = str2num(cell2mat(extractBetween(respMat.VIDEO{i},"_","_")));
        res_data_temp.trialinfo{i,1}.trialnum = i; 
        res_data_temp.trialinfo{i,1}.freq_peak = respMat.FREQPEAK{i};
        
    elseif respMat.FREQPEAK{i} ~= 0          
        res_data_temp.trialinfo{i,1}.type = 'video';
        res_data_temp.trialinfo{i,1}.trigger = TdT_triggers.trialinfo(i,1);
        res_data_temp.trialinfo{i,1}.video = str2num(cell2mat(extractBetween(respMat.VIDEO{i},"_","_")));
        res_data_temp.trialinfo{i,1}.trialnum = i; 
        res_data_temp.trialinfo{i,1}.freq_peak = respMat.FREQPEAK{i};
        
    end
    
end

    TdT_triggers.trialinfo = res_data_temp.trialinfo;
    clear res_data_temp

    %plot the signals;
    %cfg = [];
    %cfg.layout = 'biosemi128';
    %ft_databrowser(cfg, TdT_triggers);

    %create a fieldtrip structure for original soundwaves;
    load('XXX\whitenoise.mat'); 
    wav = wnseq;
    fs = 44100;

for i = 1:numel(TdT_triggers.trial)
    
    sound_wav.trial{1,i} = wav;
    sound_wav.trial{1,i}(2,:) = wav;
    sound_wav.time{1,i} = 0:1/fs:(length(sound_wav.trial{1,i})-1)/fs;
    sound_wav.trialinfo{i,1} = TdT_triggers.trialinfo{i,1}; 
    
end

    sound_wav.fsample = fs;
    sound_wav.label{1,1} = 'Audio1';
    sound_wav.label{2,1} = 'Audio2';
    sound_wav.sampleinfo(:,1) = 1:length(sound_wav.trial{1,i}):length(sound_wav.trial{1,i})*numel(TdT_triggers.trial);
    sound_wav.sampleinfo(:,2) = length(sound_wav.trial{1,i}):length(sound_wav.trial{1,i}):length(sound_wav.trial{1,i})*numel(TdT_triggers.trial)+length(sound_wav.trial{1,i})-1;

    %downsample the sounds to 2048;
    cfg = [];
    cfg.resamplefs = 2048;
    cfg.detrend = 'no';
    presented_sound = ft_resampledata(cfg,sound_wav);

    % %plot;
    %cfg = [];
    %cfg.layout = 'biosemi128';
    %ft_databrowser(cfg, res_sound);

    %select only 5 s length signals;
    cfg = [];
    cfg.channel = 'Audio2';
    cfg.latency = [0 5];
    recorded_signal = ft_selectdata(cfg, TdT_triggers);

    %cross correlation;
    sound_onset = zeros(numel(recorded_signal.trial),1);

for i = 1:numel(presented_sound.trial)
    
    sound_ch1 = presented_sound.trial{1,i}(1,:); % presented sound of the trial;
    sig_ch1 = recorded_signal.trial{1,i}(1,1:length(sound_ch1)); % recorded analogic signal of the presented sound in the trial;
    
    %calculate the envelope of the recorded and presented signals;
    recorded_sig_temp = (sig_ch1 - mean(sig_ch1))./std(sig_ch1); % normalized recorded signal;
    [up_sig,lo_sig] = envelope(recorded_sig_temp);
    enve_sig = up_sig-lo_sig;
    
    presented_sound_temp = (sound_ch1 - mean(sound_ch1))./std(sound_ch1);
    [up_sound,lo_sound] = envelope(presented_sound_temp);
    enve_sound = up_sound - lo_sound;
    
    %crosscorrelations between envelope of the presented sound and recorded signal;
    [c, lags] = xcorr(enve_sig,enve_sound); % c = crosscorrelation coefficients of each lag; lag = distances between the 2 shifted time-series;
    pos_lags = find(lags >=0); %keep lags for when envelop_sound leading only;
    max_c = max(c(pos_lags)); %find the greater cc coefficient when envelop_sound is leading;
    sound_onset(i,1) = recorded_signal.time{1,i}(lags(c == max_c)+1); % find the lag between envelop_sound and envelop_signal for which the coefficent is the greatest;
    
    %figure;
    %plot([sig_temp(lags(c == max_c)+1:end);sound_temp(1:end-lags(c == max_c))]');
    %title(['trial ',num2str(i)]);
    
    clear sound_ch* sig_ch* sig_temp sound_temp

end

    %Detect the video onset for each multimodal trial;
    video_onset = zeros(numel(TdT_triggers.trial),1);

for i = 1:numel(TdT_triggers.trial)
    
    [~, temp_pos] = findpeaks(TdT_triggers.trial{i}(2,:),'MinPeakHeight',3000);
    video_onset(i,1) = TdT_triggers.time{i}(temp_pos(1));

end

    %Calculate the lag between video and sound onset and add variable in the respMat matrix;
    video_sound_lag = sound_onset - video_onset;
    VA_LAG = video_sound_lag(:,1);
    respMat = addvars(respMat,VA_LAG,'After','SCREEN_OFFSET');
    clear VA_LAG
    
    %save the matrix containing the VA onset lag of the participant;
    cd(['XXX\subj' num2str(s)]);
    save(['subj' num2str(s) '-Tone_detectionTask_lags'],'respMat','participant_parameters','alltrials');

    keep subjects
    
end


%% MAKE SOME PLOTS OF THE LAG;
figure('Name','VA LAG - Histo');
histogram(round(respMat.VA_LAG*1000));
xlabel('Ordered trials');
ylabel('VA Lag (ms)');
title('VA LAG - Histo');


figure('Name','VA LAG - per trial');
plot(respMat.VA_LAG);
xlim([1 length(respMat.VA_LAG)]);
ylim([min(respMat.VA_LAG)-0.001 max(respMat.VA_LAG)+0.001]);
xlabel('Ordered trials');
ylabel('VA Lag (ms)');
title('VA LAG - per trial');

%%