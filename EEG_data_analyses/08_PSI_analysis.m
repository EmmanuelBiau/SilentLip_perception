clearvars; clc; addpath xxx;
addpath xxx\fieldtrip; ft_defaults

subjects = [2:9 11:17 19:26];
PSI_VA = []; psi_spectrum_T1 =[];

scount = 0;

for s = subjects
    
    scount = scount + 1;
    cd (['xxx\subj',num2str(s)]);
    load(['visual_entrainment_source_realigned_subj' num2str(s)]);

    %create data structure containing the 2 source signals;
    data_temp = source_left_audio_realigned;

    %get trial info;
    trlinfo = zeros(numel(data_temp.trialinfo),1);
    for j = 1 : numel(trlinfo)
        trlinfo(j,1) = data_temp.trialinfo{j}.condition;
        trlinfo(j,2) = data_temp.trialinfo{j}.freq_peak;
    end

    trial_temp = [];
    for ttp = 1: size(data_temp.trial,2)
        trial_temp{1,ttp} = cat(1,source_left_audio_realigned.trial{1,ttp},source_left_visual_realigned.trial{1,ttp});
    end

    %sourceAL = source in the left auditory cortex; sourceVL = source in the left visual cortex;
    data_temp.trial = trial_temp;
    data_temp.label = {'sourceAL', 'sourceVL'};

    %select the two time-windows of interest (early T1 and late T2);
    cfg = [];
    cfg.latency = [0.5 2.5];
    cfg.trials = find(cell2mat(cellfun(@(x) x.condition ==2, data_temp.trialinfo, 'UniformOutput', false)));
    data_T1 = ft_selectdata(cfg,data_temp);
    cfg = [];
    cfg.latency = [2.5 4.5];
    cfg.trials = find(cell2mat(cellfun(@(x) x.condition ==2, data_temp.trialinfo, 'UniformOutput', false)));
    data_T2 = ft_selectdata(cfg,data_temp);

    %--------------------------------------------------------------------------%
    %     Calculate the Phase Slope Index and Granger causality spectra;
    %--------------------------------------------------------------------------%

    %fit multivariate autoregressive model on data;
    cfg = [];
    cfg.order = 5;
    cfg.toolbox = 'bsmart';
    cfg.zscore = 'yes';
    mdata_T1  = ft_mvaranalysis(cfg, data_T1);
    mdata_T2  = ft_mvaranalysis(cfg, data_T2);

    %compute the spectral transfer matrix;
    %we dont use it for the psi;
    cfg=[];
    cfg.method = 'mvar';
    cfg.zscore = 'yes';
    mfreq_T1 = ft_freqanalysis(cfg,mdata_T1);
    mfreq_T2 = ft_freqanalysis(cfg,mdata_T2);

    %compute the fourier transform and realign the spectrum realigned on the freq_peak +/- 4Hz;
    cfg = [];
    cfg.method = 'mtmfft';
    cfg.taper = 'dpss';
    cfg.output = 'fourier';
    cfg.foi = 0:1:30;
    cfg.tapsmofrq = 2;

    for fp = 4:8

        if fp == 4

            cfg.trials = find(cell2mat(cellfun(@(x) x.freq_peak == fp, data_T1.trialinfo, 'UniformOutput', false)));
            freq_temp1 = ft_freqanalysis(cfg, data_T1);freq_temp2 = ft_freqanalysis(cfg, data_T2);
            freq_temp1.freq =freq_temp1.freq(1:end-4);freq_temp2.freq =freq_temp2.freq(1:end-4);
            freq_temp1.fourierspctrm = freq_temp1.fourierspctrm(:,:,1:end-4);freq_temp2.fourierspctrm = freq_temp2.fourierspctrm(:,:,1:end-4);
            freq_T1_4hz = freq_temp1;freq_T2_4hz = freq_temp2;
            clear freq_temp1 freq_temp2;

        elseif fp == 5

            cfg.trials = find(cell2mat(cellfun(@(x) x.freq_peak == fp, data_T1.trialinfo, 'UniformOutput', false)));
            freq_temp1 = ft_freqanalysis(cfg, data_T1);freq_temp2 = ft_freqanalysis(cfg, data_T2);
            freq_temp1.freq =freq_temp1.freq(2:end-3);freq_temp2.freq =freq_temp2.freq(2:end-3);
            freq_temp1.fourierspctrm = freq_temp1.fourierspctrm(:,:,2:end-3);freq_temp2.fourierspctrm = freq_temp2.fourierspctrm(:,:,2:end-3);
            freq_T1_5hz = freq_temp1;freq_T2_5hz = freq_temp2;
            clear freq_temp1 freq_temp2;

        elseif fp == 6

            cfg.trials = find(cell2mat(cellfun(@(x) x.freq_peak == fp, data_T1.trialinfo, 'UniformOutput', false)));
            freq_temp1 = ft_freqanalysis(cfg, data_T1);freq_temp2 = ft_freqanalysis(cfg, data_T2);
            freq_temp1.freq =freq_temp1.freq(3:end-2);freq_temp2.freq =freq_temp2.freq(3:end-2);
            freq_temp1.fourierspctrm = freq_temp1.fourierspctrm(:,:,3:end-2);freq_temp2.fourierspctrm = freq_temp2.fourierspctrm(:,:,3:end-2);
            freq_T1_6hz = freq_temp1;freq_T2_6hz = freq_temp2;
            clear freq_temp1 freq_temp2;

        elseif fp == 7

            cfg.trials = find(cell2mat(cellfun(@(x) x.freq_peak == fp, data_T1.trialinfo, 'UniformOutput', false)));
            freq_temp1 = ft_freqanalysis(cfg, data_T1);freq_temp2 = ft_freqanalysis(cfg, data_T2);
            freq_temp1.freq =freq_temp1.freq(4:end-1);freq_temp2.freq =freq_temp2.freq(4:end-1);
            freq_temp1.fourierspctrm = freq_temp1.fourierspctrm(:,:,4:end-1);freq_temp2.fourierspctrm = freq_temp2.fourierspctrm(:,:,4:end-1);
            freq_T1_7hz = freq_temp1;freq_T2_7hz = freq_temp2;
            clear freq_temp1 freq_temp2;

        elseif fp == 8

            cfg.trials = find(cell2mat(cellfun(@(x) x.freq_peak == fp, data_T1.trialinfo, 'UniformOutput', false)));
            freq_temp1 = ft_freqanalysis(cfg, data_T1);freq_temp2 = ft_freqanalysis(cfg, data_T2);
            freq_temp1.freq =freq_temp1.freq(5:end);freq_temp2.freq =freq_temp2.freq(5:end);
            freq_temp1.fourierspctrm = freq_temp1.fourierspctrm(:,:,5:end);freq_temp2.fourierspctrm = freq_temp2.fourierspctrm(:,:,5:end);
            freq_T1_8hz = freq_temp1;freq_T2_8hz = freq_temp2;
            clear freq_temp1 freq_temp2;

        end

    end

    %now concatenate all the trials in the same data T1 and T2; Now the freq #5 corresponds to the freq_peak of each trial;
    freq_T1 = freq_T1_4hz;
    freq_T1.fourierspctrm = cat(1,freq_T1_4hz.fourierspctrm,freq_T1_5hz.fourierspctrm,freq_T1_6hz.fourierspctrm,freq_T1_7hz.fourierspctrm,freq_T1_8hz.fourierspctrm);
    freq_T1.cumsumcnt = cat(1,freq_T1_4hz.cumsumcnt,freq_T1_5hz.cumsumcnt,freq_T1_6hz.cumsumcnt,freq_T1_7hz.cumsumcnt,freq_T1_8hz.cumsumcnt);
    freq_T1.cumtapcnt = cat(1,freq_T1_4hz.cumtapcnt,freq_T1_5hz.cumtapcnt,freq_T1_6hz.cumtapcnt,freq_T1_7hz.cumtapcnt,freq_T1_8hz.cumtapcnt);
    freq_T1.trialinfo = cat(1,freq_T1_4hz.trialinfo,freq_T1_5hz.trialinfo,freq_T1_6hz.trialinfo,freq_T1_7hz.trialinfo,freq_T1_8hz.trialinfo);
    freq_T1 = rmfield(freq_T1,'trialinfo');
    clear freq_T1_4hz freq_T1_5hz freq_T1_6hz freq_T1_7hz freq_T1_8hz
    freq_T2 = freq_T2_4hz;
    freq_T2.fourierspctrm = cat(1,freq_T2_4hz.fourierspctrm,freq_T2_5hz.fourierspctrm,freq_T2_6hz.fourierspctrm,freq_T2_7hz.fourierspctrm,freq_T2_8hz.fourierspctrm);
    freq_T2.cumsumcnt = cat(1,freq_T2_4hz.cumsumcnt,freq_T2_5hz.cumsumcnt,freq_T2_6hz.cumsumcnt,freq_T2_7hz.cumsumcnt,freq_T2_8hz.cumsumcnt);
    freq_T2.cumtapcnt = cat(1,freq_T2_4hz.cumtapcnt,freq_T2_5hz.cumtapcnt,freq_T2_6hz.cumtapcnt,freq_T2_7hz.cumtapcnt,freq_T2_8hz.cumtapcnt);
    freq_T2.trialinfo = cat(1,freq_T2_4hz.trialinfo,freq_T2_5hz.trialinfo,freq_T2_6hz.trialinfo,freq_T2_7hz.trialinfo,freq_T2_8hz.trialinfo);
    freq_T2 = rmfield(freq_T2,'trialinfo');
    clear freq_T2_4hz freq_T2_5hz freq_T2_6hz freq_T2_7hz freq_T2_8hz

    %calculate the PSI for the two time-windows on realigned spectrum;
    cfg = [];
    cfg.method = 'psi';
    cfg.bandwidth = 2;
    cf.zscore = 'yes';
    psi_T1_temp = ft_connectivityanalysis(cfg, freq_T1);
    psi_T2_temp = ft_connectivityanalysis(cfg, freq_T2);

    %select the spectrum of interest freq.peak +/- 3 bins;
    %Now freqpeak bin will be the #4;
    cfg = [];
    cfg.frequency = [2 8];
    cfg.avgoverfreq = 'no';
    psi_T1 = ft_selectdata(cfg,psi_T1_temp);
    psi_T2 = ft_selectdata(cfg,psi_T2_temp);
    psi_T1.freq = round(psi_T1.freq);
    psi_T2.freq = round(psi_T2.freq);

    %Store the PSI from VL source to AL source direction (i.e. if PSI < 0, VL leads; if PSI > 0, AL leads);
    psi_spectrum_T1(scount,:) = psi_T1.psispctrm(2,1,:);
    psi_spectrum_T2(scount,:) = psi_T2.psispctrm(2,1,:);
    PSI_VA(scount,1) = psi_T1.psispctrm(2,1,4); % select freqpeak bin;
    PSI_VA(scount,2) = psi_T2.psispctrm(2,1,4); % select freqpeak bin;

    clearvars -except scount PSI_VA psi_spectrum_T1 psi_spectrum_T2

end

%% Assess difference of PSI against 0 in the T1 and T2 windows;

mean_PSI(1,1) = mean(PSI_VA(:,1));
mean_PSI(2,2) = mean(PSI_VA(:,2));
mean_PSI(2,1) = std(PSI_VA(:,1));
mean_PSI(2,2) = std(PSI_VA(:,2));
[h_T1vs0,pvalue_T1vs0,~,~] = ttest(PSI_VA(:,1),0, 'tail','left');
[h_T2vs0,pvalue_T2vs0,~,~] = ttest(PSI_VA(:,2),0, 'tail','left');

%Plot PSI in T1 and T2 windows;
dat = PSI_VA;
col = [0 0 0];
dotsize = 15;
smooth = 0;

close all;
vert_rain_plot(dat,col,dotsize,smooth)
pause(0.01);
ax = gca;
ax.XTick = [1 2];
ax.XTickLabel = {'T1 window', 'T2 window'};
ylabel('PSI');

%%