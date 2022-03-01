%Calculate the phase information for every video stimulus (ONCE FOR ALL);
clearvars;clc; 
addpath XXX\fieldtrip\;
addpath XXX\Circular_Statistics_Toolbox\;
cd XXX\;
load Stim_signal; load data_template; 

%structure signal of each video in Fieldtrip structure;
data_mina = data_template; 
data_mina.label = {'Video'};
data_mina.fsample = 250;

for i=1:length(dat_mina(:,1))  
    
    data_mina.time{i,1} = 0:1/data_mina.fsample:5-1/data_mina.fsample;
    data_mina.trial{i,1} = dat_mina(i,:);

end

data_mina.trialinfo = [1:60]';

for i=1:length(data_mina.trialinfo)
    
        temp.trialinfo{i,1}.type = 'video';
        temp.trialinfo{i,1}.video_name = Stim_infos.VIDEO_NAME(i);
        temp.trialinfo{i,1}.name_in_task = Stim_infos.NAME_in_task(i);
        temp.trialinfo{i,1}.freqpeak = Stim_infos.FREQPEAK(i);

end

data_mina.trialinfo = temp.trialinfo;
clear temp.trialinfo

cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [2 6];        
cfg.bpfilttype = 'fir';
cfg.hilbert = 'angle';
cfg.trials = find(cell2mat(cellfun(@(x) x.freqpeak == 4, data_mina.trialinfo, 'UniformOutput', false)));
temp_phase_4 = ft_preprocessing(cfg,data_mina);

cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [3 7];        
cfg.bpfilttype = 'fir';
cfg.hilbert = 'angle';
cfg.trials = find(cell2mat(cellfun(@(x) x.freqpeak == 5, data_mina.trialinfo, 'UniformOutput', false)));
temp_phase_5 = ft_preprocessing(cfg,data_mina);

cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [4 8];        
cfg.bpfilttype = 'fir';
cfg.hilbert = 'angle';
cfg.trials = find(cell2mat(cellfun(@(x) x.freqpeak == 6, data_mina.trialinfo, 'UniformOutput', false)));
temp_phase_6 = ft_preprocessing(cfg,data_mina);

cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [5 9];        
cfg.bpfilttype = 'fir';
cfg.hilbert = 'angle';
cfg.trials = find(cell2mat(cellfun(@(x) x.freqpeak == 7, data_mina.trialinfo, 'UniformOutput', false)));
temp_phase_7 = ft_preprocessing(cfg,data_mina);

cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [6 10];        
cfg.bpfilttype = 'fir';
cfg.hilbert = 'angle';
cfg.trials = find(cell2mat(cellfun(@(x) x.freqpeak == 8, data_mina.trialinfo, 'UniformOutput', false)));
temp_phase_8 = ft_preprocessing(cfg,data_mina);

%structure again;
Phase_mina.fsample = 250;
Phase_mina.label = data_mina.label; 
Phase_mina.time = [temp_phase_4.time';temp_phase_5.time';temp_phase_6.time';temp_phase_7.time';temp_phase_8.time']';
Phase_mina.trialinfo = [temp_phase_4.trialinfo;temp_phase_5.trialinfo;temp_phase_6.trialinfo;temp_phase_7.trialinfo;temp_phase_8.trialinfo];
Phase_mina.trial = [temp_phase_4.trial';temp_phase_5.trial';temp_phase_6.trial';temp_phase_7.trial';temp_phase_8.trial']';

clear temp_phase_4 temp_phase_5 temp_phase_6 temp_phase_7 temp_phase_8 temp

%save everything alltogether;
cd XXX\Signal_Infos\;
save Phase_mina Phase_mina

%%
