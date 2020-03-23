%First, check your triggers in the raw data (nothing is applied in this first step, just control for correctness);
clearvars;clc; addpath XXX\fieldtrip; 
cd XXX\;

%suject ID;
subjects = 2;

s = subjects;
raw_dir = ['XXX\subj_', num2str(s)];
cd(raw_dir);

%check eventtype and eventvalue;
cfg = [];
cfg.dataset = ['subj_',num2str(s),'.bdf'];
cfg.trialdef.eventtype = '?';
temp = ft_definetrial(cfg);
type = {nan(length(temp.event),1)};
trig = {nan(length(temp.event),1)}; 

for t = 1:numel(temp.event) 
    
    if ~isempty(temp.event(1,t).value)
        type{t,1}=temp.event(1,t).type;
        trig{t,1}=temp.event(1,t).value; 
    end
    
end

%% filter continuous data;
clearvars;clc; addpath XXX\fieldtrip; 

%subject ID;
subjects = [2];  


for s = subjects
    
    raw_dir = ['XXX\subj_', num2str(s)];
    cd(raw_dir);

    %filter
    cfg = [];
    cfg.dataset = ['subj_',num2str(s),'.bdf']; 
    cfg.channel = 1:126; 
    
    %band pass filter [1 100] Hz with FIR
    cfg.bpfilter = 'yes';
    cfg.bpfreq = [1 100];
    cfg.bpfilttype = 'fir';
    %band stop filters
    cfg.bsfilter = 'yes';
    cfg.bsfreq = [48 52; 98 102;]; 
    
    %preprocessing
    fil_data = ft_preprocessing(cfg);
    
    %save preprocessed data
    cd('XXX\filtered\')
    save(['subj_',num2str(s),'_bpbsfiltered'], 'fil_data');
    
end


%% epoch the filtered data;
clearvars;clc; addpath XXX\fieldtrip; 

% %%%--TRIGGER's NAME--S%%%
% 65312 = practice; 
% 65282 = multimodal (120); 
% 65288 = unimodal visual (60); 
% 65284 = unimodal auditory (60);

%Subject ID;
subjects = [2];  


for s = subjects
    
    raw_dir = ['XXX\subj_', num2str(s)];
    cd(raw_dir);
    
    %define trials depending on eventtype and eventvalue
    cfg = [];
    cfg.dataset = ['subj_',num2str(s),'.bdf'];
    cfg.trialdef.eventtype = 'STATUS'; 
    cfg.trialdef.eventvalue = [65288;65284]; 
    cfg.trialdef.prestim = 2; 
    cfg.trialdef.poststim = 7;
    cfg = ft_definetrial(cfg); 
    
    cd('XXX\filtered') 
    load(['subj_',num2str(s),'_bpbsfiltered'], 'fil_data');
    
    %epoching the continuous data;
    epoch_data = ft_redefinetrial(cfg,fil_data);
    cd('XXX\unimodal')
    mkdir(['subj_',num2str(s)]);
    cd(['XXX\unimodal\subj_',num2str(s)])
    save(['subj_',num2str(s),'_epoch'], 'epoch_data');

    % %In case something happened during the experiment and we need to remove some specific trials (check it in the first section if the number of trigger is not exactly 120/60/60):
    %     epoch_data_temp = epoch_data;
    %     epoch_data_temp.trialinfo([1:1],:) = []; % select the trials'rows to be removed;
    %     epoch_data_temp.sampleinfo([1:1],:) = []; % select the trials'rows to be removed;
    %     epoch_data_temp.time(:,[1:1]) = []; % select the trials'columns to be removed;
    %     epoch_data_temp.trial(:,[1:1]) = []; % select the trials'columns to be removed;
    %     epoch_data = epoch_data_temp;
    %     clear epoch_data_temp;
    %     save(['subj_',num2str(s),'_epoch'], 'epoch_data');
    
    %downsample the raw data to 512 Hz to speed up the analyses;
    cfg = [];
    cfg.resamplefs = 512;
    cfg.detrend = 'no';
    res_data = ft_resampledata(cfg,epoch_data);
    
    %Trial Info;
    cd(raw_dir)
    soundfile = dataset('File', ['subj' num2str(s) '-Sound.csv'], 'Delimiter' ,',');
    soundtrialnum = soundfile.TRIAL;
    soundtrialstim_num = soundfile.STIM_NUMB;
    soundfreq_peak = soundfile.FREQPEAK;

    moviefile = dataset('File', ['subj' num2str(s) '-Movie.csv'], 'Delimiter' ,',');
    movietrialnum = moviefile.TRIAL;
    movietrialstim_num = moviefile.STIM_NUMB;
    moviefreq_peak = moviefile.FREQPEAK;
    
    %stick condition code to each epoch;
    triggerInfos=res_data.trialinfo;
    
    for i = 1 : numel(triggerInfos)
        
        switch triggerInfos(i)
            
            case 65284
                res_data.trialstructure{i,1}.type = 'sound';
                res_data.trialstructure{i,1}.trigger = triggerInfos (i,1);
                res_data.trialstructure{i,1}.condition = 1;
                res_data.trialstructure{i,1}.trialnum = soundtrialnum(i-60); % soundtrialnum(i); IF SOUND BLOCK FIRST THEN MOVIE BLOCK ORDER;
                res_data.trialstructure{i,1}.trialstim_num = soundtrialstim_num(i-60); % soundtrialstim_num(i); IF SOUND BLOCK FIRST THEN MOVIE BLOCK ORDER;
                res_data.trialstructure{i,1}.freq_peak = soundfreq_peak(i-60); % soundfreq_peak(i); IF SOUND BLOCK FIRST THEN MOVIE BLOCK ORDER;
                
            case 65288
                
                res_data.trialstructure{i,1}.type = 'movie';
                res_data.trialstructure{i,1}.trigger = triggerInfos (i,1);
                res_data.trialstructure{i,1}.condition = 2;
                res_data.trialstructure{i,1}.trialnum = movietrialnum(i); % movietrialnum(i-60); IF SOUND BLOCK FIRST THEN MOVIE BLOCK ORDER;
                res_data.trialstructure{i,1}.trialstim_num = movietrialstim_num(i); % movietrialstim_num(i-60); IF SOUND BLOCK FIRST THEN MOVIE BLOCK ORDER;
                res_data.trialstructure{i,1}.freq_peak = moviefreq_peak(i); % moviefreq_peak(i-60); IF SOUND BLOCK FIRST THEN MOVIE BLOCK ORDER;
        end
    end
    
    res_data.trialinfo = res_data.trialstructure;
   
    cd(['XXX\unimodal\subj_',num2str(s)])
    save(['subj_',num2str(s),'_resampled'], 'res_data');
    
end

%% Pre-processing of signal;

%subject ID;
subjects = 2; 

s = subjects;
raw_dir = ['XXX\unimodal\subj_', num2str(s)];
cd(raw_dir)
load(['subj_',num2str(s),'_resampled']);

%prepare layout;
cfg = [];
cfg.layout='biosemi128';
lay = ft_prepare_layout(cfg, []);

%visual inspection of the trials;
cfg = [];
cfg.channel = {'all', '-EXG1', '-EXG2', '-EXG3'};
cfg.viewmode='vertical';
% cfg.preproc.lpfilter = 'yes';
% cfg.preproc.lpfreq = 70; 
cfg.plotlabels = 'some'; 
cfg = ft_databrowser(cfg,res_data);

%visual inspection of the trials;
cfg = [];
cfg.viewmode='butterfly';
cfg.preproc.lpfilter = 'yes';
cfg.preproc.lpfreq = 70; 
cfg.channel = {'all', '-EXG1', '-EXG2', '-EXG3'};
cfg.ylim = [-200 200];
cfg = ft_databrowser(cfg,res_data);

%coarse artefacts rejection;
cfg = [];
cfg.method='summary';
cfg.channel = {'all'};
cfg.layout = lay;
rejvis_data = ft_rejectvisual(cfg, res_data);

cd(['XXX\unimodal\subj_',num2str(s)]);
save(['subj_',num2str(s),'_beforeICA'], 'rejvis_data');

%% ICA rejections;

%subject ID;
subjects = 2; 

%save pre-ICA data;
s = subjects;
cfg = [];
cfg.method = 'runica';
comp_data = ft_componentanalysis(cfg,rejvis_data);
save(['subj_',num2str(s),'_ICA_comp'], 'comp_data');
    
%topoplot ICA components;
cfg = [];
fourth=round(numel(comp_data.label)/4);
cfg.component = 1:fourth;
cfg.layout = lay;   
figure;ft_topoplotIC(cfg, comp_data);
cfg.component = fourth+1:fourth+fourth;
figure;ft_topoplotIC(cfg, comp_data);
cfg.component = fourth+fourth+1:fourth+fourth+fourth;    
figure;ft_topoplotIC(cfg, comp_data);
cfg.component = fourth+fourth+fourth+1:numel(comp_data.label);  
figure;ft_topoplotIC(cfg, comp_data);

%visualize the timecourse of indivdual components; 
cfg = []; 
cfg.channel = 1:8; 
cfg.layout = lay; 
cfg.viewmode = 'component'; 
ft_databrowser(cfg, comp_data);

%visualize the frequency spectrum of components;
mainfig = figure(); 
tabgroup = uitabgroup(mainfig, 'Position', [0 .05 1 .95]);
cnt = 0;
ax_ = nan(size(comp_data.label));

for step_ =1 : 10:  size(comp_data.label) 
    
    cnt = cnt+1; vec = step_:step_+9;   
    tab(cnt)=uitab(tabgroup,'Title', ['Comps_', num2str(vec(1)) '-' num2str(vec(end))]);
    axes('parent',tab(cnt))    
    ind_ = 0;
    
    for comp = vec 
        
        if comp> size(comp_data.trial{1,1},1); continue; end
        ind_ = ind_+1;
        
        % ----- here we calculate the fft to get the powerspectrum
        pspctrm = 0;
        
        for tr = 1:numel(comp_data.trial)
            
            signal  = comp_data.trial{1,tr}(comp,:);
            N       = length(signal);
            nyquist = comp_data.fsample/2;
            fourierCoefs = zeros(size(signal)); %what's this for?
            frequencies = linspace(0,nyquist,floor(N/2)+1); %get linear spaced N/2+1 freqs from 0 to nyquist
            fourierCoefsF = fft(signal) / N;
            pspctrm = pspctrm + abs(fourierCoefsF(1:length(frequencies)))*2;
            
        end
        
        pspctrm =pspctrm/numel(comp_data.trial);
        % -----------------------------------------------------------------
                
        subplot(5,4,2*(ind_-1)+1);
        hold on;
        
        cfg = [];
        cfg.component = comp; 
        cfg.layout = 'Biosemi128'; 
        cfg.comment = 'no';
        ft_topoplotIC(cfg, comp_data)
        ax_(ind_+step_-1)=subplot(5,4,2*ind_);
        set(ax_(ind_),'xlim', [0 100]);
        hold on;
        plot(frequencies,pspctrm)
        
    end
    ax_ = ax_(~isnan(ax_));
    hax = @(src, ax_) set(ax_, 'xlim', [0 get( src, 'Value' )] );
    uicontrol('Style', 'slider', 'Units','normalized', 'Min',20,'Max',300,'Value',60, 'Position', [0 0 1 0.05], 'Callback', @(src,evt) hax( src, ax_ ) );
end

%% removing components from your data;

cfg = [];
cfg.component = [1 2 5 11 14]; 

%save info of rejected components
rejcomp_data = ft_rejectcomponent(cfg, comp_data);
rejcomp_data.reject = cfg.component; 
close all;

%% interpolation of bad channels and rereference;

%check a last time if any channel looks bad before interpolarisation;
cfg = [];
cfg.channel = {'all', '-EXG1', '-EXG2', '-EXG3'};
cfg.plotlabels = 'some';
cfg.viewmode = 'vertical';
cfg.preproc.lpfilter = 'yes';
cfg.preproc.lpfreq = 40; % filters are optional, but they help to get a clean look at the data;
cfg.layout = lay;
cfg = ft_databrowser(cfg,rejcomp_data);

cfg = [];
cfg.method = 'summary';
cfg.channel = {'all'};
cfg.layout = lay;
rejcomp_data = ft_rejectvisual(cfg, rejcomp_data);

%Now we interpolarize the missing channels removed; 
allLabels = lay.label(1:128); 
misschans=find(~ismember(allLabels , rejcomp_data.label)); %indices of  the channels that are not in the data 
nmiss=length(misschans); %the amount of missing channels
cleanchans=[find(ismember(allLabels,  rejcomp_data.label));129;130;131]; %indices of labels of the clean channels
nclean=size(allLabels,1) -nmiss; %the amount of clean channels

for t=1:length(rejcomp_data.trial) %for every trial
    tmp1=rejcomp_data.trial{t};  % we write that trial into tmp1
    tmp2=[tmp1;zeros(nmiss,length(tmp1))]; % then we add zeros for every channel that is missing at the end 
    tmp3=[tmp2 [cleanchans;misschans]]; % then we write the indices of the channels in the last column...
    tmp4=sortrows(tmp3,length(tmp3)); %... and sort rows along that column
    rejcomp_data.trial{t}=tmp4(:,1:length(tmp4)-1); % finally we overwrite the trial
end
rejcomp_data.label = [allLabels;'EXG1';'EXG2';'EXG3']; % and we overwrite the labels

%use polhemus data: 
raw_dir = ['XXX\subj_', num2str(s)];
cd(raw_dir);
pos_file = (dir(['subj',num2str(s),'*.pos']));

%neighbouring the electrodes with polhemus: 
elec = ft_read_sens(pos_file.name);  
elec.label(1:128,1) = allLabels;
cfg=[];
cfg.method = 'triangulation';
cfg.elec = elec;
cfg.feedback  = 'yes';
neighbours = ft_prepare_neighbours(cfg); 

%interpolarize;
cfg = [];
%cfg.layout = 'easycapM15';
cfg.badchannel = allLabels(misschans);
cfg.neighbours = neighbours;
interp_data = ft_channelrepair(cfg, rejcomp_data);

%re-reference to a common average;
cfg = [];
cfg.reref = 'yes';
cfg.refchannel = 1:128; 
reref_data = ft_preprocessing(cfg,interp_data);
close all;

%% manually reject trials with artefacts, last check and save preprocessed data for analyses;

cfg = [];
cfg.viewmode = 'vertical';
cfg.preproc.lpfilter = 'yes';
cfg.preproc.lpfreq = 30; 
cfg.layout = lay;
cfg = ft_databrowser(cfg,reref_data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---- ONLY IF YOU FIND ARTIFACTS IN SOME TRIALS AND MARKED THEM MANUALLY
% BEFORE IN THE DATABROWSER ---------------------------------------------

% identfiy artifact trials and select them;
% artifacts = cfg.artfctdef.visual;
% for n = 1:size(artifacts.artifact,1)
%     art_samples = artifacts.artifact(n,1) : artifacts.artifact(n,2);
%     contains_art = any(bsxfun(@le, reref_data.sampleinfo(:,1), art_samples),2)...
%             & any(bsxfun(@ge, reref_data.sampleinfo(:,2), art_samples),2);    
%     artrl(n,1) = find(contains_art);
% end
% 
% % store the trials with artifacts;
% artrl = unique(artrl);
% trials = 1:numel(reref_data.trial);
% trialsel = find(~ismember(trials,artrl));

%kick out the artifacts;
%cfg = [];
%cfg.trials=trialsel;
%preproc_data = ft_preprocessing(cfg,reref_data);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Finalize he preprocessing with or without artifact before saving;
cfg = [];
preproc_data = ft_preprocessing(cfg,reref_data);

% cfg = [];
% cfg.viewmode = 'vertical';
% cfg.preproc.lpfilter = 'yes';
% cfg.preproc.lpfreq = 70; 
% cfg.layout = lay;
% cfg = ft_databrowser(cfg,preproc_data);

%SAVE YOUR CLEAN DATA, READY FOR ANALYSES;
cd(['XXX\unimodal\subj_',num2str(s)])
save(['subj_',num2str(subjects),'_after_preproc'], 'preproc_data');

%%