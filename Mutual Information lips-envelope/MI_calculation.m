%% 1 - Upsample Lips information to 250Hz for later analyses;
clearvars;close;clc;

%Configure your path;
cd XXX\stat_videos\;
out = 'XXX\Lips_envelope_signals\';
Folder = cd; Folder = fullfile(Folder);

%Prepare the list of the stimuli to analyses; 
list_v = dir('Stat_video*.mat');
clear list_v_num

for a = 1:length(list_v)  
    list_v_num(a) = str2num(list_v(a).name(11:end-4));
end

%load visual and audio signals from the lists;
[newnum, ~] = sort(list_v_num);
  
for s = newnum

    %Load valid videos (with good lip detection);
    load(fullfile(Folder,['Stat_video', num2str(s)]), 'Stat');

    %Resample area, major and minor axe information;
    fps = 25;
    ds_area = resample(Stat(1,:),250,fps); 
    ds_maja = resample(Stat(2,:),250,fps); 
    ds_mina = resample(Stat(3,:),250,fps); 

    %save all lip signals together (area,mina and maja dimensons); 
    save([out,'Info_lips_', num2str(s)], 'ds_area','ds_maja','ds_mina');
    
end

%% 2- Extract Audio Speech Envelope and downsample to 250Hz;
clearvars;close;clc;

%Configure your path;
cd XXX\audio\; 
out = 'XXX\Lips_envelope_signals\';
addpath('XXX\ChimeraSoftware'); 
Folder = cd; Folder = fullfile(Folder);

%Prepare the list of the stimuli to analyses; 
list_v = dir('audio_*.wav');
clear list_v_num

for a = 1:length(list_v)  
    list_v_num(a) = str2num(list_v(a).name(7:end-4));
end

%load visual and audio signals from the lists;
[newnum, ~] = sort(list_v_num);

for s = newnum

    %Import audio, Band pass filter, Hilbert transform, resampling audio signal;
    samples = [1,5*44100];
    [ audio, fs ] = audioread(['audio_', num2str(s),'.wav'], samples); 
    signal  = audio(:,1);
    
    %Constructs nine frequency bands in the range 100–10,000 Hz;  
    fco = equal_xbm_bands(100,10000,8); 

    for k = 1:length(fco)-1

        if (k==1)
            [b,a]=butter(3,2*[fco(k) fco(k+1)]/fs); % Create the vector of the 9 bands cutoff frequencies in Hz;
        else
            [b,a]=butter(3,2*[fco(k) fco(k+1)]/fs); % Apply a fourth-order band pass Butterworth filter on the auditory signal;
        end

        tmp = filtfilt(b,a,signal);                 % Zero-phase digital filtering;
        tmp2 = abs(hilbert(tmp));                   % Amplitude envelopes for each band are computed as absolute values of the Hilbert transform;
        env_decomp(k,:) = resample(tmp2,250,fs);    % Resample the envelopes at 250Hz to be later matched with visual information;
        
    end
    
    %Average across bands to obtain a wide-band amplitude envelope for all further analysis;
    ds_audio = mean(env_decomp,1);                  
    save([out, 'Info_audio_', num2str(s)], 'ds_audio','env_decomp');

    clear tmp tmp2 k env_decomp

end

%% Now calculate the Mutual Informnation between Lips and Envelope for each video;
clear; clc;
addpath XXX\gauss_info\; addpath XXX\gauss_info\mex\;
addpath XXX\fieldtrip\;ft_defaults 

%Prepare the list of the stimuli to analyses; 
cd XXX\Lips_envelope_signals\;
Folder = cd; Folder = fullfile(Folder);

list_v = dir('Info_lips_*.mat');
list_a = dir('Info_audio_*.mat');
clear list_v_num

for a = 1:length(list_v)  
    list_v_num(a) = str2num(list_v(a).name(11:end-4));
end

%load visual and audio signals from the lists;
[newnum, inx] = sort(list_v_num);

for a = 1:length(inx)
    load(list_v(inx(a)).name,'ds_area'); % Here I only tried area info, but you could try minor as well.
    dat_area(a,:) = ds_area; clear ds_area

    load(list_v(inx(a)).name,'ds_maja');
    dat_maja(a,:) = ds_maja; clear ds_maja

    load(list_v(inx(a)).name,'ds_mina');
    dat_mina(a,:) = ds_mina; clear ds_mina

    load(list_a(inx(a)).name);
    dat_audio(a,:) = ds_audio; clear ds_audio
end

%Structure the V and A signals to fieldtrip structure;
load data_template;
dat1 = []; dat1 = data_template;
dat2 = []; dat2 = data_template;
dat3 = []; dat3 = data_template;
dat4 = []; dat4 = data_template;
%fill with the ordered lips and audio signals;
dat1.trial{1,1} = dat_area;
dat2.trial{1,1} = dat_audio;
dat3.trial{1,1} = dat_maja;
dat4.trial{1,1} = dat_mina;
%give time info (5s at 250Hz = 1250 time-points);
dat1.time{1,1} = data_template.time{1,1}(1:length(dat_area(1,:)));
dat2.time{1,1} = data_template.time{1,1}(1:length(dat_audio(1,:)));
dat3.time{1,1} = data_template.time{1,1}(1:length(dat_maja(1,:)));
dat4.time{1,1} = data_template.time{1,1}(1:length(dat_mina(1,:)));
    
%give fake labels; 
label_temp = cell(length(newnum),1);

for m = 1:length(label_temp)
    label_temp(m) = {['A',num2str(m)]};
end

dat1.label = label_temp; 
dat2.label = label_temp;
dat3.label = label_temp;
dat4.label = label_temp;
clear label_temp m;

for f = 1:20
    
    %bandpass filter;

    if f < 3
        cfg             = [];
        cfg.lpfilter    = 'yes';
        cfg.lpfreq      = f;
    else
        cfg             = [];
        cfg.bpfilter    = 'yes';
        cfg.bpfreq      = [f-2 f+2];
    end
    
    cfg.hilbert = 'complex';
    dat5 = []; dat5 = ft_preprocessing(cfg, dat1);
    dat6 = []; dat6 = ft_preprocessing(cfg, dat2);
    dat7 = []; dat7 = ft_preprocessing(cfg, dat3);
    dat8 = []; dat8 = ft_preprocessing(cfg, dat4);
       
    %MI between V and A
    for c = 1:length(dat6.label)
        
        dat_area  = []; dat_area  = dat5.trial{1,1}(c,:);
        dat_audio = []; dat_audio = dat6.trial{1,1}(c,:);
        dat_maja  = []; dat_maja  = dat7.trial{1,1}(c,:);
        dat_mina  = []; dat_mina  = dat8.trial{1,1}(c,:);
        
        dat_area  = []; dat_area  = dat5.trial{1,1}(c,:);
        dat_audio = []; dat_audio = dat6.trial{1,1}(c,:);
        dat_maja  = []; dat_maja  = dat7.trial{1,1}(c,:);
        dat_mina  = []; dat_mina  = dat8.trial{1,1}(c,:);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Phase information;
        dat_area  = dat_area./abs(dat_area);
        dat_audio = dat_audio./abs(dat_audio);
        dat_maja  = dat_maja./abs(dat_maja);
        dat_mina  = dat_mina./abs(dat_mina);
        
        %Copula normalization;
        cdat_area  = []; cdat_area  = copnorm([real(dat_area)' imag(dat_area)']);
        cdat_audio = []; cdat_audio = copnorm([real(dat_audio)' imag(dat_audio)']);
        cdat_maja  = []; cdat_maja  = copnorm([real(dat_maja)' imag(dat_maja)']);
        cdat_mina  = []; cdat_mina  = copnorm([real(dat_mina)' imag(dat_mina)']);
        
        %Calculate Gaussian- entropy calculation and likelihood test (to describe better);     
        mi_area_original(c,f) = info_gg(cdat_area, cdat_audio, true, false, false);
        mi_maja_original(c,f) = info_gg(cdat_maja, cdat_audio, true, false, false);
        mi_mina_original(c,f) = info_gg(cdat_mina, cdat_audio, true, false, false);
    
    end
    
end
        
%Sub-zeros values to approx. zeros;
mi_area_original(mi_area_original < 0) = 0.00001;
mi_maja_original(mi_maja_original < 0) = 0.00001;
mi_mina_original(mi_mina_original < 0) = 0.00001;
    

%% PEAK DETECTION;

%Prepare the list of the stimuli to analyses; 
cd XXX\Lips_envelope_signals\;
peak_figure = 'XXX\figures_MI_peaks\'; 
Folder = cd; Folder = fullfile(Folder);

%Prepare the list of the stimuli to analyses; 
list_v = dir('Info_lips_*.mat');
clear list_v_num

for a = 1:length(list_v)  
    list_v_num(a,1) = str2num(list_v(a).name(11:end-4));
end

[newnum, inx] = sort(list_v_num);

    
%------------------------------------%
%       Detect the peaks of MI;
%------------------------------------%
    
for p = 1:size(mi_mina_original,1)
    
   
    %detect the 3 greater MI peaks in the 1-20Hz spectrum;
    [~,locs_area] = findpeaks(mi_area_original(p,1:end),'npeaks',3,'SortStr','descend');
    [~,locs_mina] = findpeaks(mi_mina_original(p,1:end),'npeaks',3,'SortStr','descend');
    [~,locs_maja] = findpeaks(mi_maja_original(p,1:end),'npeaks',3,'SortStr','descend');
    

    %Plot and save the detected peaks for area, mina and maja together for each video;
    figure; 
    ar = subplot(1,3,1);findpeaks(mi_area_original(p,1:end),'npeaks',3,'SortStr','descend');    
    mi = subplot(1,3,2);findpeaks(mi_mina_original(p,1:end),'npeaks',3,'SortStr','descend');
    mj = subplot(1,3,3);findpeaks(mi_maja_original(p,1:end),'npeaks',3,'SortStr','descend');
    title(ar, 'AREA'); title(mi, 'MINA'); title(mj, 'MAJA'); 
    xlabel(ar, 'Freq.'); xlabel(mi, 'Freq.'); xlabel(mj, 'Freq.');
    ylabel(ar, 'M.I'); ylabel(mi, 'M.I'); ylabel(mj, 'M.I');
    xticks(ar,[1 5 10 15 20]); xticks(mi,[1 5 10 15 20]); xticks(mj,[1 5 10 15 20]);
    set(gcf,'Position',[680 233 860 745]);
    a = axes; 
    a.Visible = 'off'; 
    t1 = title(['video ' num2str(newnum(p))], 'FontSize',15, 'EdgeColor', [0 0 0]);
    t1.Visible = 'on';
    t1.Position = t1.Position + [-0.57 0.03 0];
    set(gcf,'color','w');
    print([peak_figure, 'peaks_video_' num2str(newnum(p))],'-dpng');
    close all

    %--- PEAKS AREA ---%

    if locs_area(1) > 3 &&  locs_area(1) < 18        
        Peaks_MI_area(p,1) = newnum(p); 
        Peaks_MI_area(p,2:8) = locs_area(1)-3:locs_area(1)+3;
        Peaks_MI_area(p,9:15) = mi_area_original(p, Peaks_MI_area(p,2:8));
      
        elseif locs_area(1) == 3        
            Peaks_MI_area(p,1) = newnum(p); 
            Peaks_MI_area(p,2:8) = [0 locs_area(1)-2:locs_area(1)+3];
            Peaks_MI_area(p,9:15) = [0 mi_area_original(p, Peaks_MI_area(p,3:8))];

        elseif locs_area(1) == 2                 
            Peaks_MI_area(p,1) = newnum(p); 
            Peaks_MI_area(p,2:8) = [0 0 locs_area(1)-1:locs_area(1)+3];
            Peaks_MI_area(p,9:15) = [0 0 mi_area_original(p, Peaks_MI_area(p,4:8))];        

        elseif locs_area(1) == 1              
            Peaks_MI_area(p,1) = newnum(p); 
            Peaks_MI_area(p,2:8) = [0 0 0 locs_area(1):locs_area(1)+3];
            Peaks_MI_area(p,9:15) = [0 0 0 mi_area_original(p, Peaks_MI_area(p,5:8))];         

        elseif locs_area(1) == 18        
            Peaks_MI_area(p,1) = newnum(p); 
            Peaks_MI_area(p,2:8) = [locs_area(1)-3:locs_area(1)+2 0];
            Peaks_MI_area(p,9:15) = [mi_area_original(p, Peaks_MI_area(p,2:7)) 0];

        elseif locs_area(1) == 19                
            Peaks_MI_area(p,1) = newnum(p); 
            Peaks_MI_area(p,2:8) = [locs_area(1)-3:locs_area(1)+1 0 0];
            Peaks_MI_area(p,9:15) = [mi_area_original(p, Peaks_MI_area(p,2:6)) 0 0];       

        elseif locs_area(1) == 20             
            Peaks_MI_area(p,1) = newnum(p); 
            Peaks_MI_area(p,2:8) = [locs_area(1)-3:locs_area(1) 0 0 0];
            Peaks_MI_area(p,9:15) = [mi_area_original(p, Peaks_MI_area(p,2:5)) 0 0 0];      
    end
    
    
    %--- PEAKS MINA ---%
    
    if locs_mina(1) > 3 &&  locs_mina(1) < 18        
        Peaks_MI_mina(p,1) = newnum(p); 
        Peaks_MI_mina(p,2:8) = locs_mina(1)-3:locs_mina(1)+3;
        Peaks_MI_mina(p,9:15) = mi_mina_original(p, Peaks_MI_mina(p,2:8));

        elseif locs_mina(1) == 3        
            Peaks_MI_mina(p,1) = newnum(p); 
            Peaks_MI_mina(p,2:8) = [0 locs_mina(1)-2:locs_mina(1)+3];
            Peaks_MI_mina(p,9:15) = [0 mi_mina_original(p, Peaks_MI_mina(p,3:8))];

        elseif locs_mina(1) == 2                 
            Peaks_MI_mina(p,1) = newnum(p); 
            Peaks_MI_mina(p,2:8) = [0 0 locs_mina(1)-1:locs_mina(1)+3];
            Peaks_MI_mina(p,9:15) = [0 0 mi_mina_original(p, Peaks_MI_mina(p,4:8))];        

        elseif locs_mina(1) == 1              
            Peaks_MI_mina(p,1) = newnum(p); 
            Peaks_MI_mina(p,2:8) = [0 0 0 locs_mina(1):locs_mina(1)+3];
            Peaks_MI_mina(p,9:15) = [0 0 0 mi_mina_original(p, Peaks_MI_mina(p,5:8))];         

        elseif locs_mina(1) == 18        
            Peaks_MI_mina(p,1) = newnum(p); 
            Peaks_MI_mina(p,2:8) = [locs_mina(1)-3:locs_mina(1)+2 0];
            Peaks_MI_mina(p,9:15) = [mi_mina_original(p, Peaks_MI_mina(p,2:7)) 0];

        elseif locs_mina(1) == 19                
            Peaks_MI_mina(p,1) = newnum(p); 
            Peaks_MI_mina(p,2:8) = [locs_mina(1)-3:locs_mina(1)+1 0 0];
            Peaks_MI_mina(p,9:15) = [mi_mina_original(p, Peaks_MI_mina(p,2:6)) 0 0];       

        elseif locs_mina(1) == 20             
            Peaks_MI_mina(p,1) = newnum(p); 
            Peaks_MI_mina(p,2:8) = [locs_mina(1)-3:locs_mina(1) 0 0 0];
            Peaks_MI_mina(p,9:15) = [mi_mina_original(p, Peaks_MI_mina(p,2:5)) 0 0 0];     
    end
    
    
    %--- PEAKS MAJA ---%
    
    if locs_maja(1) > 3 &&  locs_maja(1) < 18        
        Peaks_MI_maja(p,1) = newnum(p); 
        Peaks_MI_maja(p,2:8) = locs_maja(1)-3:locs_maja(1)+3;
        Peaks_MI_maja(p,9:15) = mi_maja_original(p, Peaks_MI_maja(p,2:8));
      
        elseif locs_maja(1) == 3        
            Peaks_MI_maja(p,1) = newnum(p); 
            Peaks_MI_maja(p,2:8) = [0 locs_maja(1)-2:locs_maja(1)+3];
            Peaks_MI_maja(p,9:15) = [0 mi_maja_original(p, Peaks_MI_maja(p,3:8))];

        elseif locs_maja(1) == 2                 
            Peaks_MI_maja(p,1) = newnum(p); 
            Peaks_MI_maja(p,2:8) = [0 0 locs_maja(1)-1:locs_maja(1)+3];
            Peaks_MI_maja(p,9:15) = [0 0 mi_maja_original(p, Peaks_MI_maja(p,4:8))];        

        elseif locs_maja(1) == 1              
            Peaks_MI_maja(p,1) = newnum(p); 
            Peaks_MI_maja(p,2:8) = [0 0 0 locs_maja(1):locs_maja(1)+3];
            Peaks_MI_maja(p,9:15) = [0 0 0 mi_maja_original(p, Peaks_MI_maja(p,5:8))];         

        elseif locs_maja(1) == 18        
            Peaks_MI_maja(p,1) = newnum(p); 
            Peaks_MI_maja(p,2:8) = [locs_maja(1)-3:locs_maja(1)+2 0];
            Peaks_MI_maja(p,9:15) = [mi_maja_original(p, Peaks_MI_maja(p,2:7)) 0];

        elseif locs_maja(1) == 19                
            Peaks_MI_maja(p,1) = newnum(p); 
            Peaks_MI_maja(p,2:8) = [locs_maja(1)-3:locs_maja(1)+1 0 0];
            Peaks_MI_maja(p,9:15) = [mi_maja_original(p, Peaks_MI_maja(p,2:6)) 0 0];       

        elseif locs_maja(1) == 20             
            Peaks_MI_maja(p,1) = newnum(p); 
            Peaks_MI_maja(p,2:8) = [locs_maja(1)-3:locs_maja(1) 0 0 0];
            Peaks_MI_maja(p,9:15) = [mi_maja_original(p, Peaks_MI_maja(p,2:5)) 0 0 0];      
    end
    
    clear locs_area locs_mina locs_maja 
    
end

%save the peak info in 3 separate tables;
peak_MI_AREA = array2table(Peaks_MI_area(:,[1 5 9:15]),'VariableNames',{'Video','Freq_Peak_area','minus_3hz','minus_2hz','minus_1hz','MI_FqPeak_area','plus_1hz','plus_2hz','plus_3hz'});
peak_MI_MINA = array2table(Peaks_MI_mina(:,[1 5 9:15]),'VariableNames',{'Video','Freq_Peak_mina','minus_3hz','minus_2hz','minus_1hz','MI_FqPeak_mina','plus_1hz','plus_2hz','plus_3hz'});
peak_MI_MAJA = array2table(Peaks_MI_maja(:,[1 5 9:15]),'VariableNames',{'Video','Freq_Peak_maja','minus_3hz','minus_2hz','minus_1hz','MI_FqPeak_maja','plus_1hz','plus_2hz','plus_3hz'});
 
%sort the videos in ascend order according to mina info;
[~,idx] = sortrows(peak_MI_MINA,{'Freq_Peak_mina','Video'},{'ascend','ascend'});   
peak_MI_MINA = sortrows(peak_MI_MINA,{'Freq_Peak_mina','Video'},{'ascend','ascend'}); 
mi_area_original = mi_area_original(idx,:);
mi_mina_original = mi_mina_original(idx,:);
mi_maja_original = mi_maja_original(idx,:);
save(fullfile(Folder,['MI_allstimuli', num2str(length(list_a))]), 'mi_area_original','mi_maja_original','mi_mina_original','c','f');
save(fullfile(Folder,['Peaks_allstimuli', num2str(length(list_a))]), 'peak_MI_AREA','peak_MI_MINA','peak_MI_MAJA');
clear Peaks_MI_area Peaks_MI_mina Peaks_MI_maja;

%%