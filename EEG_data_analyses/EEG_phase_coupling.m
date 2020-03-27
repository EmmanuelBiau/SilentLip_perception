%Now let's do Phase coupling analysis between the left audio and left visual sources reconstructed (and flipped) signals;
clearvars;clc; 
addpath XXX\Circular_Statistics_Toolbox\;
cd XXX;

%subjects ID;
subjects = [2:9 11:17 19:26];

scount = 0;
for s = subjects
    
    scount = scount + 1;
    cd (['XXX\subj',num2str(s)]);
    load(['visual_entrainment_source_realigned_subj' num2str(s)]);
    
    source_audio_left = source_left_audio_realigned; 
    source_visual_left = source_left_visual_realigned;
    
    %GRam-Schmidt orthogonality correction;
    for i = 1:length(source_audio_left.trial) 

        X = source_visual_left.trial{1,i}; 
        Y = source_audio_left.trial{1,i};

        %project left_audio signal orthogonally onto the line spanned by visual signal.
        Yp = (sum(X.*Y )./sum(X.*X ))*X; 
        Yo = Y - Yp;
        source_audio_left_temp.trial{1,i} = Yo;

        %This value should be zero (or very close to zero);
        check_left(i,1)= dot(Yo,X); 

        clear X Y Z Yp Zp Yo Zo

    end

    %Now replace the orthogonal signal in the original source data;
    source_audio_left.trial = source_audio_left_temp.trial;

    %Get the phase angle in each freq.peak trials with adaptative filter band according to the frequency of MI peak in the videos; 
    cfg = [];
    cfg.bpfilter = 'yes';
    cfg.bpfreq = [2 6];        
    cfg.bpfilttype = 'fir';
    cfg.hilbert = 'angle';
    cfg.trials = find(cell2mat(cellfun(@(x) x.freq_peak == 4, source_audio_left.trialinfo, 'UniformOutput', false)));
    source_audio_lefth_4 = ft_preprocessing(cfg, source_audio_left);
    source_visual_lefth_4 = ft_preprocessing(cfg, source_visual_left);
    
    cfg = [];
    cfg.bpfilter = 'yes';
    cfg.bpfreq = [3 7];        
    cfg.bpfilttype = 'fir';
    cfg.hilbert = 'angle';
    cfg.trials = find(cell2mat(cellfun(@(x) x.freq_peak == 5, source_audio_left.trialinfo, 'UniformOutput', false)));
    source_audio_lefth_5 = ft_preprocessing(cfg, source_audio_left);
    source_visual_lefth_5 = ft_preprocessing(cfg, source_visual_left);
    
    cfg = [];
    cfg.bpfilter = 'yes';
    cfg.bpfreq = [4 8];        
    cfg.bpfilttype = 'fir';
    cfg.hilbert = 'angle';
    cfg.trials = find(cell2mat(cellfun(@(x) x.freq_peak == 6, source_audio_left.trialinfo, 'UniformOutput', false)));
    source_audio_lefth_6 = ft_preprocessing(cfg, source_audio_left);
    source_visual_lefth_6 = ft_preprocessing(cfg, source_visual_left);
    
    cfg = [];
    cfg.bpfilter = 'yes';
    cfg.bpfreq = [5 9];        
    cfg.bpfilttype = 'fir';
    cfg.hilbert = 'angle';
    cfg.trials = find(cell2mat(cellfun(@(x) x.freq_peak == 7, source_audio_left.trialinfo, 'UniformOutput', false)));
    source_audio_lefth_7 = ft_preprocessing(cfg, source_audio_left);
    source_visual_lefth_7 = ft_preprocessing(cfg, source_visual_left);
    
    cfg = [];
    cfg.bpfilter = 'yes';
    cfg.bpfreq = [6 10];        
    cfg.bpfilttype = 'fir';
    cfg.hilbert = 'angle';
    cfg.trials = find(cell2mat(cellfun(@(x) x.freq_peak == 8, source_audio_left.trialinfo, 'UniformOutput', false)));
    source_audio_lefth_8 = ft_preprocessing(cfg, source_audio_left);
    source_visual_lefth_8 = ft_preprocessing(cfg, source_visual_left);
            
    %Structure again to a fieldtrip structure;
    source_audio_lefth.fsample = 512;
    source_audio_lefth.label = source_audio_lefth_4.label; 
    source_audio_lefth.time = [source_audio_lefth_4.time';source_audio_lefth_5.time';source_audio_lefth_6.time';source_audio_lefth_7.time';source_audio_lefth_8.time']';
    source_audio_lefth.sampleinfo = [source_audio_lefth_4.sampleinfo;source_audio_lefth_5.sampleinfo;source_audio_lefth_6.sampleinfo;source_audio_lefth_7.sampleinfo;source_audio_lefth_8.sampleinfo];
    source_audio_lefth.trialinfo = [source_audio_lefth_4.trialinfo;source_audio_lefth_5.trialinfo;source_audio_lefth_6.trialinfo;source_audio_lefth_7.trialinfo;source_audio_lefth_8.trialinfo];
    source_audio_lefth.trial = [source_audio_lefth_4.trial';source_audio_lefth_5.trial';source_audio_lefth_6.trial';source_audio_lefth_7.trial';source_audio_lefth_8.trial']';
    
    source_visual_lefth.fsample = 512;
    source_visual_lefth.label = source_visual_lefth_4.label; 
    source_visual_lefth.time = [source_visual_lefth_4.time';source_visual_lefth_5.time';source_visual_lefth_6.time';source_visual_lefth_7.time';source_visual_lefth_8.time']';
    source_visual_lefth.sampleinfo = [source_visual_lefth_4.sampleinfo;source_visual_lefth_5.sampleinfo;source_visual_lefth_6.sampleinfo;source_visual_lefth_7.sampleinfo;source_visual_lefth_8.sampleinfo];
    source_visual_lefth.trialinfo = [source_visual_lefth_4.trialinfo;source_visual_lefth_5.trialinfo;source_visual_lefth_6.trialinfo;source_visual_lefth_7.trialinfo;source_visual_lefth_8.trialinfo];
    source_visual_lefth.trial = [source_visual_lefth_4.trial';source_visual_lefth_5.trial';source_visual_lefth_6.trial';source_visual_lefth_7.trial';source_visual_lefth_8.trial']';
    
    clear source_audio_lefth_4 source_audio_lefth_5 source_audio_lefth_6 source_audio_lefth_7 source_audio_lefth_8
    
    
    %--------------------------------------------------------------------------------
   
    %Now calculate the phase angle difference between left visual and audio sources for each time-point of the chose time-window;

    cfg = [];
    source_audio_lefth0 = ft_selectdata(cfg, source_audio_lefth);
    source_visual_lefth0 = ft_selectdata(cfg, source_visual_lefth);   
    lah = source_audio_lefth0.trial;
    lvh = source_visual_lefth0.trial;

    %Chose the time-windows of interest:
    
    t1_start = 0.5;
    t1_end = 2;
    
    t2_start = 3;
    t2_end = 4.5;
    
    fs = source_audio_lefth.fsample; %sample of data (512Hz);  
    tw1 = fs*(2+t1_start)+1:fs*(2+t1_end)+1; 
    tw2 = fs*(2+t2_start)+1:fs*(2+t2_end)+1;  

    
    for t = 1:numel(lah)
        
        %calculate the phase in audio and visual ROI signals;
        ulah0{1,t} = unwrap(lah{1,t});
        ulvh0{1,t} = unwrap(lvh{1,t});
        
        %calculate the difference between the two phases (phaseA- phaseV);
        hpda_T1{1,t} = ulah0{1,t}(tw1) - ulvh0{1,t}(tw1);
        hpda_T2{1,t} = ulah0{1,t}(tw2) - ulvh0{1,t}(tw2);
             
        %Compare the mean diff.Phase to a particular angle;
        %chose the reference angle (Here 0 because we want to compare the distance of the distribution to 0 angle);
        ang1 = 0;  
        ang2 = 0;  
        pT1 = ones(1,length(hpda_T1{1,t}))*(ang1)*pi/180;  
        pT2 = ones(1,length(hpda_T2{1,t}))*(ang2)*pi/180;  
        
        hda_T1(t,:) = circ_r([hpda_T1{1,t};pT1],[],[],1);    
        hda_T2(t,:) = circ_r([hpda_T2{1,t};pT2],[],[],1);
        
    end
    
    %Resultant mean phase difference;
    mhd_T1_T2(:,1) = mean(hda_T1,2);
    mhd_T1_T2(:,2) = mean(hda_T2,2);
    gavar_T1_T2(scount,[1 2]) = mean(mhd_T1_T2,1);

    %Bin phase angle diff for each trial;
    for i = 1:length(hpda_T1)
    phase_count_T1(i,:) = hpda_T1{i};
    end
    
    for i = 1:length(hpda_T2)
    phase_count_T2(i,:) = hpda_T2{i};
    end
    
    %average the diff of phase at each time-point of T1 and T2 across trials for each participants;
    count_T1(scount,:) = circ_mean(phase_count_T1); 
    count_T2(scount,:) = circ_mean(phase_count_T2); 
    
    %store orthogonality for all subjects;
    check_ortho{1,scount}  = check_left;
    keep scount fs gavar_T1_T2 phase_count_T1 phase_count_T2 count_T1 count_T2 check_ortho
  
end

%Average across participants and get a mean phase difference for each time-point in T1 and T2 time-windows separately;
GA_T1_T2 = [];
GA_T1_T2(1,:) = circ_mean(count_T1);
GA_T1_T2(2,:) = circ_mean(count_T2);

%Get the mean angle difference in the two Time-windows (degree);
patience = msgbox('wait a second');

alpha1 = circ_rad2ang(GA_T1_T2(1,:)); 
alpha2 = circ_rad2ang(GA_T1_T2(2,:));    
nBins = 500;
mean_alpha1 = CircHist(alpha1, nBins);
mean_angleT1 = mean_alpha1.avgAng;
mean_alpha2 = CircHist(alpha2, nBins);
mean_angleT2 = mean_alpha2.avgAng;
delete(patience);
close all; clear mean_alpha1 mean_alpha2 fH subAx1 subAx2;

%% Plot Phase coupling in T1 and T2 time windows between left audio and left visual in Silent movie condition;
addpath XXX\CircHist-master; %you need this toolbox for ploting the phase distribution;

%convert phase radians units in degrees;
alpha1 = circ_rad2ang(GA_T1_T2(1,:)); 
alpha2 = circ_rad2ang(GA_T1_T2(2,:)); 

%Choose number of bins for plotting;
nBins = 500;

%T1 Tone;
close all;
fH = figure('Visible', 'on');
subAx1 = subplot(1, 2, 1, polaraxes);
condition_0 = CircHist(alpha1, nBins, 'parent', subAx1);
thetaticks(condition_0.polarAxs, [0, 360]);
condition_0.polarAxs.ThetaAxis.MinorTickValues = [];
%Make rho-axes equal for both diagrams
maxRho = max([max(rlim(subAx1))]);
newLimits = [min(rlim(subAx1)), maxRho];
condition_0.setRLim(newLimits);
condition_0.polarAxs.ThetaZeroLocation = 'right';
%sync condition;
condition_0.polarAxs.ThetaTick = [0 90 180 270 360]; % change major ticks
condition_0.polarAxs.ThetaTickLabel = {'\bf0','\bf\pi/2','\bf+/-\pi','\bf-\pi/2'};
condition_0.polarAxs.FontSize = 14;
condition_0.avgAngH.LineStyle = 'none';
%draw resultant vector r as arrow
delete(condition_0.rH)
r1 = rlim(condition_0.polarAxs); % get current limits
mean_vect1 = condition_0.drawArrow(condition_0.avgAng, condition_0.r * range(r1), 'HeadWidth', 10, 'LineWidth', 2, 'Color', 'r');
%T2 Tone;
subAx2 = subplot(1, 2, 2, polaraxes);
condition_180 = CircHist(alpha2, nBins, 'parent', subAx2);
thetaticks(condition_180.polarAxs, [0, 360]);
condition_180.setRLim(newLimits);
condition_180.polarAxs.ThetaZeroLocation = 'right';
%async condition;
condition_180.polarAxs.ThetaAxis.MinorTickValues = [];
thetaticks(condition_180.polarAxs, [0, 360]);
condition_180.polarAxs.ThetaTick = [0 90 180 270 360]; % change major ticks
condition_180.polarAxs.ThetaTickLabel = {'\bf0','\bf\pi/2','\bf+/-\pi','\bf-\pi/2'};
condition_180.polarAxs.FontSize = 14;
condition_180.avgAngH.LineStyle = 'none';
%draw resultant vector r as arrow
delete(condition_180.rH)
r2 = rlim(condition_180.polarAxs);
mean_vect2 = condition_180.drawArrow(condition_180.avgAng, condition_180.r * range(r2), 'HeadWidth', 10, 'LineWidth', 2, 'Color', 'r');
subAx1.Title.String = 'Tone T1';
subAx2.Title.String = 'Tone T2';
set(gcf,'color','w','Position', [2200 536 892 390]);

%Test the difference of mean angle between the 2 distributions;
[pval_kuipt_T1vsT2, k_T1vsT2, K_T1vsT2] = circ_kuipertest(GA_T1_T2(1,:), GA_T1_T2(2,:));

%% Plot Vector length entrainment to zeros in T1 and T2 conditions;

%Inputs for raincloud plots;
dat = gavar_T1_T2(:,1:2);          
col = [0 0 0];                        
dotsize = 15;                              
smooth = 200;    

%Generate your plot by calling the function;
close all;
vert_rain_plot(dat,col,dotsize,smooth)
pause(0.01);

%configure axis;
ax = gca;
ax.Box = 'off';
ax.LineWidth = 2;
ax.TickDir = 'out';
ax.XTick = [1 2];
ax.XTickLabel = {'Tone T1','Tone T2'};
ax.XAxis.TickLength = [0 0];
ax.XRuler.Axle.LineStyle = 'none'; 
ax.XAxis.FontSize = 14;
ax.XAxis.FontWeight = 'bold';
ax.YAxis.Limits = [0.55 0.75];
ax.YTick = [0.55 0.65 0.75];
ax.YAxis.FontSize = 12;
ax.YAxis.FontWeight = 'bold';
ylabel('Mean vector length (r)','FontSize',14, 'Layer','front', 'FontWeight', 'bold');
hold on 
x2 = linspace(1,1.6);
y2 = linspace(0.74,0.74);
plot(x2,y2, 'Color', [0 0 0],'LineWidth',0.5);
text(1.25,0.745,'*','FontSize', 20,'FontWeight','light');
set(gcf,'color','w','Positio',[674 546 329 420]);

%% Convert the mean distance phase-coupling to zeros in time;
cd XXX\;
load video_info; %This contains the same info of the simuli used for the MI calculations (i.e. stim_freqpeak);

%mean frequency across stimuli;
mean_Fq = mean(Stim_infos.FREQPEAK);

%Convert mean angle in T1 and T2 in corresponding ms for the mean_Fq; 
pd_mean_fq = (1000/mean_Fq)/(2*pi);

%convert difference angle in radian;
diff_angT1 = (360 - mean_angleT1)*pi/180;
diff_angT2 = (360 - mean_angleT2)*pi/180;
lagVA_T1 = pd_mean_fq*diff_angT1;
lagVA_T2 = pd_mean_fq*diff_angT2;

%% 