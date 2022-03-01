clearvars;clc;
addpath C:\toolbox\CircHist-master; addpath C:\toolbox\Circular_Statistics_Toolbox\; cd C:\ANALYSES_DATA\;

%subjects ID;
subjects = 1:29;

condition = 21;

for s = subjects
    
    cd(['C:\ANALYSES_DATA\subj' num2str(s)]);
    load(['subj' num2str(s) '-Phase_analyses']);

    switch condition
    
        case 11    
            
            %T1 window of TwoTones Condition; 
            ori_hit(s,:) = [cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.CONDITION==1 & PHASE_ANALYSES.TONE_POSITION==1 & PHASE_ANALYSES.ONSET/44100 <= 2.5)); NaN(100 - length(cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.CONDITION==1 & PHASE_ANALYSES.TONE_POSITION==1 & PHASE_ANALYSES.ONSET/44100 <= 2.5))),1)]';
            ori_miss(s,:) = [cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==0 & PHASE_ANALYSES.CONDITION==1 & PHASE_ANALYSES.TONE_POSITION==1 & PHASE_ANALYSES.ONSET/44100 <= 2.5)); NaN(100 - length(cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==0 & PHASE_ANALYSES.CONDITION==1 & PHASE_ANALYSES.TONE_POSITION==1 & PHASE_ANALYSES.ONSET/44100 <= 2.5))),1)]';
            ori_mean_hit(s,1) = circ_mean(ori_hit(s,~isnan(ori_hit(s,:)))');
            ori_mean_miss(s,1) = circ_mean(ori_miss(s,~isnan(ori_miss(s,:)))');
            title_plot = 'Single Tone - First';
            
        case 12
            
            %T2 window of TwoTones Condition; 
            ori_hit(s,:) = [cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.CONDITION==1 & PHASE_ANALYSES.TONE_POSITION==1 & PHASE_ANALYSES.ONSET/44100 > 2.5)); NaN(100 - length(cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.CONDITION==1 & PHASE_ANALYSES.TONE_POSITION==1 & PHASE_ANALYSES.ONSET/44100 > 2.5))),1)]';
            ori_miss(s,:) = [cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==0 & PHASE_ANALYSES.CONDITION==1 & PHASE_ANALYSES.TONE_POSITION==1 & PHASE_ANALYSES.ONSET/44100 > 2.5)); NaN(100 - length(cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==0 & PHASE_ANALYSES.CONDITION==1 & PHASE_ANALYSES.TONE_POSITION==1 & PHASE_ANALYSES.ONSET/44100 > 2.5))),1)]';
            ori_mean_hit(s,1) = circ_mean(ori_hit(s,~isnan(ori_hit(s,:)))');
            ori_mean_miss(s,1) = circ_mean(ori_miss(s,~isnan(ori_miss(s,:)))');
            title_plot = 'Single Tone - Second';
            
        case 21    
            
            %T1 window of TwoTones Condition; 
            ori_hit(s,:) = [cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==1)); NaN(100 - length(cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==1))),1)]';
            ori_miss(s,:) = [cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==0 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==1)); NaN(100 - length(cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==0 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==1))),1)]';
            ori_mean_hit(s,1) = circ_mean(ori_hit(s,~isnan(ori_hit(s,:)))');
            ori_mean_miss(s,1) = circ_mean(ori_miss(s,~isnan(ori_miss(s,:)))');
            title_plot = 'Two Tones - First';
            
        case 22
            
            %T2 window of TwoTones Condition; 
            ori_hit(s,:) = [cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==2)); NaN(100 - length(cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==2))),1)]';
            ori_miss(s,:) = [cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==0 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==2)); NaN(100 - length(cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==0 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==2))),1)]';
            ori_mean_hit(s,1) = circ_mean(ori_hit(s,~isnan(ori_hit(s,:)))');
            ori_mean_miss(s,1) = circ_mean(ori_miss(s,~isnan(ori_miss(s,:)))');
            title_plot = 'Two Tones - Second';
    end
    
    clear PHASE_ANALYSES alltrials participant_parameters

end
   
%Remove the bad subjects to exclude from statistic analyses;
bad_subjects = [1 4 25 26 27];
ori_mean_hit(bad_subjects,:) = [];
ori_mean_miss(bad_subjects,:) = [];   

%Plot the mean angle of each participant;
close all;
figure;
%hits (green line);
x1 = circ_plot(ori_mean_hit(:,1),'pretty', [],[],false,true,'linewidth',3,'color','g');
pause(0.001);
ax1 = gca;
ax1.XRuler.Axle.LineStyle = 'none'; 
ax1.YRuler.Axle.LineStyle = 'none'; 
hold on
%misses (red line);
x2 = circ_plot(ori_mean_miss(:,1),'pretty', [],[],false,true,'linewidth',3,'color','r');
pause(0.001);
ax2 = gca;
ax2.XRuler.Axle.LineStyle = 'none'; 
ax2.YRuler.Axle.LineStyle = 'none'; 
title(title_plot,'Position', [-1.2, 1, 0], 'FontSize',15);
set(gcf,'color','w','Position',[ 425.0000  256.2000  560.0000  388.0000]);

%Rayleigh's Test on the T1 and T2 mean phase distributions; 
[p,~] = circ_rtest(ori_mean_hit);
r = circ_r(ori_mean_hit);
disp(['condition: ' title_plot]);
disp(['p-value: ' num2str(p)]);
disp(['r_length: ' num2str(r)]);

%%