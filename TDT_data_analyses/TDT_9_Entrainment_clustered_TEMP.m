%%PHASE ANALYSES WITHOUT REALIGNING;
clearvars;clc;
addpath C:\toolbox\CircHist-master;
addpath C:\toolbox\Circular_Statistics_Toolbox\;
addpath D:\Birmingham_2018\BEHAVIOR_project\ANALYSES_DATA\;

%subjects ID;
subjects = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29];

% I- Compile the hit phase of all participants together;
hit_angle_GA = zeros(length(subjects),300);
miss_angle_GA = zeros(length(subjects),300);
hit_angle_1Tone = zeros(length(subjects),100);
miss_angle_1Tone = zeros(length(subjects),100);
hit_angle_2Tones_first = zeros(length(subjects),100);
miss_angle_2Tones_first = zeros(length(subjects),100);
hit_angle_2Tones_second = zeros(length(subjects),100);
miss_angle_2Tones_second = zeros(length(subjects),100);

for s = subjects
    
    cd(['D:\Birmingham_2018\BEHAVIOR_project\ANALYSES_DATA\subj' num2str(s)]);
    load(['subj' num2str(s) '-Phase_analyses']) 
    
    hit_angle_GA(s,:) = [cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1)); NaN(300 - length(cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1))),1)]';
    miss_angle_GA(s,:) = [cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==0)); NaN(300 - length(cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==0))),1)]';
    hit_angle_1Tone(s,:) = [cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.CONDITION==1)); NaN(100 - length(cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.CONDITION==1))),1)]';
    miss_angle_1Tone(s,:) = [cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==0 & PHASE_ANALYSES.CONDITION==1)); NaN(100 - length(cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==0 & PHASE_ANALYSES.CONDITION==1))),1)]';
    hit_angle_2Tones_first(s,:) = [cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==1)); NaN(100 - length(cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==1))),1)]';
    miss_angle_2Tones_first(s,:) = [cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==0 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==1)); NaN(100 - length(cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==0 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==1))),1)]';
    hit_angle_2Tones_second(s,:) = [cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==2)); NaN(100 - length(cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==2))),1)]';
    miss_angle_2Tones_second(s,:) = [cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==0 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==2)); NaN(100 - length(cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==0 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==2))),1)]';

    clear PHASE_ANALYSES alltrials participant_parameters

end


%% PLOT phase average for each participant;
close all

condition = 1;

%Calculate the mean angle (rad) for each participant;
mean_hit = zeros(length(subjects),1);
mean_miss = zeros(length(subjects),1);
         
for s = subjects
    
    if condition == 3
        mean_hit(s,1) = circ_mean(hit_angle_GA(s,~isnan(hit_angle_GA(s,:)))');
        mean_miss(s,1) = circ_mean(miss_angle_GA(s,~isnan(miss_angle_GA(s,:)))');
        title_plot = 'Average';
    elseif condition == 1
        mean_hit(s,1) = circ_mean(hit_angle_1Tone(s,~isnan(hit_angle_1Tone(s,:)))');
        mean_miss(s,1) = circ_mean(miss_angle_1Tone(s,~isnan(miss_angle_1Tone(s,:)))');
        title_plot = '1Tone';
    elseif condition == 21
        mean_hit(s,1) = circ_mean(hit_angle_2Tones_first(s,~isnan(hit_angle_2Tones_first(s,:)))');
        mean_hit(s,2) = s;
        mean_miss(s,1) = circ_mean(miss_angle_2Tones_first(s,~isnan(miss_angle_2Tones_first(s,:)))');
        mean_miss(s,2) = s;        
        title_plot = '2Tones-T1';
    elseif condition == 22
        mean_hit(s,1) = circ_mean(hit_angle_2Tones_second(s,~isnan(hit_angle_2Tones_second(s,:)))');
        mean_hit(s,2) = s;
        mean_miss(s,1) = circ_mean(miss_angle_2Tones_second(s,~isnan(miss_angle_2Tones_second(s,:)))');
        mean_miss(s,2) = s; 
        title_plot = '2Tones-T2';
    end
    
end


%------- PLOT THE MEAN VECTOR OF HIT/MISS IN EACH CONDITION --------%

%Remove the bad subjects to exclude from statistic analyses;
bad_subjects = [1 4 25 26 27];
mean_hit(bad_subjects,:) = [];
mean_miss(bad_subjects,:) = [];

%Plot the mean angle of each participant;
close all;
figure;
%hits;
x1 = subplot(1, 2, 1);
circ_plot(mean_hit(:,1),'pretty', [],[],false,true,'linewidth',2,'color','r');
pause(0.001);
ax1 = gca;
ax1.XRuler.Axle.LineStyle = 'none'; 
ax1.YRuler.Axle.LineStyle = 'none'; 
%miss;
subplot(1, 2, 2);
x2 = circ_plot(mean_miss(:,1),'pretty', [],[],false,true,'linewidth',2,'color','r');
pause(0.001);
ax2 = gca;
ax2.XRuler.Axle.LineStyle = 'none'; 
ax2.YRuler.Axle.LineStyle = 'none'; 
title(x1,'HITS','Units', 'normalized', 'Position', [0.5, 1.2, 0], 'FontSize',15);
title(x2, 'MISS','Units', 'normalized', 'Position', [0.5, 1.2, 0], 'FontSize',15);
%final format and put mean title( ridiculously complicated;
axes('Units','Normal');
m_title = title(title_plot);
set(gca,'visible','off')
set(m_title,'visible','on','FontSize',20)
set(gcf,'color','w');

%% PLOT THE HITS in T1 and T2 conditions (without separating the clusters);
    
%Calculate the mean angle (rad) for each participant;
mean_hit_T1 = zeros(length(subjects),1);
mean_hit_T2 = zeros(length(subjects),1);

for s = subjects 

mean_hit_T1(s,1) = circ_mean(hit_angle_2Tones_first(s,~isnan(hit_angle_2Tones_first(s,:)))');
mean_hit_T1(s,2) = s;      
mean_hit_T2(s,1) = circ_mean(hit_angle_2Tones_second(s,~isnan(hit_angle_2Tones_second(s,:)))');
mean_hit_T2(s,2) = s;

end

[r_hit_T1, z_hit_T1] = circ_rtest(mean_hit_T1(:,1));
[r_hit_T2, z_hit_T2] = circ_rtest(mean_hit_T2(:,1));



%------- PLOT THE MEAN VECTOR OF HITS IN T1 and T2 CONDITIONS --------%

%Remove the bad subjects to exclude from statistic analyses;
bad_subjects = [1 4 25 26 27];
mean_hit_T1(bad_subjects,:) = [];
mean_hit_T2(bad_subjects,:) = [];

%Plot the mean angle of each participant;
close all;
figure;
%hits T1;
subplot(1, 2, 1);
x1 = circ_plot2(mean_hit_T1(:,1),'pretty', [],[],false,true,'linewidth',3,'color','r');
pause(0.001);
ax1 = gca;
ax1.XRuler.Axle.LineStyle = 'none'; 
ax1.YRuler.Axle.LineStyle = 'none'; 
%hits T2;
subplot(1, 2, 2);
x2 = circ_plot2(mean_hit_T2(:,1),'pretty', [],[],false,true,'linewidth',3,'color','r');
pause(0.001);
ax2 = gca;
ax2.XRuler.Axle.LineStyle = 'none'; 
ax2.YRuler.Axle.LineStyle = 'none'; 
title(x1,'Tone T1','Units', 'normalized', 'Position', [0.5, 1.2, 0], 'FontSize',15);
title(x2, 'Tone T2','Units', 'normalized', 'Position', [0.5, 1.2, 0], 'FontSize',15);
set(gcf,'color','w','Position',[674 546 560 420]);
%     axes('Units','Normal');
%     m_title = title('HITS T1 vs. T2');
%     set(gca,'visible','off')
%     set(m_title,'visible','on','FontSize',20)    
    
%% find polar coordinate of each participant; 
  
%Calculate the polar coordinates of each participant;
for s = 1:length(mean_hit_T2(:,1))

    phi = mean_hit_T2(s,1);
    exp(1i*phi);  
    zm(s,1) = mean_hit_T2(s,2);
    zm(s,2) = real(exp(1i*phi)); % X;
    zm(s,3) = imag(exp(1i*phi)); % Y;

end

%Separate the 2 clusters os participants based on real/imag coordinates; 
X1 = 0.59;
Y1 = 0.81;
Y2 = -0.81;   
X3 = -0.98;
Y3 = -0.25;

%Discriminate the 2 cluster populations;
cluster_1 = zeros(length(zm(:,1)),1);
cluster_2 = zeros(length(zm(:,1)),1);

for i= 1:length(zm(:,1))

    if zm(i,2) > X1 && zm(i,3) < Y1 && zm(i,3) > Y2
        cluster_1(i,1) = zm(i,1);
    elseif zm(i,2) > X3 && zm(i,2) < X1 &&  zm(i,3) < Y3
        cluster_2(i,1) = zm(i,1);
    end

end

cluster_1 = cluster_1(cluster_1~=0);
cluster_2 = cluster_2(cluster_2~=0);
idx_cluster_1 = ~ismember(mean_hit_T2(:,2),cluster_1);
idx_cluster_2 = ~ismember(mean_hit_T2(:,2),cluster_2);

%Store the mean phase T1/T2 separately for the 2 groups;
mean_hit_T1_cluster1 = mean_hit_T1;
mean_hit_T1_cluster2 = mean_hit_T1;
mean_hit_T2_cluster1 = mean_hit_T2;
mean_hit_T2_cluster2 = mean_hit_T2;
mean_hit_T1_cluster1(idx_cluster_1,:) = [];
mean_hit_T1_cluster2(idx_cluster_2,:) = [];
mean_hit_T2_cluster1(idx_cluster_1,:) = [];
mean_hit_T2_cluster2(idx_cluster_2,:) = [];


%------- PLOT THE MEAN VECTOR OF HITS IN T1 and T2 CONDITIONS FOR THE 2 CLUSTER POPULATIONS  --------%

%Plot the mean angle of each participant;
close all;
figure;
%hits;
subplot(1, 4, 1);
x1 = circ_plot3(mean_hit_T1_cluster1(:,1),'pretty', [],[],false,true,'linewidth',2,'color','r');
pause(0.001);
ax1 = gca;
ax1.XRuler.Axle.LineStyle = 'none'; 
ax1.YRuler.Axle.LineStyle = 'none'; 
%miss;
subplot(1, 4, 2);
x2 = circ_plot3(mean_hit_T1_cluster2(:,1),'pretty', [],[],false,true,'linewidth',2,'color','r');
pause(0.001);
ax2 = gca;
ax2.XRuler.Axle.LineStyle = 'none'; 
ax2.YRuler.Axle.LineStyle = 'none'; 
pause(1);
subplot(1, 4, 3);
x3 = circ_plot3(mean_hit_T2_cluster1(:,1),'pretty', [],[],false,true,'linewidth',2,'color','r');
pause(0.001);    
ax3 = gca;
ax3.XRuler.Axle.LineStyle = 'none'; 
ax3.YRuler.Axle.LineStyle = 'none'; 
subplot(1, 4, 4);
x4 = circ_plot3(mean_hit_T2_cluster2(:,1),'pretty', [],[],false,true,'linewidth',2,'color','r');
pause(0.001);    
ax4 = gca;
ax4.XRuler.Axle.LineStyle = 'none'; 
ax4.YRuler.Axle.LineStyle = 'none';     
%     title(x1,'T1 Clust.1','Units', 'normalized', 'Position', [0.5, 1.2, 0], 'FontSize',15);
%     title(x2, 'T1 Clust.2','Units', 'normalized', 'Position', [0.5, 1.2, 0], 'FontSize',15);
%     title(x3,'T2 Clust.1','Units', 'normalized', 'Position', [0.5, 1.2, 0], 'FontSize',15);
%     title(x4, 'T2 Clust.2','Units', 'normalized', 'Position', [0.5, 1.2, 0], 'FontSize',15);
set(gcf,'color','w');


%     axes('Units','Normal');
%     m_title = title('HITS T1 vs. T2 Clustered');
%     set(gca,'visible','off')
%     set(m_title,'visible','on','FontSize',20)
    
%%