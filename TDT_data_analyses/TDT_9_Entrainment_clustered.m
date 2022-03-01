%%Here, we will segregate the two populations of entrainers based on entrainment at T2 tone in the Two Tone condition (where entrainment took place);
clearvars;clc;
addpath C:\toolbox\CircHist-master; addpath C:\toolbox\Circular_Statistics_Toolbox\; 
cd D:\Birmingham_2018\BEHAVIOR_project\ANALYSES_DATA\;

%subjects ID;
subjects = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29];

% I- Compile the hit phases of all participants together;
clear hit_angle_2Tones_first hit_angle_2Tones_second

for s = subjects
    
    cd(['D:\Birmingham_2018\BEHAVIOR_project\ANALYSES_DATA\subj' num2str(s)]);
    load(['subj' num2str(s) '-Phase_analyses']) 

    hit_angle_2Tones_first(s,:) = [cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==1)); NaN(100 - length(cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==1))),1)]';
    hit_angle_2Tones_second(s,:) = [cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==2)); NaN(100 - length(cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==2))),1)]';

    clear PHASE_ANALYSES alltrials participant_parameters

end

%------- PLOT THE HITS in T1 and T2 Windows (without separating the 2 populations of entrainers) --------%
    
%Calculate the mean angle (rad) for each participant;
clear mean_hit_T1 mean_hit_T2

for s = subjects 

mean_hit_T1(s,1) = circ_mean(hit_angle_2Tones_first(s,~isnan(hit_angle_2Tones_first(s,:)))');
mean_hit_T1(s,2) = s;      
mean_hit_T2(s,1) = circ_mean(hit_angle_2Tones_second(s,~isnan(hit_angle_2Tones_second(s,:)))');
mean_hit_T2(s,2) = s;

end

%Remove the bad subjects to exclude from statistic analyses;
bad_subjects = [1 4 25 26 27];
mean_hit_T1(bad_subjects,:) = [];
mean_hit_T2(bad_subjects,:) = [];

%Rayleigh's tests on mean phase distributions;
[r_hit_T1, z_hit_T1] = circ_rtest(mean_hit_T1(:,1));
[r_hit_T2, z_hit_T2] = circ_rtest(mean_hit_T2(:,1));

%Plot the mean angle of each participant;
close all;
figure;
%hits T1;
subplot(1, 2, 1);
x1 = circ_plot(mean_hit_T1(:,1),'pretty', [],[],false,true,'linewidth',3,'color','r');
pause(0.001);
ax1 = gca;
ax1.XRuler.Axle.LineStyle = 'none'; 
ax1.YRuler.Axle.LineStyle = 'none'; 
%hits T2;
subplot(1, 2, 2);
x2 = circ_plot(mean_hit_T2(:,1),'pretty', [],[],false,true,'linewidth',3,'color','r');
pause(0.001);
ax2 = gca;
ax2.XRuler.Axle.LineStyle = 'none'; 
ax2.YRuler.Axle.LineStyle = 'none'; 
title(x1,'Tone T1','Units', 'normalized', 'Position', [0.5, 1.2, 0], 'FontSize',15);
title(x2, 'Tone T2','Units', 'normalized', 'Position', [0.5, 1.2, 0], 'FontSize',15);
set(gcf,'color','w','Position',[674 546 560 420]);

%% Separate the two populations of entrainers based on their main phase; 
  
%Calculate the polar coordinates of each participant;
clear phi zm

for s = 1:length(mean_hit_T2(:,1))

    phi = mean_hit_T2(s,1);
    exp(1i*phi);  
    zm(s,1) = mean_hit_T2(s,2);
    zm(s,2) = real(exp(1i*phi)); % X;
    zm(s,3) = imag(exp(1i*phi)); % Y;

end

%Separate the 2 clusters of participants based on real/imag coordinates; To adapt to the range of phases we want to separate;
X1 = 0.59;
Y1 = 0.81;
Y2 = -0.81;   
X3 = -0.98;
Y3 = -0.25;

%Segragate the 2 clusters of populations;
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

clear idx_cluster_1 idx_cluster_2


%------- PLOT THE MEAN VECTOR OF HITS IN T1 and T2 CONDITIONS FOR THE 2 ENTRAINER POPULATIONS  --------%

%Plot the mean angle of each participant;
close all;
figure;
subplot(1, 4, 1);
x1 = circ_plot(mean_hit_T1_cluster1(:,1),'pretty', [],[],false,true,'linewidth',2,'color','r');
pause(0.001);
ax1 = gca;
ax1.XRuler.Axle.LineStyle = 'none'; 
ax1.YRuler.Axle.LineStyle = 'none'; 
subplot(1, 4, 2);
x2 = circ_plot(mean_hit_T1_cluster2(:,1),'pretty', [],[],false,true,'linewidth',2,'color','r');
pause(0.001);
ax2 = gca;
ax2.XRuler.Axle.LineStyle = 'none'; 
ax2.YRuler.Axle.LineStyle = 'none'; 
pause(1);
subplot(1, 4, 3);
x3 = circ_plot(mean_hit_T2_cluster1(:,1),'pretty', [],[],false,true,'linewidth',2,'color','r');
pause(0.001);    
ax3 = gca;
ax3.XRuler.Axle.LineStyle = 'none'; 
ax3.YRuler.Axle.LineStyle = 'none'; 
subplot(1, 4, 4);
x4 = circ_plot(mean_hit_T2_cluster2(:,1),'pretty', [],[],false,true,'linewidth',2,'color','r');
pause(0.001);    
ax4 = gca;
ax4.XRuler.Axle.LineStyle = 'none'; 
ax4.YRuler.Axle.LineStyle = 'none';     
title(x1,'T1 Entr.1','Units', 'normalized', 'Position', [0.5, 1.2, 0], 'FontSize',15);
title(x2, 'T1 Entr.2','Units', 'normalized', 'Position', [0.5, 1.2, 0], 'FontSize',15);
title(x3,'T2 Entr.1','Units', 'normalized', 'Position', [0.5, 1.2, 0], 'FontSize',15);
title(x4, 'T2 Entr.2','Units', 'normalized', 'Position', [0.5, 1.2, 0], 'FontSize',15);
set(gcf,'color','w','Position',[676 434 865 420]);


%% Plot Dprimes/RTs Clustered data IN TWO TONES CONDITION;
clear;clc;
cd D:\Birmingham_2018\BEHAVIOR_project\ANALYSES_DATA\;
load allsubj_T1vsT2_performances;

cluster1 = [2 3 5 6 9 10 11 13 15 19 22]; 
cluster2 = [7 8 12 14 17 20 21 23 24 29];

%First separate the data from cluster1 and cluster2 (because different
%lengths);
MEAN_HIT_cluster1 = TwoTones_T1vsT2_performances(cluster1,:);
MEAN_HIT_cluster2 = TwoTones_T1vsT2_performances(cluster2,:);
clear dprimes_clustered rt_clustered
dprimes_clustered(:,1)= MEAN_HIT_cluster1.dp_T1;
dprimes_clustered(:,2)= [MEAN_HIT_cluster2.dp_T1;NaN*(length(MEAN_HIT_cluster1.dp_T1)-length(MEAN_HIT_cluster2.dp_T1))]; 
dprimes_clustered(:,3)= MEAN_HIT_cluster1.dp_T2;
dprimes_clustered(:,4)= [MEAN_HIT_cluster2.dp_T2;NaN*(length(MEAN_HIT_cluster1.dp_T2)-length(MEAN_HIT_cluster2.dp_T2))];     
rt_clustered(:,1)= MEAN_HIT_cluster1.rt_T1*1000;
rt_clustered(:,2)= [MEAN_HIT_cluster2.rt_T1;NaN*(length(MEAN_HIT_cluster1.rt_T1)-length(MEAN_HIT_cluster2.rt_T1))]*1000; 
rt_clustered(:,3)= MEAN_HIT_cluster1.rt_T2*1000;
rt_clustered(:,4)= [MEAN_HIT_cluster2.rt_T2;NaN*(length(MEAN_HIT_cluster1.rt_T2)-length(MEAN_HIT_cluster2.rt_T2))]*1000; 


%------------------------------------------------------------------------%
%                   Plot Dprime Clustered data;
%------------------------------------------------------------------------%

%select the columns of dprimes;
dat = dprimes_clustered;          
col = [0 0 0];                        
dotsize = 15;                              
smooth = 0;    

%Generate your plot by calling the function;
close all;
vert_rain_plot(dat,col,dotsize,smooth)
pause(0.01);

%configure axis;
set(gcf,'color','w');
ax = gca;
ax.Box = 'off';
ax.LineWidth = 2;
ax.TickDir = 'out';
ax.XTick = [1 1.5 2 3 3.5 4];
labels = {'\fontsize{12}Pop.1','\fontsize{17}  \fontsize{14}\bfTone T1','\fontsize{12}Pop.2','\fontsize{12}Pop.1','\fontsize{17}  \fontsize{14}\bfTone T2','\fontsize{12}Pop.2'};
labels = cellfun(@(x) strrep(x,'  ','\newline'), labels,'UniformOutput',false);
ax.XTickLabel = labels;
ax.XAxis.TickLength = [0 0];
ax.XRuler.Axle.LineStyle = 'none'; 
ax.XAxis.FontSize = 14;
ax.YAxis.Limits = [1.5 7];
ax.YTick = [1.5 3.5 5.5 7];
ax.YAxis.FontSize = 12;
ax.YAxis.FontWeight = 'bold';
ax.YColor = [0 0 0];
ylabel('Dprime','FontSize',14, 'Layer','front', 'FontWeight', 'bold');
%add significance;
hold on 
x1 = linspace(3.1,3.7);
y1 = linspace(6.3,6.3);
plot(x1,y1, 'Color', [0 0 0],'LineWidth',0.5);
text(3.3,6.4,'*','FontSize', 20,'FontWeight','light');
hold on 
x1 = linspace(1.5,3.5);
y1 = linspace(6.8,6.8);
plot(x1,y1, 'Color', [0 0 0],'LineWidth',0.5);
text(2.4,6.9,'*','FontSize', 20,'FontWeight','light');
set(gcf,'color','w','Positio',[727 375 432 517]);

%------------------------------------------------------------------------%
%                Plot Reaction Times Clustered data;
%------------------------------------------------------------------------%

%select the columns of reaction times;
dat = rt_clustered;          
col = [0 0 0];                        
dotsize = 15;                              
smooth = 0;    

%Generate your plot by calling the function;
close all;
vert_rain_plot(dat,col,dotsize,smooth)
pause(0.01);

set(gcf,'color','w');
ax = gca;
ax.Box = 'off';
ax.LineWidth = 2;
ax.TickDir = 'out';
ax.XTick = [1 1.5 2 3 3.5 4];
labels = {'\fontsize{12}Pop.1','\fontsize{17}  \fontsize{14}\bfTone T1','\fontsize{12}Pop.2','\fontsize{12}Pop.1','\fontsize{17}  \fontsize{14}\bfTone T2','\fontsize{12}Pop.2'};
labels = cellfun(@(x) strrep(x,'  ','\newline'), labels,'UniformOutput',false);
ax.XTickLabel = labels;
ax.XAxis.TickLength = [0 0];
ax.XRuler.Axle.LineStyle = 'none'; 
ax.XAxis.FontSize = 14;
ax.YAxis.Limits = [300 700];
ax.YTick = [300 400 500 600 700];
ax.YAxis.FontSize = 12;
ax.YAxis.FontWeight = 'bold';
ax.YColor = [0 0 0];
ylabel('Reaction Times (ms)','FontSize',14, 'Layer','front', 'FontWeight', 'bold');
%add significance;
hold on 
x1 = linspace(3.1,3.7);
y1 = linspace(6.3,6.3);
plot(x1,y1, 'Color', [0 0 0],'LineWidth',0.5);
text(3.3,6.4,'*','FontSize', 20,'FontWeight','light');
hold on 
x1 = linspace(1.5,3.5);
y1 = linspace(670,670);
plot(x1,y1, 'Color', [0 0 0],'LineWidth',0.5);
text(2.4,675,'*','FontSize', 20,'FontWeight','light');
set(gcf,'color','w','Positio',[727 375 432 517]);

%% Plot Dprimes/RTs Clustered data IN ONE TONE CONDITION;

cluster1 = [2 3 5 6 9 10 11 13 15 19 22]; 
cluster2 = [7 8 12 14 17 20 21 23 24 29];

%First separate the data from cluster1 and cluster2 (because different
%lengths);
MEAN_HIT_cluster1 = OneTone_T1vsT2_performances(cluster1,:);
MEAN_HIT_cluster2 = OneTone_T1vsT2_performances(cluster2,:);
clear dprimes_clustered rt_clustered
dprimes_clustered(:,1)= MEAN_HIT_cluster1.dp_T1;
dprimes_clustered(:,2)= [MEAN_HIT_cluster2.dp_T1;NaN*(length(MEAN_HIT_cluster1.dp_T1)-length(MEAN_HIT_cluster2.dp_T1))]; 
dprimes_clustered(:,3)= MEAN_HIT_cluster1.dp_T2;
dprimes_clustered(:,4)= [MEAN_HIT_cluster2.dp_T2;NaN*(length(MEAN_HIT_cluster1.dp_T2)-length(MEAN_HIT_cluster2.dp_T2))];     
rt_clustered(:,1)= MEAN_HIT_cluster1.rt_T1*1000;
rt_clustered(:,2)= [MEAN_HIT_cluster2.rt_T1;NaN*(length(MEAN_HIT_cluster1.rt_T1)-length(MEAN_HIT_cluster2.rt_T1))]*1000; 
rt_clustered(:,3)= MEAN_HIT_cluster1.rt_T2*1000;
rt_clustered(:,4)= [MEAN_HIT_cluster2.rt_T2;NaN*(length(MEAN_HIT_cluster1.rt_T2)-length(MEAN_HIT_cluster2.rt_T2))]*1000; 


%------------------------------------------------------------------------%
%                   Plot Dprime Clustered data;
%------------------------------------------------------------------------%

%select the columns of dprimes;
dat = dprimes_clustered;          
col = [0 0 0];                        
dotsize = 15;                              
smooth = 0;    

%Generate your plot by calling the function;
close all;
vert_rain_plot(dat,col,dotsize,smooth)
pause(0.01);

%configure axis;
set(gcf,'color','w');
ax = gca;
ax.Box = 'off';
ax.LineWidth = 2;
ax.TickDir = 'out';
ax.XTick = [1 1.5 2 3 3.5 4];
labels = {'\fontsize{12}Pop.1','\fontsize{17}  \fontsize{14}\bfTone T1','\fontsize{12}Pop.2','\fontsize{12}Pop.1','\fontsize{17}  \fontsize{14}\bfTone T2','\fontsize{12}Pop.2'};
labels = cellfun(@(x) strrep(x,'  ','\newline'), labels,'UniformOutput',false);
ax.XTickLabel = labels;
ax.XAxis.TickLength = [0 0];
ax.XRuler.Axle.LineStyle = 'none'; 
ax.XAxis.FontSize = 14;
ax.YAxis.Limits = [1.5 7];
ax.YTick = [1.5 3.5 5.5 7];
ax.YAxis.FontSize = 12;
ax.YAxis.FontWeight = 'bold';
ax.YColor = [0 0 0];
ylabel('Dprime','FontSize',14, 'Layer','front', 'FontWeight', 'bold');
%add significance;
hold on 
x1 = linspace(3.1,3.7);
y1 = linspace(6.3,6.3);
plot(x1,y1, 'Color', [0 0 0],'LineWidth',0.5);
text(3.3,6.4,'*','FontSize', 20,'FontWeight','light');
hold on 
x1 = linspace(1.5,3.5);
y1 = linspace(6.8,6.8);
plot(x1,y1, 'Color', [0 0 0],'LineWidth',0.5);
text(2.4,6.9,'*','FontSize', 20,'FontWeight','light');
set(gcf,'color','w','Positio',[727 375 432 517]);


%------------------------------------------------------------------------%
%                Plot Reaction Times Clustered data;
%------------------------------------------------------------------------%

%select the columns of reaction times;
dat = rt_clustered;          
col = [0 0 0];                        
dotsize = 15;                              
smooth = 0;    

%Generate your plot by calling the function;
close all;
vert_rain_plot(dat,col,dotsize,smooth)
pause(0.01);

set(gcf,'color','w');
ax = gca;
ax.Box = 'off';
ax.LineWidth = 2;
ax.TickDir = 'out';
ax.XTick = [1 1.5 2 3 3.5 4];
labels = {'\fontsize{12}Pop.1','\fontsize{17}  \fontsize{14}\bfTone T1','\fontsize{12}Pop.2','\fontsize{12}Pop.1','\fontsize{17}  \fontsize{14}\bfTone T2','\fontsize{12}Pop.2'};
labels = cellfun(@(x) strrep(x,'  ','\newline'), labels,'UniformOutput',false);
ax.XTickLabel = labels;
ax.XAxis.TickLength = [0 0];
ax.XRuler.Axle.LineStyle = 'none'; 
ax.XAxis.FontSize = 14;
ax.YAxis.Limits = [300 700];
ax.YTick = [300 400 500 600 700];
ax.YAxis.FontSize = 12;
ax.YAxis.FontWeight = 'bold';
ax.YColor = [0 0 0];
ylabel('Reaction Times (ms)','FontSize',14, 'Layer','front', 'FontWeight', 'bold');
%add significance;
hold on 
x1 = linspace(3.1,3.7);
y1 = linspace(6.3,6.3);
plot(x1,y1, 'Color', [0 0 0],'LineWidth',0.5);
text(3.3,6.4,'*','FontSize', 20,'FontWeight','light');
hold on 
x1 = linspace(1.5,3.5);
y1 = linspace(670,670);
plot(x1,y1, 'Color', [0 0 0],'LineWidth',0.5);
text(2.4,675,'*','FontSize', 20,'FontWeight','light');
set(gcf,'color','w','Positio',[727 375 432 517]);


%%