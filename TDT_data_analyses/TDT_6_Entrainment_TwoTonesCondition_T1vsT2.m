%%PHASE ANALYSES WITHOUT REALIGNING;
clearvars;clc;
addpath XXX\CircHist-master; addpath XXX\Circular_Statistics_Toolbox\;
cd XXX;

%subjects ID;
subjects = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29];

% I- Compile the hit phases of all participants together;
hit_angle_GA = zeros(length(subjects),300);
miss_angle_GA = zeros(length(subjects),300);
hit_angle_1Tone = zeros(length(subjects),100);
miss_angle_1Tone = zeros(length(subjects),100);
hit_angle_2Tones_first = zeros(length(subjects),100);
miss_angle_2Tones_first = zeros(length(subjects),100);
hit_angle_2Tones_second = zeros(length(subjects),100);
miss_angle_2Tones_second = zeros(length(subjects),100);

for s = subjects
    
    cd(['XXX\subj' num2str(s)]);
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


%% Plot the mean theta phase of subjects at T1 and T2 Tones in the Two Tones conditions (hits only;
    
%Calculate the mean angle for each participant;
clear mean_hit_T1 mean_hit_T2 

for s = subjects 

mean_hit_T1(s,1) = circ_mean(hit_angle_2Tones_first(s,~isnan(hit_angle_2Tones_first(s,:)))');
mean_hit_T1(s,2) = s;      
mean_hit_T2(s,1) = circ_mean(hit_angle_2Tones_second(s,~isnan(hit_angle_2Tones_second(s,:)))');
mean_hit_T2(s,2) = s;

end

%Plot the resultant vector at T1 and T2 Tones in the two Tones conditions;

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

%Rayleigh's Test on the T1 and T2 mean phase distributions; 
[r_hit_T1, z_hit_T1] = circ_rtest(mean_hit_T1(:,1));
[r_hit_T2, z_hit_T2] = circ_rtest(mean_hit_T2(:,1));

%% Now we want to compare the resultant vector length between T1 and T2 Tones to test whether the mean entrainment is significantly greater in at tone T2 in TwoTone ondition;
%We'll do this with a permutation test on z-value differences from Rayleigh tests. 
clearvars;clc; 
addpath XXX\CircHist-master; addpath XXX\Circular_Statistics_Toolbox\; cd XXX\;

%subjects ID;
subjects = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29];

%new ranperm for each participant;
rng('Shuffle');

%Chose the number of permutations;
permutations = 10000;

%Compute the 10000 permutations;
clear perm_hit_2Tones_first perm_hit_2Tones_second 
for s = subjects

    cd(['XXX\subj' num2str(s)]);
    load(['subj' num2str(s) '-Phase_analyses']);
    
    %keep only TwoTones condition;
    new_phase_analysis = PHASE_ANALYSES;    
    new_phase_analysis(new_phase_analysis.CONDITION == 1,:) = [];
    
    for j = 1:permutations
        
        new_conlab = randperm(size(new_phase_analysis,1));
        temp_position = new_phase_analysis.TONE_POSITION;
        new_phase_analysis.TONE_POSITION = temp_position(new_conlab);
        
        %match the number of trials in the subsamples for each permutation;
        if length(cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 1 & new_phase_analysis.CONDITION== 2 & new_phase_analysis.TONE_POSITION== 1))) < length(cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 1 & new_phase_analysis.CONDITION== 2 & new_phase_analysis.TONE_POSITION== 2)))
             perm_hit_2Tones_first(s,j,:) = [datasample(cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 1 & new_phase_analysis.CONDITION== 2 & new_phase_analysis.TONE_POSITION== 1)),length(cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 1 & new_phase_analysis.CONDITION== 2 & new_phase_analysis.TONE_POSITION== 1))),'Replace',false); NaN(300-length(cell2mat(new_phase_analysis.THETA_PHASE((new_phase_analysis.HIT== 1 & new_phase_analysis.CONDITION== 2 & new_phase_analysis.TONE_POSITION== 1)))),1)]';
             perm_hit_2Tones_second(s,j,:) = [datasample(cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 1 & new_phase_analysis.CONDITION== 2 & new_phase_analysis.TONE_POSITION== 2)),length(cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 1 & new_phase_analysis.CONDITION== 2 & new_phase_analysis.TONE_POSITION== 1))),'Replace',false); NaN(300-length(cell2mat(new_phase_analysis.THETA_PHASE((new_phase_analysis.HIT== 1 & new_phase_analysis.CONDITION== 2 & new_phase_analysis.TONE_POSITION== 1)))),1)]';

        elseif length(cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 1 & new_phase_analysis.CONDITION== 2 & new_phase_analysis.TONE_POSITION== 1))) > length(cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 1 & new_phase_analysis.CONDITION== 2 & new_phase_analysis.TONE_POSITION== 2)))
             perm_hit_2Tones_first(s,j,:) = [datasample(cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 1 & new_phase_analysis.CONDITION== 2 & new_phase_analysis.TONE_POSITION== 1)),length(cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 1 & new_phase_analysis.CONDITION== 2 & new_phase_analysis.TONE_POSITION== 2))),'Replace',false); NaN(300-length(cell2mat(new_phase_analysis.THETA_PHASE((new_phase_analysis.HIT== 1 & new_phase_analysis.CONDITION== 2 & new_phase_analysis.TONE_POSITION== 2)))),1)]';
             perm_hit_2Tones_second(s,j,:) = [datasample(cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 1 & new_phase_analysis.CONDITION== 2 & new_phase_analysis.TONE_POSITION== 2)),length(cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 1 & new_phase_analysis.CONDITION== 2 & new_phase_analysis.TONE_POSITION== 2))),'Replace',false); NaN(300-length(cell2mat(new_phase_analysis.THETA_PHASE((new_phase_analysis.HIT== 1 & new_phase_analysis.CONDITION== 2 & new_phase_analysis.TONE_POSITION== 2)))),1)]';
        end  

    end       
        
end

% cd XXX\;
% save permutations_10000_2Tones_T1vsT2 perm_hit_2Tones_first perm_hit_2Tones_second

%% Compare statistically the entrainment between fT1 and T2 tones in Two Tones condition: Z-values from Rayleigh's tests;
clearvars; clc;
addpath XXX\CircHist-master; addpath XXX\Circular_Statistics_Toolbox\; cd XXX\;
load permutations_10000_2Tones_T1vsT2;

%subjects ID;
subjects = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29];

%Calculate the mean phase T1/T2 hits of each of the 10000 permutations in every participants;
clear mean_hit_2Tones_first mean_hit_2Tones_second 

for s = subjects
try
    
    for j = 1:size(perm_hit_2Tones_first(:,:,:),2)

        mean_hit_2Tones_first(s,j) = circ_mean(permute(perm_hit_2Tones_first(s,j,~isnan(perm_hit_2Tones_first(s,j,:))),[3 1 2]));
        mean_hit_2Tones_second(s,j) = circ_mean(permute(perm_hit_2Tones_second(s,j,~isnan(perm_hit_2Tones_second(s,j,:))),[3 1 2]));
    
    end
    
catch
end
    
end

%Remove the bad subjects to exclude from statistic analyses;
bad_subjects = [1 4 25 26 27];
mean_hit_2Tones_first(bad_subjects,:) = [];
mean_hit_2Tones_second(bad_subjects,:) = [];

%now, for each mean phase permuted, perform a Rayleigh's test to get the resultant vector lengths of T1/T2 hits(z_value);
clear Z_2Tones_first Z_2Tones_second    

for t = 1:length(mean_hit_2Tones_second(1,:)) 
    
    [p_mean_hit_2Tones_first, z_2Tones_first] = circ_rtest(mean_hit_2Tones_first(:,t));
    [p_mean_hit_2Tones_second, z_2Tones_second] = circ_rtest(mean_hit_2Tones_second(:,t));
    
    Z_2Tones_first(t,1)= p_mean_hit_2Tones_first;
    Z_2Tones_second(t,1)= p_mean_hit_2Tones_second;
    Z_2Tones_first(t,2)= z_2Tones_first;
    Z_2Tones_second(t,2)= z_2Tones_second;
     
    clear p_mean_hit_2Tones_second z_2Tones_second
  
end  
    
%Now calculate the mean phase on the original data set to get real Zvalue difference between T1/T2 Tones;
clear ori_hit_2Tones_first ori_hit_2Tones_second mean_ori_hit_2Tones_first mean_ori_hit_2Tones_second

for s = subjects
    
    cd(['XXX\subj' num2str(s)]);
    load(['subj' num2str(s) '-Phase_analyses']) 

    ori_hit_2Tones_first(s,:) = [cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==1)); NaN(100 - length(cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==1))),1)]';
    ori_hit_2Tones_second(s,:) = [cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==2)); NaN(100 - length(cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==2))),1)]';
    mean_ori_hit_2Tones_first(s,1) = circ_mean(ori_hit_2Tones_first(s,~isnan(ori_hit_2Tones_first(s,:)))');
    mean_ori_hit_2Tones_second(s,1) = circ_mean(ori_hit_2Tones_second(s,~isnan(ori_hit_2Tones_second(s,:)))');
    
    clear PHASE_ANALYSES alltrials participant_parameters

end
   
%Remove the bad subjects to exclude from statistic analyses;
bad_subjects = [1 4 25 26 27];
mean_ori_hit_2Tones_first(bad_subjects,:) = [];
mean_ori_hit_2Tones_second(bad_subjects,:) = [];   

[~, ori_z_2Tones_first] = circ_rtest(mean_ori_hit_2Tones_first(:,1));
[~, ori_z_2Tones_second] = circ_rtest(mean_ori_hit_2Tones_second(:,1));
       
%Calculate the difference of vector lengths (z_values) T2-T1 hits for each permutation and for the original data;
%Then sort in descend order;
diff_Z_T2minusT1(:,1) = Z_2Tones_second(:,2) - Z_2Tones_first(:,2);
diff_Z_T2minusT1(:,1)= sortrows(diff_Z_T2minusT1(:,1),'descend');
diff_Z_orig = ori_z_2Tones_second - ori_z_2Tones_first;
P_value_perm = sum(diff_Z_T2minusT1(:,1)> diff_Z_orig)/(length(diff_Z_T2minusT1(:,1)));

display(['P_value_perm: ' num2str(P_value_perm)]);
display(['Effect size: ' num2str(circ_r(mean_ori_hit_2Tones_second(:,1))- circ_r(mean_ori_hit_2Tones_first(:,1)))]); 

%%
