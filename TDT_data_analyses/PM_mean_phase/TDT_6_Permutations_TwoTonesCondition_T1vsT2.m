%% Compare the resultant vector length between T1 and T2 Tones to test whether the mean entrainment is significantly greater in at tone T2 in TwoTone ondition;
% We'll do this with a permutation test on z-value differences from Rayleigh tests. 
clearvars;clc; 
addpath C:\toolbox\CircHist-master; addpath C:\toolbox\Circular_Statistics_Toolbox\; cd C:\ANALYSES_DATA\;

%subjects ID;
subjects = 1:20;

%new ranperm for each participant;
rng('Shuffle');

%Chose the number of permutations;
permutations = 10000;

%Compute the 10000 permutations;
clear perm_hit_2Tones_first perm_hit_2Tones_second 
for s = subjects

    cd(['C:\ANALYSES_DATA\subj' num2str(s)]);
    load(['subj' num2str(s) '-Phase_analyses']);
    
    %keep only TwoTones condition;
    new_phase_analysis = PHASE_ANALYSES;    
    new_phase_analysis(new_phase_analysis.CONDITION == 1,:) = [];
    %remove T2 occurring in the T1 time-window;
    new_phase_analysis(new_phase_analysis.TONE_POSITION==2,:) = [];
    
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
% save permutations_10000_TwoTones_T1vsT2 perm_hit_2Tones_first perm_hit_2Tones_second

%% Compare statistically the entrainment between fT1 and T2 tones in Two Tones condition: Z-values from Rayleigh's tests;
clearvars; clc;
addpath C:\toolbox\CircHist-master; addpath C:\toolbox\Circular_Statistics_Toolbox\; cd C:\ANALYSES_DATA\scripts_Github\;
load permutations_10000_TwoTones_T1vsT2;

%subjects ID;
subjects = 1:29;

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
    
    cd(['C:\ANALYSES_DATA\subj' num2str(s)]);
    load(['subj' num2str(s) '-Phase_analyses']) 

    ori_hit_2Tones_first(s,:) = [cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==1)); NaN(100 - length(cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==1))),1)]';
    ori_hit_2Tones_second(s,:) = [cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==2 & PHASE_ANALYSES.ONSET/44100 > 2.5)); NaN(100 - length(cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==2 & PHASE_ANALYSES.ONSET/44100 > 2.5))),1)]';
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
P_value_perm = (sum(diff_Z_T2minusT1(:,1)> diff_Z_orig)+1)/((length(diff_Z_T2minusT1(:,1)))+1);

display(['P_value_perm: ' num2str(P_value_perm)]);
display(['Effect size: ' num2str(circ_r(mean_ori_hit_2Tones_second(:,1))- circ_r(mean_ori_hit_2Tones_first(:,1)))]); 

%%
