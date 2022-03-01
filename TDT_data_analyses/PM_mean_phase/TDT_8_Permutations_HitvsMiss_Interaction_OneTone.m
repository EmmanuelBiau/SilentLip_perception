%% Calculate the difference of Vector length hit-miss in T1 and T2 tones windows of the original data (single tone condition);
clearvars;clc; 
addpath C:\toolbox\CircHist-master; addpath C:\toolbox\Circular_Statistics_Toolbox\; cd C:\ANALYSES_DATA\;

%subjects ID;
subjects = 1:29;

%Calculate the mean phase on the original data set to get real Zvalue difference between hits and miss;
for s = subjects
    
    cd(['C:\ANALYSES_DATA\subj' num2str(s)]);
    load(['subj' num2str(s) '-Phase_analyses']);
    
    %sort single tones in T1 or T2;
    PHASE_ANALYSES(PHASE_ANALYSES.CONDITION ==2,:) = []; 
    t1 = PHASE_ANALYSES.ONSET/44100 <= 2.5;
    t2 = PHASE_ANALYSES.ONSET/44100 > 2.5;    
    PHASE_ANALYSES.TONE_POSITION(t1,:) = 1;
    PHASE_ANALYSES.TONE_POSITION(t2,:) = 2;
                 
    %T1 window of OneTone Condition; 
    ori_hit_T1_first(s,:) = [cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.TONE_POSITION==1)); NaN(100 - length(cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.TONE_POSITION==1))),1)]';
    ori_miss_T1_first(s,:) = [cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==0 & PHASE_ANALYSES.TONE_POSITION==1)); NaN(100 - length(cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==0 & PHASE_ANALYSES.TONE_POSITION==1))),1)]';
    ori_mean_hit_T1_first(s,1) = circ_mean(ori_hit_T1_first(s,~isnan(ori_hit_T1_first(s,:)))');
    ori_mean_miss_T1_first(s,1) = circ_mean(ori_miss_T1_first(s,~isnan(ori_miss_T1_first(s,:)))');

    %T2 window of OneTone Condition; 
    ori_hit_T1_second(s,:) = [cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.TONE_POSITION==2)); NaN(100 - length(cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.TONE_POSITION==2))),1)]';
    ori_miss_T1_second(s,:) = [cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==0 & PHASE_ANALYSES.TONE_POSITION==2)); NaN(100 - length(cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==0 & PHASE_ANALYSES.TONE_POSITION==2))),1)]';
    ori_mean_hit_T1_second(s,1) = circ_mean(ori_hit_T1_second(s,~isnan(ori_hit_T1_second(s,:)))');
    ori_mean_miss_T1_second(s,1) = circ_mean(ori_miss_T1_second(s,~isnan(ori_miss_T1_second(s,:)))');
            
    clear PHASE_ANALYSES alltrials participant_parameters
end


%Remove the bad subjects to exclude from statistic analyses;
bad_subjects = [1 4 25 26 27];
ori_mean_hit_T1_first(bad_subjects,:) = [];
ori_mean_miss_T1_first(bad_subjects,:) = [];
ori_mean_hit_T1_second(bad_subjects,:) = [];
ori_mean_miss_T1_second(bad_subjects,:) = [];


[~, ori_z_hit_T1_first] = circ_rtest(ori_mean_hit_T1_first(:,1));
[~, ori_z_miss_T1_first] = circ_rtest(ori_mean_miss_T1_first(:,1));
[~, ori_z_hit_T1_second] = circ_rtest(ori_mean_hit_T1_second(:,1));
[~, ori_z_miss_T1_second] = circ_rtest(ori_mean_miss_T1_second(:,1));

diff_Z_ori_T1 = ori_z_hit_T1_first - ori_z_miss_T1_first;
diff_Z_ori_T2 = ori_z_hit_T1_second - ori_z_miss_T1_second;

%Get the original difference of size effect between T2 and T1 tones from the two tones condition;
diff_Z_T2vsT1_ori = diff_Z_ori_T2 - diff_Z_ori_T1;

%% Generate the permuted data by shuffling the tone position condition; 
cd C:\ANALYSES_DATA\scripts_Github\;

%subject ID;
subjects = 1:29;

%new ranperm for each participant;
rng('shuffle');

%chose the number of permutations;
permutations = 10000;

%Generate the permutations in the 3 conditions;
clear perm_hit_2Tones_first perm_miss_2Tones_first perm_hit_2Tones_second perm_miss_2Tones_second

for s = subjects

    cd(['C:\ANALYSES_DATA\subj' num2str(s)]);
    load(['subj' num2str(s) '-Phase_analyses']) 
    
    %sort single tones in T1 or T2;
    PHASE_ANALYSES(PHASE_ANALYSES.CONDITION ==2,:) = []; 
    t1 = PHASE_ANALYSES.ONSET/44100 <= 2.5;
    t2 = PHASE_ANALYSES.ONSET/44100 > 2.5;    
    PHASE_ANALYSES.TONE_POSITION(t1,:) = 1;
    PHASE_ANALYSES.TONE_POSITION(t2,:) = 2;
    
    new_phase_analysis = PHASE_ANALYSES;
       
    for j = 1:permutations
        
    %Shuffle the hit/miss labels in each conditions;
    new_conlab1 = randperm(size(new_phase_analysis,1));  
    temp_hit1 = new_phase_analysis.TONE_POSITION;
    new_phase_analysis.TONE_POSITION = temp_hit1(new_conlab1);
    
        %Generate the permutations in the T1 window of TwoTones Condition;    
        if length(cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 0 & new_phase_analysis.CONDITION== 1 & new_phase_analysis.TONE_POSITION== 1))) < length(cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 1 & new_phase_analysis.CONDITION== 1 & new_phase_analysis.TONE_POSITION== 1))) 
             perm_hit_1Tone_first(s,j,:) = [datasample(cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 1 & new_phase_analysis.CONDITION== 1 & new_phase_analysis.TONE_POSITION== 1)),length(cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 0 & new_phase_analysis.CONDITION== 1 & new_phase_analysis.TONE_POSITION== 1))),'Replace',false); NaN(300-length(cell2mat(new_phase_analysis.THETA_PHASE((new_phase_analysis.HIT== 0 & new_phase_analysis.CONDITION== 1 & new_phase_analysis.TONE_POSITION== 1)))),1)]';
             perm_miss_1Tone_first(s,j,:) = [cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 0 & new_phase_analysis.CONDITION== 1 & new_phase_analysis.TONE_POSITION== 1)); NaN(300-length(cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 0 & new_phase_analysis.CONDITION== 1 & new_phase_analysis.TONE_POSITION== 1))),1)]';  
        elseif length(cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 0 & new_phase_analysis.CONDITION== 1 & new_phase_analysis.TONE_POSITION== 1))) > length(cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 1 & new_phase_analysis.CONDITION== 1 & new_phase_analysis.TONE_POSITION== 1)))
             perm_hit_1Tone_first(s,j,:) = [datasample(cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 1 & new_phase_analysis.CONDITION== 1 & new_phase_analysis.TONE_POSITION== 1)),length(cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 1 & new_phase_analysis.CONDITION== 1 & new_phase_analysis.TONE_POSITION== 1))),'Replace',false); NaN(300-length(cell2mat(new_phase_analysis.THETA_PHASE((new_phase_analysis.HIT== 1 & new_phase_analysis.CONDITION== 1 & new_phase_analysis.TONE_POSITION== 1)))),1)]';
             perm_miss_1Tone_first(s,j,:) = [cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 0 & new_phase_analysis.CONDITION== 1 & new_phase_analysis.TONE_POSITION== 1)); NaN(300-length(cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 0 & new_phase_analysis.CONDITION== 1 & new_phase_analysis.TONE_POSITION== 1))),1)]';
        end

        %Generate the permutations in the T2 window of TwoTones Condition;            
         if length(cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 0 & new_phase_analysis.CONDITION== 1 & new_phase_analysis.TONE_POSITION== 2))) < length(cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 1 & new_phase_analysis.CONDITION== 1 & new_phase_analysis.TONE_POSITION== 2))) 
             perm_hit_1Tone_second(s,j,:) = [datasample(cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 1 & new_phase_analysis.CONDITION== 1 & new_phase_analysis.TONE_POSITION== 2)),length(cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 0 & new_phase_analysis.CONDITION== 1 & new_phase_analysis.TONE_POSITION== 2))),'Replace',false); NaN(300-length(cell2mat(new_phase_analysis.THETA_PHASE((new_phase_analysis.HIT== 0 & new_phase_analysis.CONDITION== 1 & new_phase_analysis.TONE_POSITION== 2)))),1)]';
             perm_miss_1Tone_second(s,j,:) = [cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 0 & new_phase_analysis.CONDITION== 1 & new_phase_analysis.TONE_POSITION== 2)); NaN(300-length(cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 0 & new_phase_analysis.CONDITION== 1 & new_phase_analysis.TONE_POSITION== 2))),1)]';  
        elseif length(cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 0 & new_phase_analysis.CONDITION== 1 & new_phase_analysis.TONE_POSITION== 2))) > length(cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 1 & new_phase_analysis.CONDITION== 1 & new_phase_analysis.TONE_POSITION== 2)))
             perm_hit_1Tone_second(s,j,:) = [datasample(cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 1 & new_phase_analysis.CONDITION== 1 & new_phase_analysis.TONE_POSITION== 2)),length(cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 1 & new_phase_analysis.CONDITION== 1 & new_phase_analysis.TONE_POSITION== 2))),'Replace',false); NaN(300-length(cell2mat(new_phase_analysis.THETA_PHASE((new_phase_analysis.HIT== 1 & new_phase_analysis.CONDITION== 1 & new_phase_analysis.TONE_POSITION== 2)))),1)]';
             perm_miss_1Tone_second(s,j,:) = [cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 0 & new_phase_analysis.CONDITION== 1 & new_phase_analysis.TONE_POSITION== 2)); NaN(300-length(cell2mat(new_phase_analysis.THETA_PHASE(new_phase_analysis.HIT== 0 & new_phase_analysis.CONDITION== 1 & new_phase_analysis.TONE_POSITION== 2))),1)]';
         end
         
    end

    clear new_phase_analysis
    
end

% cd C:\ANALYSES_DATA\scripts_Github\;
% save permutations_10000_HitvsMiss_T1vsT2_OneTone perm_hit_1Tone_first perm_hit_1Tone_second perm_miss_1Tone_first perm_miss_1Tone_second

%% Now perform the permutation test between original difference of effect size and permuted data; 
clc; addpath C:\toolbox\CircHist-master; addpath C:\toolbox\Circular_Statistics_Toolbox\;
cd C:\ANALYSES_DATA\scripts_Github\;
load permutations_10000_HitvsMiss_T1vsT2_OneTone;

%subjects ID;
subjects = 1:29;

%Calculate the mean phase of each of the 10000 permutations in every participants;
clear mean_hits_1T_first mean_miss_1T_first mean_hits_1T_second mean_miss_1T_seconds

for s = subjects
try
    
    for j = 1:size(perm_hit_1Tone_first(:,:,:),2)
                
        %T1 window of TwoTones Condition; 
        mean_hits_1T_first(s,j) = circ_mean(permute(perm_hit_1Tone_first(s,j,~isnan(perm_hit_1Tone_first(s,j,:))),[3 1 2]));
        mean_miss_1T_first(s,j) = circ_mean(permute(perm_miss_1Tone_first(s,j,~isnan(perm_miss_1Tone_first(s,j,:))),[3 1 2]));

        %T2 window of TwoTones Condition; 
        mean_hits_1T_second(s,j) = circ_mean(permute(perm_hit_1Tone_second(s,j,~isnan(perm_hit_1Tone_second(s,j,:))),[3 1 2]));
        mean_miss_1T_second(s,j) = circ_mean(permute(perm_miss_1Tone_second(s,j,~isnan(perm_miss_1Tone_second(s,j,:))),[3 1 2]));
                
    end    
              
catch
end
    
end

%Remove the bad subjects to exclude from statistic analyses;
bad_subjects = [1 4 25 26 27];
mean_hits_1T_first(bad_subjects,:) = [];
mean_miss_1T_first(bad_subjects,:) = [];
mean_hits_1T_second(bad_subjects,:) = [];
mean_miss_1T_second(bad_subjects,:) = [];

%now, for each mean phase permuted, perform a Rayleigh's test to get the hit/miss vector lengths(z_values);
clear Z_hit Z_miss   

for t = 1:length(mean_hits_1T_first(1,:)) 
    
    [p_mean_hit_1T_first, z_hit_1T_first] = circ_rtest(mean_hits_1T_first(:,t));
    [p_mean_miss_1T_first, z_miss_1T_first] = circ_rtest(mean_miss_1T_first(:,t));
    [p_mean_hit_1T_second, z_hit_1T_second] = circ_rtest(mean_hits_1T_second(:,t));
    [p_mean_miss_1T_second, z_miss_1T_second] = circ_rtest(mean_miss_1T_second(:,t));
    
    Z_hit_first(t,1)= p_mean_hit_1T_first;
    Z_miss_first(t,1)= p_mean_miss_1T_first;
    Z_hit_first(t,2)= z_hit_1T_first;
    Z_miss_first(t,2)= z_miss_1T_first;
    
    Z_hit_second(t,1)= p_mean_hit_1T_second;
    Z_miss_second(t,1)= p_mean_miss_1T_second;
    Z_hit_second(t,2)= z_hit_1T_second;
    Z_miss_second(t,2)= z_miss_1T_second;
    
    diff_Z_perm_T1(t,1) = Z_hit_first(t,2) - Z_miss_first(t,2); 
    diff_Z_perm_T2(t,1) = Z_hit_second(t,2) - Z_miss_second(t,2);
    
    clear p_mean_hit_2T_first z_hit_2T_first p_mean_miss_2T_first z_miss_2T_first
    clear p_mean_hit_2T_second z_hit_2T_second p_mean_miss_2T_second z_miss_2T_second
  
end  

%Calculate the difference of vector lengths (z_values) hit-miss for each permutation and for the original data;
%Then sort in descend order;
diff_Z_T2vsT1_perm = diff_Z_perm_T2 - diff_Z_perm_T1;
diff_Z_T2vsT1_perm(:,1)= sortrows(diff_Z_T2vsT1_perm,'descend');

%Now compare the original size effect to the permuted data and get p-value;
P_value_perm = (sum(diff_Z_T2vsT1_perm > diff_Z_T2vsT1_ori)+1)/(length(diff_Z_T2vsT1_perm)+1);
display(['P_value_perm: ' num2str(P_value_perm)]); 

%display the original effect size;
r_hit_T1_first = circ_r(ori_mean_hit_T1_first(:,1));
r_miss_T1_first = circ_r(ori_mean_miss_T1_first(:,1));
r_hit_T1_second = circ_r(ori_mean_hit_T1_second(:,1));
r_miss_T1_second = circ_r(ori_mean_miss_T1_second(:,1));
es_T1 = r_hit_T1_first - r_miss_T1_first;
es_T2 = r_hit_T1_second - r_miss_T1_second;
display(['Effect size: ' num2str(es_T2-es_T1)]); 

%%
