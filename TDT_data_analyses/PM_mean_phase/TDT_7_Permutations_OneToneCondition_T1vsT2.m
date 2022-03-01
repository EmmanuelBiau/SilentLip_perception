%%Generate N-permutations of hit trials from first and second tone of OneTone condition; 
clearvars;clc; 
addpath C:\toolbox\CircHist-master; addpath C:\toolbox\Circular_Statistics_Toolbox\; cd C:\ANALYSES_DATA\;
load ONE_TONE_ANALYSES; %Where we sorted the single tones according to their onset to attribute them to T1 or T2 equivalent window;

%subjects ID;
subjects = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29];

%new ranperm for each participant;
rng('Shuffle');

%Chose the number of permutations;
permutations = 10000;


clear perm_hit_1Tone_T1 perm_hit_1Tone_T2
for s = subjects

    new_phase_analysis = ONE_TONE_ANALYSES(ONE_TONE_ANALYSES.PARTICIPANT == s,:);   
            
    for j = 1:permutations
        
        new_conlab = randperm(size(new_phase_analysis,1));
        temp_position = new_phase_analysis.OneTone_position;
        new_phase_analysis.OneTone_position = temp_position(new_conlab);

        if length(new_phase_analysis.OneTone_phase(new_phase_analysis.Hit== 1 & new_phase_analysis.OneTone_position== 1)) < length(new_phase_analysis.OneTone_phase(new_phase_analysis.Hit== 1 & new_phase_analysis.OneTone_position== 2))
             perm_hit_1Tone_T1(s,j,:) = [datasample(new_phase_analysis.OneTone_phase(new_phase_analysis.Hit== 1 & new_phase_analysis.OneTone_position== 1),length(new_phase_analysis.OneTone_phase(new_phase_analysis.Hit== 1 & new_phase_analysis.OneTone_position== 1)),'Replace',false); NaN(300-length(new_phase_analysis.OneTone_phase((new_phase_analysis.Hit== 1 & new_phase_analysis.OneTone_position== 1))),1)]';
             perm_hit_1Tone_T2(s,j,:) = [datasample(new_phase_analysis.OneTone_phase(new_phase_analysis.Hit== 1 & new_phase_analysis.OneTone_position== 2),length(new_phase_analysis.OneTone_phase(new_phase_analysis.Hit== 1 & new_phase_analysis.OneTone_position== 1)),'Replace',false); NaN(300-length(new_phase_analysis.OneTone_phase((new_phase_analysis.Hit== 1 & new_phase_analysis.OneTone_position== 1))),1)]';

        elseif length(new_phase_analysis.OneTone_phase(new_phase_analysis.Hit== 1 & new_phase_analysis.OneTone_position== 1)) > length(new_phase_analysis.OneTone_phase(new_phase_analysis.Hit== 1 & new_phase_analysis.OneTone_position== 2))
             perm_hit_1Tone_T1(s,j,:) = [datasample(new_phase_analysis.OneTone_phase(new_phase_analysis.Hit== 1 & new_phase_analysis.OneTone_position== 1),length(new_phase_analysis.OneTone_phase(new_phase_analysis.Hit== 1 & new_phase_analysis.OneTone_position== 2)),'Replace',false); NaN(300-length(new_phase_analysis.OneTone_phase((new_phase_analysis.Hit== 1 & new_phase_analysis.OneTone_position== 2))),1)]';
             perm_hit_1Tone_T2(s,j,:) = [datasample(new_phase_analysis.OneTone_phase(new_phase_analysis.Hit== 1 & new_phase_analysis.OneTone_position== 2),length(new_phase_analysis.OneTone_phase(new_phase_analysis.Hit== 1 & new_phase_analysis.OneTone_position== 2)),'Replace',false); NaN(300-length(new_phase_analysis.OneTone_phase((new_phase_analysis.Hit== 1 & new_phase_analysis.OneTone_position== 2))),1)]';
        end  

    end       
        
end

% cd XXX\;
% save permutations_10000_OneTone_T1vsT2 perm_hit_1Tone_T1 perm_hit_1Tone_T2

%% Compare statistically the entrainment between first and second tone conditions: Z-values from Rayleigh's tests;
clearvars;clc; 
addpath C:\toolbox\CircHist-master; addpath C:\toolbox\Circular_Statistics_Toolbox\; cd C:\ANALYSES_DATA\scripts_Github\;
load permutations_10000_OneTone_T1vsT2; load ONE_TONE_ANALYSES;

%subjects ID;
subjects = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29];

%Calculate the mean phase T1/T2 hits of each of the 10000 permutations in every participants;
clear mean_hit_2Tones_first mean_hit_2Tones_second 

for s = subjects
try
    
    for j = 1:size(perm_hit_1Tone_T1(:,:,:),2)

        mean_hit_1Tone_T1(s,j) = circ_mean(permute(perm_hit_1Tone_T1(s,j,~isnan(perm_hit_1Tone_T1(s,j,:))),[3 1 2]));
        mean_hit_1Tone_T2(s,j) = circ_mean(permute(perm_hit_1Tone_T2(s,j,~isnan(perm_hit_1Tone_T2(s,j,:))),[3 1 2]));
    
    end
    
catch
end
    
end

%Remove the bad subjects to exclude from statistic analyses;
bad_subjects = [1 4 25 26 27];
mean_hit_1Tone_T1(bad_subjects,:) = [];
mean_hit_1Tone_T2(bad_subjects,:) = [];

%now, for each mean phase permuted, perform a Rayleigh's test to get the resultant vector lengths of T1/T2 hits(z_value);
clear Z_2Tones_first Z_2Tones_second    

for t = 1:length(mean_hit_1Tone_T2(1,:)) 
    
    [p_mean_hit_1Tone_T1, z_1Tone_T1] = circ_rtest(mean_hit_1Tone_T1(:,t));
    [p_mean_hit_1Tone_T2, z_1Tone_T2] = circ_rtest(mean_hit_1Tone_T2(:,t));
    
    Z_1Tone_T1(t,1)= p_mean_hit_1Tone_T1;
    Z_1Tone_T2(t,1)= p_mean_hit_1Tone_T2;
    Z_1Tone_T1(t,2)= z_1Tone_T1;
    Z_1Tone_T2(t,2)= z_1Tone_T2;
     
    clear p_mean_hit_1Tone_T1 z_1Tone_T2
  
end  
    
%Now calculate the mean phase on the original data set to get real Zvalue difference between T1/T2 hits;
clear ori_hit_2Tones_first ori_hit_2Tones_second mean_ori_hit_2Tones_first mean_ori_hit_2Tones_second

for s = subjects
    
    ori_phase_analysis = ONE_TONE_ANALYSES(ONE_TONE_ANALYSES.PARTICIPANT == s,:); 

    ori_hit_1Tone_T1(s,:) = [ori_phase_analysis.OneTone_phase(ori_phase_analysis.Hit==1 & ori_phase_analysis.OneTone_position==1); NaN(100 - length(ori_phase_analysis.OneTone_phase(ori_phase_analysis.Hit==1 & ori_phase_analysis.OneTone_position==1)),1)]';
    ori_hit_1Tone_T2(s,:) = [ori_phase_analysis.OneTone_phase(ori_phase_analysis.Hit==1 & ori_phase_analysis.OneTone_position==2); NaN(100 - length(ori_phase_analysis.OneTone_phase(ori_phase_analysis.Hit==1 & ori_phase_analysis.OneTone_position==2)),1)]';
    mean_ori_hit_1Tone_T1(s,1) = circ_mean(ori_hit_1Tone_T1(s,~isnan(ori_hit_1Tone_T1(s,:)))');
    mean_ori_hit_1Tone_T2(s,1) = circ_mean(ori_hit_1Tone_T2(s,~isnan(ori_hit_1Tone_T2(s,:)))');
    
    clear PHASE_ANALYSES alltrials participant_parameters

end
   
%Remove the bad subjects to exclude from statistic analyses;
bad_subjects = [1 4 25 26 27];
mean_ori_hit_1Tone_T1(bad_subjects,:) = [];
mean_ori_hit_1Tone_T2(bad_subjects,:) = [];   

[~, ori_z_2Tones_first] = circ_rtest(mean_ori_hit_1Tone_T1(:,1));
[~, ori_z_2Tones_second] = circ_rtest(mean_ori_hit_1Tone_T2(:,1));
       
%Calculate the difference of vector lengths (z_values) T2-T1 hits for each permutation and for the original data;
%Then sort in descend order;
diff_Z_T2minusT1(:,1) = Z_1Tone_T2(:,2) - Z_1Tone_T1(:,2);
diff_Z_T2minusT1(:,1)= sortrows(diff_Z_T2minusT1(:,1),'descend');
diff_Z_orig = ori_z_2Tones_second - ori_z_2Tones_first;
P_value_perm = (sum(diff_Z_T2minusT1(:,1)> diff_Z_orig)+1)/((length(diff_Z_T2minusT1(:,1)))+1);

display(['P_value_perm: ' num2str(P_value_perm)]);
display(['Effect size: ' num2str(circ_r(mean_ori_hit_1Tone_T2(:,1))- circ_r(mean_ori_hit_1Tone_T1(:,1)))]);
  
%%