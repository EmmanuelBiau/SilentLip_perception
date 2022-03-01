%% GENERATE PERMUTATIONS ON PHASE DATA WITH SUBSAMPLING THE NUMBER OF HITS TO THE MISS TRIALS IN EACH CONDITIONS;
clearvars;clc;
addpath C:\toolbox\CircHist-master; addpath C:\toolbox\Circular_Statistics_Toolbox\; addpath D:\Birmingham_2018\BEHAVIOR_project\ANALYSES_DATA\;

%subjects ID;
subjects = 1:29;

%new ranperm for each participant;
rng('Shuffle');

%chose the number of permutations;
permutations = 10000;

%Generate the permutations in the 3 conditions;
clear perm_hit_GA perm_miss_GA perm_hit_1Tone perm_miss_1Tone perm_hit_2Tones_first perm_miss_2Tones_first perm_hit_2Tones_second perm_miss_2Tones_second

for s = subjects

    cd(['D:\Birmingham_2018\BEHAVIOR_project\ANALYSES_DATA\subj' num2str(s)]);
    load(['subj' num2str(s) '-Phase_analyses']) 
    new_phase_analysis_1Tone = PHASE_ANALYSES(PHASE_ANALYSES.CONDITION ==1,:);  
    new_phase_analysis_2Tones_first = PHASE_ANALYSES(PHASE_ANALYSES.CONDITION ==2 & PHASE_ANALYSES.TONE_POSITION ==1,:);
    new_phase_analysis_2Tones_second = PHASE_ANALYSES(PHASE_ANALYSES.CONDITION ==2 & PHASE_ANALYSES.TONE_POSITION ==2,:);
        
    for j = 1:permutations
        
        %Shuffle the hit/miss labels in each conditions;
        new_conlab1 = randperm(size(new_phase_analysis_1Tone,1));
        new_conlab2 = randperm(size(new_phase_analysis_2Tones_first,1));
        new_conlab3 = randperm(size(new_phase_analysis_2Tones_second,1));
        temp_hit1 = new_phase_analysis_1Tone.HIT;
        temp_hit2 = new_phase_analysis_2Tones_first.HIT;
        temp_hit3 = new_phase_analysis_2Tones_second.HIT;
        new_phase_analysis_1Tone.HIT = temp_hit1(new_conlab1);
        new_phase_analysis_2Tones_first.HIT = temp_hit2(new_conlab2);
        new_phase_analysis_2Tones_second.HIT = temp_hit3(new_conlab3);
        
        %Generate the permutations in the 1Tone Condition;        
        if length(cell2mat(new_phase_analysis_1Tone.THETA_PHASE(new_phase_analysis_1Tone.HIT== 0 & new_phase_analysis_1Tone.CONDITION== 1))) < length(cell2mat(new_phase_analysis_1Tone.THETA_PHASE(new_phase_analysis_1Tone.HIT== 1 & new_phase_analysis_1Tone.CONDITION== 1)))
             perm_hit_1Tone(s,j,:) = [datasample(cell2mat(new_phase_analysis_1Tone.THETA_PHASE(new_phase_analysis_1Tone.HIT==1 & new_phase_analysis_1Tone.CONDITION== 1)),length(cell2mat(new_phase_analysis_1Tone.THETA_PHASE(new_phase_analysis_1Tone.HIT==0 & new_phase_analysis_1Tone.CONDITION== 1))),'Replace',false); NaN(300-length(cell2mat(new_phase_analysis_1Tone.THETA_PHASE((new_phase_analysis_1Tone.HIT==0 & new_phase_analysis_1Tone.CONDITION== 1)))),1)]';
             perm_miss_1Tone(s,j,:) = [cell2mat(new_phase_analysis_1Tone.THETA_PHASE(new_phase_analysis_1Tone.HIT==0 & new_phase_analysis_1Tone.CONDITION== 1)); NaN(300-length(cell2mat(new_phase_analysis_1Tone.THETA_PHASE(new_phase_analysis_1Tone.HIT==0 & new_phase_analysis_1Tone.CONDITION== 1))),1)]';
        elseif length(cell2mat(new_phase_analysis_1Tone.THETA_PHASE(new_phase_analysis_1Tone.HIT== 0 & new_phase_analysis_1Tone.CONDITION== 1))) > length(cell2mat(new_phase_analysis_1Tone.THETA_PHASE(new_phase_analysis_1Tone.HIT== 1 & new_phase_analysis_1Tone.CONDITION== 1)))
             perm_hit_1Tone(s,j,:) = [datasample(cell2mat(new_phase_analysis_1Tone.THETA_PHASE(new_phase_analysis_1Tone.HIT==1 & new_phase_analysis_1Tone.CONDITION== 1)),length(cell2mat(new_phase_analysis_1Tone.THETA_PHASE(new_phase_analysis_1Tone.HIT==1 & new_phase_analysis_1Tone.CONDITION== 1))),'Replace',false); NaN(300-length(cell2mat(new_phase_analysis_1Tone.THETA_PHASE((new_phase_analysis_1Tone.HIT==1 & new_phase_analysis_1Tone.CONDITION== 1)))),1)]';
             perm_miss_1Tone(s,j,:) = [cell2mat(new_phase_analysis_1Tone.THETA_PHASE(new_phase_analysis_1Tone.HIT==0 & new_phase_analysis_1Tone.CONDITION== 1)); NaN(300-length(cell2mat(new_phase_analysis_1Tone.THETA_PHASE(new_phase_analysis_1Tone.HIT==0 & new_phase_analysis_1Tone.CONDITION== 1))),1)]';
        end 

        %Generate the permutations in the 2Tones-T1 Condition;    
        if length(cell2mat(new_phase_analysis_2Tones_first.THETA_PHASE(new_phase_analysis_2Tones_first.HIT== 0 & new_phase_analysis_2Tones_first.CONDITION== 2 & new_phase_analysis_2Tones_first.TONE_POSITION== 1))) < length(cell2mat(new_phase_analysis_2Tones_first.THETA_PHASE(new_phase_analysis_2Tones_first.HIT== 1 & new_phase_analysis_2Tones_first.CONDITION== 2 & new_phase_analysis_2Tones_first.TONE_POSITION== 1))) 
             perm_hit_2Tones_first(s,j,:) = [datasample(cell2mat(new_phase_analysis_2Tones_first.THETA_PHASE(new_phase_analysis_2Tones_first.HIT== 1 & new_phase_analysis_2Tones_first.CONDITION== 2 & new_phase_analysis_2Tones_first.TONE_POSITION== 1)),length(cell2mat(new_phase_analysis_2Tones_first.THETA_PHASE(new_phase_analysis_2Tones_first.HIT== 0 & new_phase_analysis_2Tones_first.CONDITION== 2 & new_phase_analysis_2Tones_first.TONE_POSITION== 1))),'Replace',false); NaN(300-length(cell2mat(new_phase_analysis_2Tones_first.THETA_PHASE((new_phase_analysis_2Tones_first.HIT== 0 & new_phase_analysis_2Tones_first.CONDITION== 2 & new_phase_analysis_2Tones_first.TONE_POSITION== 1)))),1)]';
             perm_miss_2Tones_first(s,j,:) = [cell2mat(new_phase_analysis_2Tones_first.THETA_PHASE(new_phase_analysis_2Tones_first.HIT== 0 & new_phase_analysis_2Tones_first.CONDITION== 2 & new_phase_analysis_2Tones_first.TONE_POSITION== 1)); NaN(300-length(cell2mat(new_phase_analysis_2Tones_first.THETA_PHASE(new_phase_analysis_2Tones_first.HIT== 0 & new_phase_analysis_2Tones_first.CONDITION== 2 & new_phase_analysis_2Tones_first.TONE_POSITION== 1))),1)]';  
        elseif length(cell2mat(new_phase_analysis_2Tones_first.THETA_PHASE(new_phase_analysis_2Tones_first.HIT== 0 & new_phase_analysis_2Tones_first.CONDITION== 2 & new_phase_analysis_2Tones_first.TONE_POSITION== 1))) > length(cell2mat(new_phase_analysis_2Tones_first.THETA_PHASE(new_phase_analysis_2Tones_first.HIT== 1 & new_phase_analysis_2Tones_first.CONDITION== 2 & new_phase_analysis_2Tones_first.TONE_POSITION== 1)))
             perm_hit_2Tones_first(s,j,:) = [datasample(cell2mat(new_phase_analysis_2Tones_first.THETA_PHASE(new_phase_analysis_2Tones_first.HIT== 1 & new_phase_analysis_2Tones_first.CONDITION== 2 & new_phase_analysis_2Tones_first.TONE_POSITION== 1)),length(cell2mat(new_phase_analysis_2Tones_first.THETA_PHASE(new_phase_analysis_2Tones_first.HIT== 1 & new_phase_analysis_2Tones_first.CONDITION== 2 & new_phase_analysis_2Tones_first.TONE_POSITION== 1))),'Replace',false); NaN(300-length(cell2mat(new_phase_analysis_2Tones_first.THETA_PHASE((new_phase_analysis_2Tones_first.HIT== 1 & new_phase_analysis_2Tones_first.CONDITION== 2 & new_phase_analysis_2Tones_first.TONE_POSITION== 1)))),1)]';
             perm_miss_2Tones_first(s,j,:) = [cell2mat(new_phase_analysis_2Tones_first.THETA_PHASE(new_phase_analysis_2Tones_first.HIT== 0 & new_phase_analysis_2Tones_first.CONDITION== 2 & new_phase_analysis_2Tones_first.TONE_POSITION== 1)); NaN(300-length(cell2mat(new_phase_analysis_2Tones_first.THETA_PHASE(new_phase_analysis_2Tones_first.HIT== 0 & new_phase_analysis_2Tones_first.CONDITION== 2 & new_phase_analysis_2Tones_first.TONE_POSITION== 1))),1)]';
        end
        
        %Generate the permutations in the 2Tones-T2 Condition;            
         if length(cell2mat(new_phase_analysis_2Tones_second.THETA_PHASE(new_phase_analysis_2Tones_second.HIT== 0 & new_phase_analysis_2Tones_second.CONDITION== 2 & new_phase_analysis_2Tones_second.TONE_POSITION== 2))) < length(cell2mat(new_phase_analysis_2Tones_second.THETA_PHASE(new_phase_analysis_2Tones_second.HIT== 1 & new_phase_analysis_2Tones_second.CONDITION== 2 & new_phase_analysis_2Tones_second.TONE_POSITION== 2))) 
             perm_hit_2Tones_second(s,j,:) = [datasample(cell2mat(new_phase_analysis_2Tones_second.THETA_PHASE(new_phase_analysis_2Tones_second.HIT== 1 & new_phase_analysis_2Tones_second.CONDITION== 2 & new_phase_analysis_2Tones_second.TONE_POSITION== 2)),length(cell2mat(new_phase_analysis_2Tones_second.THETA_PHASE(new_phase_analysis_2Tones_second.HIT== 0 & new_phase_analysis_2Tones_second.CONDITION== 2 & new_phase_analysis_2Tones_second.TONE_POSITION== 2))),'Replace',false); NaN(300-length(cell2mat(new_phase_analysis_2Tones_second.THETA_PHASE((new_phase_analysis_2Tones_second.HIT== 0 & new_phase_analysis_2Tones_second.CONDITION== 2 & new_phase_analysis_2Tones_second.TONE_POSITION== 2)))),1)]';
             perm_miss_2Tones_second(s,j,:) = [cell2mat(new_phase_analysis_2Tones_second.THETA_PHASE(new_phase_analysis_2Tones_second.HIT== 0 & new_phase_analysis_2Tones_second.CONDITION== 2 & new_phase_analysis_2Tones_second.TONE_POSITION== 2)); NaN(300-length(cell2mat(new_phase_analysis_2Tones_second.THETA_PHASE(new_phase_analysis_2Tones_second.HIT== 0 & new_phase_analysis_2Tones_second.CONDITION== 2 & new_phase_analysis_2Tones_second.TONE_POSITION== 2))),1)]';  
        elseif length(cell2mat(new_phase_analysis_2Tones_second.THETA_PHASE(new_phase_analysis_2Tones_second.HIT== 0 & new_phase_analysis_2Tones_second.CONDITION== 2 & new_phase_analysis_2Tones_second.TONE_POSITION== 2))) > length(cell2mat(new_phase_analysis_2Tones_second.THETA_PHASE(new_phase_analysis_2Tones_second.HIT== 1 & new_phase_analysis_2Tones_second.CONDITION== 2 & new_phase_analysis_2Tones_second.TONE_POSITION== 2)))
             perm_hit_2Tones_second(s,j,:) = [datasample(cell2mat(new_phase_analysis_2Tones_second.THETA_PHASE(new_phase_analysis_2Tones_second.HIT== 1 & new_phase_analysis_2Tones_second.CONDITION== 2 & new_phase_analysis_2Tones_second.TONE_POSITION== 2)),length(cell2mat(new_phase_analysis_2Tones_second.THETA_PHASE(new_phase_analysis_2Tones_second.HIT== 1 & new_phase_analysis_2Tones_second.CONDITION== 2 & new_phase_analysis_2Tones_second.TONE_POSITION== 2))),'Replace',false); NaN(300-length(cell2mat(new_phase_analysis_2Tones_second.THETA_PHASE((new_phase_analysis_2Tones_second.HIT== 1 & new_phase_analysis_2Tones_second.CONDITION== 2 & new_phase_analysis_2Tones_second.TONE_POSITION== 2)))),1)]';
             perm_miss_2Tones_second(s,j,:) = [cell2mat(new_phase_analysis_2Tones_second.THETA_PHASE(new_phase_analysis_2Tones_second.HIT== 0 & new_phase_analysis_2Tones_second.CONDITION== 2 & new_phase_analysis_2Tones_second.TONE_POSITION== 2)); NaN(300-length(cell2mat(new_phase_analysis_2Tones_second.THETA_PHASE(new_phase_analysis_2Tones_second.HIT== 0 & new_phase_analysis_2Tones_second.CONDITION== 2 & new_phase_analysis_2Tones_second.TONE_POSITION== 2))),1)]';
        end      
    end

    clear new_phase_analysis_1Tone new_phase_analysis_2Tones_first new_phase_analysis_2Tones_second
    
end

% cd D:\Birmingham_2018\BEHAVIOR_project\ANALYSES_DATA\;
% save permutations_10000_HitvsMiss perm_hit_1Tone perm_miss_1Tone perm_hit_2Tones_first perm_miss_2Tones_first perm_hit_2Tones_second perm_miss_2Tones_second

%% Now perform test statistically whether mean hit phase is greater than misses for the three conditions;
clearvars;clc;
addpath C:\toolbox\CircHist-master; addpath C:\toolbox\Circular_Statistics_Toolbox\; addpath D:\Birmingham_2018\BEHAVIOR_project\ANALYSES_DATA\;
load permutations_10000_HitvsMiss;


%subjects ID;
subjects = 1:29;

%Make sure it is the same number used to generate the permutations before; 
permutations = 10000;

%condition to test;
condition  = 22;

%Calculate the mean phase of each of the 10000 permutations in every participants;
clear mean_hits mean_miss

for s = subjects
try
    
    for j = 1:size(perm_hit_1Tone(:,:,:),2)
 
        switch condition
            
            case 1
                mean_hits(s,j) = circ_mean(permute(perm_hit_1Tone(s,j,~isnan(perm_hit_1Tone(s,j,:))),[3 1 2]));
                mean_miss(s,j) = circ_mean(permute(perm_miss_1Tone(s,j,~isnan(perm_miss_1Tone(s,j,:))),[3 1 2]));
    
            case 21
                mean_hits(s,j) = circ_mean(permute(perm_hit_2Tones_first(s,j,~isnan(perm_hit_2Tones_first(s,j,:))),[3 1 2]));
                mean_miss(s,j) = circ_mean(permute(perm_miss_2Tones_first(s,j,~isnan(perm_miss_2Tones_first(s,j,:))),[3 1 2]));

            case 22        
                mean_hits(s,j) = circ_mean(permute(perm_hit_2Tones_second(s,j,~isnan(perm_hit_2Tones_second(s,j,:))),[3 1 2]));
                mean_miss(s,j) = circ_mean(permute(perm_miss_2Tones_second(s,j,~isnan(perm_miss_2Tones_second(s,j,:))),[3 1 2]));
        end    
          
    end
    
catch
end
    
end

%Remove the bad subjects to exclude from statistic analyses;
bad_subjects = [1 4 25 26 27];
mean_hits(bad_subjects,:) = [];
mean_miss(bad_subjects,:) = [];

%now, for each mean phase permuted, perform a Rayleigh's test to get the hit/miss vector lengths(z_value);
clear Z_hit Z_miss   

for t = 1:length(mean_hits(1,:)) 
    
    [p_mean_hit, z_hit] = circ_rtest(mean_hits(:,t));
    [p_mean_miss, z_miss] = circ_rtest(mean_miss(:,t));
    
    Z_hit(t,1)= p_mean_hit;
    Z_miss(t,1)= p_mean_miss;
    Z_hit(t,2)= z_hit;
    Z_miss(t,2)= z_miss;
     
    clear p_mean_hit p_mean_miss z_hit z_miss
  
end  
    
%Now calculate the mean phase on the original data set to get real Zvalue difference between hits and miss;
clear ori_hit ori_miss ori_mean_hit ori_mean_miss

for s = subjects
    
    cd(['C:\ANALYSES_DATA\subj' num2str(s)]);
    load(['subj' num2str(s) '-Phase_analyses']);

    switch condition
    
        case 1
            ori_hit(s,:) = [cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.CONDITION==1)); NaN(100 - length(cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.CONDITION==1))),1)]';
            ori_miss(s,:) = [cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==0 & PHASE_ANALYSES.CONDITION==1)); NaN(100 - length(cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==0 & PHASE_ANALYSES.CONDITION==1))),1)]';
            ori_mean_hit(s,1) = circ_mean(ori_hit(s,~isnan(ori_hit(s,:)))');
            ori_mean_miss(s,1) = circ_mean(ori_miss(s,~isnan(ori_miss(s,:)))');
            
        case 21    
            ori_hit(s,:) = [cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==1)); NaN(100 - length(cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==1))),1)]';
            ori_miss(s,:) = [cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==0 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==1)); NaN(100 - length(cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==0 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==1))),1)]';
            ori_mean_hit(s,1) = circ_mean(ori_hit(s,~isnan(ori_hit(s,:)))');
            ori_mean_miss(s,1) = circ_mean(ori_miss(s,~isnan(ori_miss(s,:)))');
            
        case 22
            ori_hit(s,:) = [cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==2)); NaN(100 - length(cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==1 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==2))),1)]';
            ori_miss(s,:) = [cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==0 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==2)); NaN(100 - length(cell2mat(PHASE_ANALYSES.THETA_PHASE(PHASE_ANALYSES.HIT==0 & PHASE_ANALYSES.CONDITION==2 & PHASE_ANALYSES.TONE_POSITION==2))),1)]';
            ori_mean_hit(s,1) = circ_mean(ori_hit(s,~isnan(ori_hit(s,:)))');
            ori_mean_miss(s,1) = circ_mean(ori_miss(s,~isnan(ori_miss(s,:)))');
            
    end
    
    clear PHASE_ANALYSES alltrials participant_parameters

end
   
%Remove the bad subjects to exclude from statistic analyses;
bad_subjects = [1 4 25 26 27];
ori_mean_hit(bad_subjects,:) = [];
ori_mean_miss(bad_subjects,:) = [];   

[~, ori_z_hit] = circ_rtest(ori_mean_hit(:,1));
[~, ori_z_miss] = circ_rtest(ori_mean_miss(:,1));

%Calculate the difference of vector lengths (z_values) hit-miss for each permutation and for the original data;
%Then sort in descend order;
diff_Z_HitvsMiss(:,1) = Z_hit(:,2) - Z_miss(:,2);
diff_Z_HitvsMiss(:,1)= sortrows(diff_Z_HitvsMiss(:,1),'descend');
diff_Z_orig = ori_z_hit - ori_z_miss;
P_value_perm = sum(diff_Z_HitvsMiss(:,1)> diff_Z_orig)/(length(diff_Z_HitvsMiss(:,1)));

display(['P_value_perm: ' num2str(P_value_perm)]);
display(['Effect size: ' num2str(circ_r(ori_mean_hit(:,1))- circ_r(ori_mean_miss(:,1)))]);    

%%