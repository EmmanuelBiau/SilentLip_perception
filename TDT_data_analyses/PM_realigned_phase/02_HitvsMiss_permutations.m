%% First, we have to split hit/miss trials and realign the visual phase on preferred bin for each participant;
clear;clc; addpath C:\toolbox\CircHist-master; addpath C:\toolbox\Circular_Statistics_Toolbox\;  
cd C:\ANALYSES_DATA\analyses_realign_phase\;
load ONE_TONE_ANALYSES;load TWO_TONES_ANALYSES;

% Determine the number of bins in theta phase;
bin = pi/10;
theta_phase = -pi:bin:pi-bin;
centerbin = length(theta_phase)/2+1;

% number of permutations;
permutations = 10000;

% subjects IDs;
subjects = 1:29;

% I- Compile the hits from real data of all participants together;
ori_OneTone_first_hit = []; ori_OneTone_second_hit = []; ori_OneTone_first_miss = []; ori_OneTone_second_miss = [];
ori_TwoTones_first_hit = []; ori_TwoTones_second_hit = []; ori_TwoTones_first_miss = []; ori_TwoTones_second_miss = [];

for s = subjects
    
    OneTone_tmp = ONE_TONE_ANALYSES(ismember(ONE_TONE_ANALYSES.PARTICIPANT,s),:); TwoTones_tmp = TWO_TONES_ANALYSES(ismember(TWO_TONES_ANALYSES.PARTICIPANT,s),:);     
    ori_OneTone_first_hit(s,:) = [OneTone_tmp.OneTone_phase(OneTone_tmp.Hit == 1 & OneTone_tmp.OneTone_position == 1); NaN(100 - length(OneTone_tmp.OneTone_phase(OneTone_tmp.Hit == 1 & OneTone_tmp.OneTone_position == 1)),1)]';
    ori_OneTone_second_hit(s,:) = [OneTone_tmp.OneTone_phase(OneTone_tmp.Hit == 1 & OneTone_tmp.OneTone_position == 2); NaN(100 - length(OneTone_tmp.OneTone_phase(OneTone_tmp.Hit == 1 & OneTone_tmp.OneTone_position == 2)),1)]';  
    ori_OneTone_first_miss(s,:) = [OneTone_tmp.OneTone_phase(OneTone_tmp.Hit == 0 & OneTone_tmp.OneTone_position == 1); NaN(100 - length(OneTone_tmp.OneTone_phase(OneTone_tmp.Hit == 0 & OneTone_tmp.OneTone_position == 1)),1)]';
    ori_OneTone_second_miss(s,:) = [OneTone_tmp.OneTone_phase(OneTone_tmp.Hit == 0 & OneTone_tmp.OneTone_position == 2); NaN(100 - length(OneTone_tmp.OneTone_phase(OneTone_tmp.Hit == 0 & OneTone_tmp.OneTone_position == 2)),1)]';  
    ori_TwoTones_first_hit(s,:) = [TwoTones_tmp.Tone_phase(TwoTones_tmp.Hit == 1 & TwoTones_tmp.Tone_position == 1); NaN(100 - length(TwoTones_tmp.Tone_phase(TwoTones_tmp.Hit == 1 & TwoTones_tmp.Tone_position == 1)),1)]';
    ori_TwoTones_second_hit(s,:) = [TwoTones_tmp.Tone_phase(TwoTones_tmp.Hit == 1 & TwoTones_tmp.Tone_position == 2); NaN(100 - length(TwoTones_tmp.Tone_phase(TwoTones_tmp.Hit == 1 & TwoTones_tmp.Tone_position == 2)),1)]';
    ori_TwoTones_first_miss(s,:) = [TwoTones_tmp.Tone_phase(TwoTones_tmp.Hit == 0 & TwoTones_tmp.Tone_position == 1); NaN(100 - length(TwoTones_tmp.Tone_phase(TwoTones_tmp.Hit == 0 & TwoTones_tmp.Tone_position == 1)),1)]';
    ori_TwoTones_second_miss(s,:) = [TwoTones_tmp.Tone_phase(TwoTones_tmp.Hit == 0 & TwoTones_tmp.Tone_position == 2); NaN(100 - length(TwoTones_tmp.Tone_phase(TwoTones_tmp.Hit == 0 & TwoTones_tmp.Tone_position == 2)),1)]';

    % II- Calculate the number of p(hit) for each phase bin;   
    if nansum(ori_OneTone_first_hit(s,:)) ~= 0       
      ori_bin_OneTone_first_hit(s,:) = histcounts(ori_OneTone_first_hit(s,:),-pi:bin:pi);
    else
      ori_bin_OneTone_first_hit(s,:) = zeros(1,length(theta_phase)); 
    end

    if nansum(ori_OneTone_second_hit(s,:)) ~= 0       
      ori_bin_OneTone_second_hit(s,:) = histcounts(ori_OneTone_second_hit(s,:),-pi:bin:pi);
    else
      ori_bin_OneTone_second_hit(s,:) = zeros(1,length(theta_phase)); 
    end      

    if nansum(ori_TwoTones_first_hit(s,:)) ~= 0       
      ori_bin_TwoTones_first_hit(s,:) = histcounts(ori_TwoTones_first_hit(s,:),-pi:bin:pi);
    else
      ori_bin_TwoTones_first_hit(s,:) = zeros(1,length(theta_phase)); 
    end

    if nansum(ori_TwoTones_second_hit(s,:)) ~= 0       
      ori_bin_TwoTones_second_hit(s,:) = histcounts(ori_TwoTones_second_hit(s,:),-pi:bin:pi);
    else
      ori_bin_TwoTones_second_hit(s,:) = zeros(1,length(theta_phase)); 
    end            

    if nansum(ori_OneTone_first_miss(s,:)) ~= 0       
      ori_bin_OneTone_first_miss(s,:) = histcounts(ori_OneTone_first_miss(s,:),-pi:bin:pi);
    else
      ori_bin_OneTone_first_miss(s,:) = zeros(1,length(theta_phase)); 
    end

    if nansum(ori_OneTone_second_miss(s,:)) ~= 0       
      ori_bin_OneTone_second_miss(s,:) = histcounts(ori_OneTone_second_miss(s,:),-pi:bin:pi);
    else
      ori_bin_OneTone_second_miss(s,:) = zeros(1,length(theta_phase)); 
    end      

    if nansum(ori_TwoTones_first_miss(s,:)) ~= 0       
      ori_bin_TwoTones_first_miss(s,:) = histcounts(ori_TwoTones_first_miss(s,:),-pi:bin:pi);
    else
      ori_bin_TwoTones_first_miss(s,:) = zeros(1,length(theta_phase)); 
    end

    if nansum(ori_TwoTones_second_miss(s,:)) ~= 0       
      ori_bin_TwoTones_second_miss(s,:) = histcounts(ori_TwoTones_second_miss(s,:),-pi:bin:pi);
    else
      ori_bin_TwoTones_second_miss(s,:) = zeros(1,length(theta_phase)); 
    end 
   
    % III- Now realign the phase at the maximum hit phase bin for each participants and conditions;
    [~,max_phase_hit] = max(ori_bin_OneTone_first_hit(s,:));
    max_toshift_hit = length(theta_phase)/2+1 - max_phase_hit;
    ori_maxPhase_realign.OneTone_first_hit(s,:) = circshift(ori_bin_OneTone_first_hit(s,:),max_toshift_hit);
    [~,max_phase_hit] = max(ori_bin_OneTone_second_hit(s,:));
    max_toshift_hit = length(theta_phase)/2+1 - max_phase_hit;
    ori_maxPhase_realign.OneTone_second_hit(s,:) = circshift(ori_bin_OneTone_second_hit(s,:),max_toshift_hit);
    [~,max_phase_hit] = max(ori_bin_TwoTones_first_hit(s,:));
    max_toshift_hit = length(theta_phase)/2+1 - max_phase_hit;
    ori_maxPhase_realign.TwoTones_first_hit(s,:) = circshift(ori_bin_TwoTones_first_hit(s,:),max_toshift_hit);
    [~,max_phase_hit] = max(ori_bin_TwoTones_second_hit(s,:));
    max_toshift_hit = length(theta_phase)/2+1 - max_phase_hit;
    ori_maxPhase_realign.TwoTones_second_hit(s,:) = circshift(ori_bin_TwoTones_second_hit(s,:),max_toshift_hit);
    [~,max_phase_miss] = max(ori_bin_OneTone_first_miss(s,:));
    max_toshift_miss = length(theta_phase)/2+1 - max_phase_miss;
    ori_maxPhase_realign.OneTone_first_miss(s,:) = circshift(ori_bin_OneTone_first_miss(s,:),max_toshift_miss);
    [~,max_phase_miss] = max(ori_bin_OneTone_second_miss(s,:));
    max_toshift_miss = length(theta_phase)/2+1 - max_phase_miss;
    ori_maxPhase_realign.OneTone_second_miss(s,:) = circshift(ori_bin_OneTone_second_miss(s,:),max_toshift_miss);
    [~,max_phase_miss] = max(ori_bin_TwoTones_first_miss(s,:));
    max_toshift_miss = length(theta_phase)/2+1 - max_phase_miss;
    ori_maxPhase_realign.TwoTones_first_miss(s,:) = circshift(ori_bin_TwoTones_first_miss(s,:),max_toshift_miss);
    [~,max_phase_miss] = max(ori_bin_TwoTones_second_miss(s,:));
    max_toshift_miss = length(theta_phase)/2+1 - max_phase_miss;
    ori_maxPhase_realign.TwoTones_second_miss(s,:) = circshift(ori_bin_TwoTones_second_miss(s,:),max_toshift_miss);

    clear max_phase_hit max_phase_miss max_toshift_hit max_toshift_miss ...
          max_phase_miss max_phase_miss max_toshift_miss max_toshift_miss 


    % compute phase bin entrainment hit/miss and size effect at first/second tones in OneTone and Two Tones conditions;
    centerbin = length(theta_phase)/2+1;
    max_phs_OneTone_first_hit = ori_maxPhase_realign.OneTone_first_hit(s,:);
    max_phs_OneTone_second_hit = ori_maxPhase_realign.OneTone_second_hit(s,:);
    max_phs_TwoTones_first_hit = ori_maxPhase_realign.TwoTones_first_hit(s,:);
    max_phs_TwoTones_second_hit = ori_maxPhase_realign.TwoTones_second_hit(s,:);
    max_phs_OneTone_first_miss = ori_maxPhase_realign.OneTone_first_miss(s,:);
    max_phs_OneTone_second_miss = ori_maxPhase_realign.OneTone_second_miss(s,:);
    max_phs_TwoTones_first_miss = ori_maxPhase_realign.TwoTones_first_miss(s,:);
    max_phs_TwoTones_second_miss = ori_maxPhase_realign.TwoTones_second_miss(s,:);  

    % hit-miss: OneTone first;
    tocomp1 = ((max_phs_OneTone_first_hit(:,centerbin-1)+ max_phs_OneTone_first_hit(:,centerbin+1))./2) - ((max_phs_OneTone_first_hit(:,2) + max_phs_OneTone_first_hit(:,end))./2);        
    tocomp2 = ((max_phs_OneTone_first_miss(:,centerbin-1)+ max_phs_OneTone_first_miss(:,centerbin+1))./2) - ((max_phs_OneTone_first_miss(:,2) + max_phs_OneTone_first_miss(:,end))./2);
    ori_diff_phs_OneTone_first(s,1) = tocomp1 - tocomp2;
    % hit-miss: OneTone second;
    tocomp3 = ((max_phs_OneTone_second_hit(:,centerbin-1)+ max_phs_OneTone_second_hit(:,centerbin+1))./2) - ((max_phs_OneTone_second_hit(:,2) + max_phs_OneTone_second_hit(:,end))./2);        
    tocomp4 = ((max_phs_OneTone_second_miss(:,centerbin-1)+ max_phs_OneTone_second_miss(:,centerbin+1))./2) - ((max_phs_OneTone_second_miss(:,2) + max_phs_OneTone_second_miss(:,end))./2);
    ori_diff_phs_OneTone_second(s,1) = tocomp3 - tocomp4;
    % hit-miss: TwoTones first;
    tocomp5 = ((max_phs_TwoTones_first_hit(:,centerbin-1)+ max_phs_TwoTones_first_hit(:,centerbin+1))./2) - ((max_phs_TwoTones_first_hit(:,2) + max_phs_TwoTones_first_hit(:,end))./2);        
    tocomp6 = ((max_phs_TwoTones_first_miss(:,centerbin-1)+ max_phs_TwoTones_first_miss(:,centerbin+1))./2) - ((max_phs_TwoTones_first_miss(:,2) + max_phs_TwoTones_first_miss(:,end))./2);        
    ori_diff_phs_TwoTones_first(s,1) = tocomp5 - tocomp6;    
    % hit-miss: TwoTones second;
    tocomp7 = ((max_phs_TwoTones_second_hit(:,centerbin-1)+ max_phs_TwoTones_second_hit(:,centerbin+1))./2) - ((max_phs_TwoTones_second_hit(:,2) + max_phs_TwoTones_second_hit(:,end))./2);        
    tocomp8 = ((max_phs_TwoTones_second_miss(:,centerbin-1)+ max_phs_TwoTones_second_miss(:,centerbin+1))./2) - ((max_phs_TwoTones_second_miss(:,2) + max_phs_TwoTones_second_miss(:,end))./2);
    ori_diff_phs_TwoTones_second(s,1) = tocomp7 - tocomp8;
    % effect size;
    ori_effect_size_TwoTones(s,1) = tocomp7 - tocomp8 - tocomp5 - tocomp6; 
    ori_effect_size_OneTone(s,1) = tocomp3 - tocomp4 - tocomp1 - tocomp2; 

    clear tocomp1 tocomp2 tocomp3 tocomp4 tocomp5 tocomp6 tocomp7 tocomp8 

end

% Remove the bad subjects to exclude from statistic analyses;
bad_subjects = [1 4 25 26 27];
ori_diff_phs_OneTone_first(bad_subjects,:) = []; ori_diff_phs_OneTone_second(bad_subjects,:) = [];
ori_diff_phs_TwoTones_first(bad_subjects,:) = []; ori_diff_phs_TwoTones_second(bad_subjects,:) = [];
ori_effect_size_TwoTones(bad_subjects,:) = []; ori_effect_size_OneTone(bad_subjects,:) = [];
ori_diff_phs_OneTone_first(ori_diff_phs_OneTone_first==0,:) = 0.00001; ori_diff_phs_OneTone_second(ori_diff_phs_OneTone_second==0,:) = 0.00001;
ori_diff_phs_TwoTones_first(ori_diff_phs_TwoTones_first==0,:) = 0.00001; ori_diff_phs_TwoTones_second(ori_diff_phs_TwoTones_second==0,:) = 0.00001;
ori_effect_size_TwoTones(ori_effect_size_TwoTones==0,:) = 0.00001; ori_effect_size_OneTone(ori_effect_size_OneTone==0,:) = 0.00001;

% run the ttest2 on original data to get the t-value;
% TwoTones condition;
[~,ori_TwoTones_first_p_max_adj(1,2),~,stats1] = ttest(ori_diff_phs_TwoTones_first(:,1),0,'tail','right');
ori_TwoTones_first_p_max_adj(1,1) = stats1.tstat;
clear stats1
[~,ori_TwoTones_second_p_max_adj(1,2),~,stats2] = ttest(ori_diff_phs_TwoTones_second(:,1),0,'tail','right');
ori_TwoTones_second_p_max_adj(1,1) = stats2.tstat;
clear stats2
[~,ori_TwoTones_p_max_adj(1,2),~,stats3] = ttest2(ori_diff_phs_TwoTones_second(:,1),ori_diff_phs_TwoTones_first(:,1),'tail','right');
ori_TwoTones_p_max_adj(1,1) = stats3.tstat;
clear stat3
% OneTone condition;
[~,ori_OneTone_first_p_max_adj(1,2),~,stats4] = ttest(ori_diff_phs_OneTone_first(:,1),0,'tail','right');
ori_OneTone_first_p_max_adj(1,1) = stats4.tstat;
clear stats4
[~,ori_OneTone_second_p_max_adj(1,2),~,stats5] = ttest(ori_diff_phs_OneTone_second(:,1),0,'tail','right');
ori_OneTone_second_p_max_adj(1,1) = stats5.tstat;
clear stats5 
[~,ori_OneTone_p_max_adj(1,2),~,stats6] = ttest2(ori_diff_phs_OneTone_second(:,1),ori_diff_phs_OneTone_first(:,1),'tail','right');
ori_OneTone_p_max_adj(1,1) = stats6.tstat;
clear stat6
[~,ori_EffectSize_p_max_adj(1,2),~,stats7] = ttest2(ori_effect_size_TwoTones(:,1),ori_effect_size_OneTone(:,1),'tail','right');
ori_EffectSize_p_max_adj(1,1) = stats7.tstat;
clear stat7

% store T-values as ori_Effect Sizes;
Ori_effect_size.ES_OneTone_first = ori_OneTone_first_p_max_adj(1,1);
Ori_effect_size.ES_OneTone_second = ori_OneTone_second_p_max_adj(1,1);
Ori_effect_size.ES_TwoTones_first = ori_TwoTones_first_p_max_adj(1,1);
Ori_effect_size.ES_TwoTones_second = ori_TwoTones_second_p_max_adj(1,1);
Ori_effect_size.ES_TwoTones = ori_TwoTones_p_max_adj(1,1);
Ori_effect_size.ES_OneTone = ori_OneTone_p_max_adj(1,1);
Ori_effect_size.ES_T1vsT2 = ori_EffectSize_p_max_adj(1,1);

ori_EffectSize_OneTone_second = mean(ori_diff_phs_OneTone_second);
ori_EffectSize_OneTone_first = mean(ori_diff_phs_OneTone_first);
ori_EffectSize_Two_tones_second = mean(ori_diff_phs_TwoTones_second);
ori_EffectSize_Two_tones_first = mean(ori_diff_phs_TwoTones_first);
ori_EffectSize_TwoTones = mean(ori_diff_phs_TwoTones_second - ori_diff_phs_TwoTones_first);
ori_EffectSize_OneTone = mean(ori_diff_phs_OneTone_second - ori_diff_phs_OneTone_first);

%% Compute the permutations;

% I- Compile the hit/miss phase angle of all participants together;
perm_OneTone_first_miss = []; perm_OneTone_first_hit = []; perm_OneTone_second_miss = []; perm_OneTone_second_hit = [];
perm_TwoTones_first_miss = []; perm_TwoTones_first_hit = []; perm_TwoTones_second_miss = []; perm_TwoTones_second_hit = [];

f = waitbar(0,'Running permutations now','Name','Work in Progress');

for s = subjects
    
    OneTone_tmp = ONE_TONE_ANALYSES(ismember(ONE_TONE_ANALYSES.PARTICIPANT,s),:); 
    TwoTones_tmp = TWO_TONES_ANALYSES(ismember(TWO_TONES_ANALYSES.PARTICIPANT,s),:); 
    
    new_conlab1 = randperm(size(OneTone_tmp,1));  
    temp_hit1 = OneTone_tmp.OneTone_position;
    OneTone_tmp.OneTone_position = temp_hit1(new_conlab1);
    
    new_conlab2 = randperm(size(TwoTones_tmp,1));  
    temp_hit2 = TwoTones_tmp.Tone_position;
    TwoTones_tmp.Tone_position = temp_hit2(new_conlab2);
   
    for j = 1:permutations
    
        % OneTone first hit/miss permutations;
        if length(OneTone_tmp.OneTone_phase(OneTone_tmp.OneTone_position == 1 & OneTone_tmp.Hit== 0)) < length(OneTone_tmp.OneTone_phase(OneTone_tmp.OneTone_position == 1 & OneTone_tmp.Hit== 1))
            perm_OneTone_first_miss(s,j,:) = [OneTone_tmp.OneTone_phase(OneTone_tmp.OneTone_position == 1 & OneTone_tmp.Hit== 0); NaN(100 - length(OneTone_tmp.OneTone_phase(OneTone_tmp.OneTone_position == 1 & OneTone_tmp.Hit== 0)),1)]';
            perm_OneTone_first_hit(s,j,:) = [datasample(OneTone_tmp.OneTone_phase(OneTone_tmp.OneTone_position == 1 & OneTone_tmp.Hit==1),length(OneTone_tmp.OneTone_phase(OneTone_tmp.OneTone_position == 1 & OneTone_tmp.Hit==0)),'Replace',false); NaN(100 - length(OneTone_tmp.OneTone_phase(OneTone_tmp.OneTone_position == 1 & OneTone_tmp.Hit== 0)),1)]';   
        elseif length(OneTone_tmp.OneTone_phase(OneTone_tmp.OneTone_position == 1 & OneTone_tmp.Hit== 0)) > length(OneTone_tmp.OneTone_phase(OneTone_tmp.OneTone_position == 1 & OneTone_tmp.Hit== 1))
            perm_OneTone_first_hit(s,j,:) = [OneTone_tmp.OneTone_phase(OneTone_tmp.OneTone_position == 1 & OneTone_tmp.Hit== 1); NaN(100 - length(OneTone_tmp.OneTone_phase(OneTone_tmp.OneTone_position == 1 & OneTone_tmp.Hit== 1)),1)]';
            perm_OneTone_first_miss(s,j,:) = [datasample(OneTone_tmp.OneTone_phase(OneTone_tmp.OneTone_position == 1 & OneTone_tmp.Hit==0),length(OneTone_tmp.OneTone_phase(OneTone_tmp.OneTone_position == 1 & OneTone_tmp.Hit==1)),'Replace',false); NaN(100 - length(OneTone_tmp.OneTone_phase(OneTone_tmp.OneTone_position == 1 & OneTone_tmp.Hit== 1)),1)]';   
        end     
        % OneTone second hit/miss permutations;       
        if length(OneTone_tmp.OneTone_phase(OneTone_tmp.OneTone_position == 2 & OneTone_tmp.Hit== 0)) < length(OneTone_tmp.OneTone_phase(OneTone_tmp.OneTone_position == 2 & OneTone_tmp.Hit== 1))
            perm_OneTone_second_miss(s,j,:) = [OneTone_tmp.OneTone_phase(OneTone_tmp.OneTone_position == 2 & OneTone_tmp.Hit== 0); NaN(100 - length(OneTone_tmp.OneTone_phase(OneTone_tmp.OneTone_position == 2 & OneTone_tmp.Hit== 0)),1)]';
            perm_OneTone_second_hit(s,j,:) = [datasample(OneTone_tmp.OneTone_phase(OneTone_tmp.OneTone_position == 2 & OneTone_tmp.Hit==1),length(OneTone_tmp.OneTone_phase(OneTone_tmp.OneTone_position == 2 & OneTone_tmp.Hit==0)),'Replace',false); NaN(100 - length(OneTone_tmp.OneTone_phase(OneTone_tmp.OneTone_position == 2 & OneTone_tmp.Hit== 0)),1)]';   
        elseif length(OneTone_tmp.OneTone_phase(OneTone_tmp.OneTone_position == 2 & OneTone_tmp.Hit== 0)) > length(OneTone_tmp.OneTone_phase(OneTone_tmp.OneTone_position == 2 & OneTone_tmp.Hit== 1))
            perm_OneTone_second_hit(s,j,:) = [OneTone_tmp.OneTone_phase(OneTone_tmp.OneTone_position == 2 & OneTone_tmp.Hit== 1); NaN(100 - length(OneTone_tmp.OneTone_phase(OneTone_tmp.OneTone_position == 2 & OneTone_tmp.Hit== 1)),1)]';
            perm_OneTone_second_miss(s,j,:) = [datasample(OneTone_tmp.OneTone_phase(OneTone_tmp.OneTone_position == 2 & OneTone_tmp.Hit==0),length(OneTone_tmp.OneTone_phase(OneTone_tmp.OneTone_position == 2 & OneTone_tmp.Hit==1)),'Replace',false); NaN(100 - length(OneTone_tmp.OneTone_phase(OneTone_tmp.OneTone_position == 2 & OneTone_tmp.Hit== 1)),1)]';   
        end      
        % TwoTones first hit/miss permutations;
        if length(TwoTones_tmp.Tone_phase(TwoTones_tmp.Tone_position == 1 & TwoTones_tmp.Hit== 0)) < length(TwoTones_tmp.Tone_phase(TwoTones_tmp.Tone_position == 1 & TwoTones_tmp.Hit== 1))
            perm_TwoTones_first_miss(s,j,:) = [TwoTones_tmp.Tone_phase(TwoTones_tmp.Tone_position == 1 & TwoTones_tmp.Hit== 0); NaN(100 - length(TwoTones_tmp.Tone_phase(TwoTones_tmp.Tone_position == 1 & TwoTones_tmp.Hit== 0)),1)]';
            perm_TwoTones_first_hit(s,j,:) = [datasample(TwoTones_tmp.Tone_phase(TwoTones_tmp.Tone_position == 1 & TwoTones_tmp.Hit==1),length(TwoTones_tmp.Tone_phase(TwoTones_tmp.Tone_position == 1 & TwoTones_tmp.Hit==0)),'Replace',false); NaN(100 - length(TwoTones_tmp.Tone_phase(TwoTones_tmp.Tone_position == 1 & TwoTones_tmp.Hit== 0)),1)]';   
        elseif length(TwoTones_tmp.Tone_phase(TwoTones_tmp.Tone_position == 1 & TwoTones_tmp.Hit== 0)) > length(TwoTones_tmp.Tone_phase(TwoTones_tmp.Tone_position == 1 & TwoTones_tmp.Hit== 1))
            perm_TwoTones_first_hit(s,j,:) = [TwoTones_tmp.Tone_phase(TwoTones_tmp.Tone_position == 1 & TwoTones_tmp.Hit== 1); NaN(100 - length(TwoTones_tmp.Tone_phase(TwoTones_tmp.Tone_position == 1 & TwoTones_tmp.Hit== 1)),1)]';
            perm_TwoTones_first_miss(s,j,:) = [datasample(TwoTones_tmp.Tone_phase(TwoTones_tmp.Tone_position == 1 & TwoTones_tmp.Hit==0),length(TwoTones_tmp.Tone_phase(TwoTones_tmp.Tone_position == 1 & TwoTones_tmp.Hit==1)),'Replace',false); NaN(100 - length(TwoTones_tmp.Tone_phase(TwoTones_tmp.Tone_position == 1 & TwoTones_tmp.Hit== 1)),1)]';   
        end                
        % TwoTones second hit/miss permutations;
        if length(TwoTones_tmp.Tone_phase(TwoTones_tmp.Tone_position == 2 & TwoTones_tmp.Hit== 0)) < length(TwoTones_tmp.Tone_phase(TwoTones_tmp.Tone_position == 2 & TwoTones_tmp.Hit== 1))
            perm_TwoTones_second_miss(s,j,:) = [TwoTones_tmp.Tone_phase(TwoTones_tmp.Tone_position == 2 & TwoTones_tmp.Hit== 0); NaN(100 - length(TwoTones_tmp.Tone_phase(TwoTones_tmp.Tone_position == 2 & TwoTones_tmp.Hit== 0)),1)]';
            perm_TwoTones_second_hit(s,j,:) = [datasample(TwoTones_tmp.Tone_phase(TwoTones_tmp.Tone_position == 2 & TwoTones_tmp.Hit==1),length(TwoTones_tmp.Tone_phase(TwoTones_tmp.Tone_position == 2 & TwoTones_tmp.Hit==0)),'Replace',false); NaN(100 - length(TwoTones_tmp.Tone_phase(TwoTones_tmp.Tone_position == 2 & TwoTones_tmp.Hit== 0)),1)]';   
        elseif length(TwoTones_tmp.Tone_phase(TwoTones_tmp.Tone_position == 2 & TwoTones_tmp.Hit== 0)) > length(TwoTones_tmp.Tone_phase(TwoTones_tmp.Tone_position == 2 & TwoTones_tmp.Hit== 1))
            perm_TwoTones_second_hit(s,j,:) = [TwoTones_tmp.Tone_phase(TwoTones_tmp.Tone_position == 2 & TwoTones_tmp.Hit== 1); NaN(100 - length(TwoTones_tmp.Tone_phase(TwoTones_tmp.Tone_position == 2 & TwoTones_tmp.Hit== 1)),1)]';
            perm_TwoTones_second_miss(s,j,:) = [datasample(TwoTones_tmp.Tone_phase(TwoTones_tmp.Tone_position == 2 & TwoTones_tmp.Hit==0),length(TwoTones_tmp.Tone_phase(TwoTones_tmp.Tone_position == 2 & TwoTones_tmp.Hit==1)),'Replace',false); NaN(100 - length(TwoTones_tmp.Tone_phase(TwoTones_tmp.Tone_position == 2 & TwoTones_tmp.Hit== 1)),1)]';   
        end 
    end

    % Update waitbar and message;
    waitbar(s/length(subjects),f,sprintf('Running permutations now'));

end

delete(f);
disp('ALL PERMUTATIONS GENERATED')

% %save permutations\;
% cd xxx;
% save (['permutations_10000_RealignPhs_HitvsMiss_' num2str(pi/bin) 'bins'], 'perm_OneTone_first_miss','perm_OneTone_first_hit','perm_OneTone_second_miss','perm_OneTone_second_hit', ...
%                                                         'perm_TwoTones_first_miss','perm_TwoTones_first_hit','perm_TwoTones_second_miss','perm_TwoTones_second_hit', ...
%                                                         'bin','theta_phase','centerbin','permutations','subjects','bad_subjects', ...
%                                                         'ori_TwoTones_first_p_max_adj', 'ori_TwoTones_second_p_max_adj', 'ori_TwoTones_p_max_adj', 'ori_OneTone_first_p_max_adj', 'ori_OneTone_second_p_max_adj',...
%                                                         'ori_OneTone_p_max_adj', 'ori_EffectSize_p_max_adj');

% II- Calculate the number of p(hit) for each phase bin;   
clc; % load permutations_10000_RealignPhs_HitvsMiss_10bins;

% Realign the visual phase on preferred bin for hit/miss permutations as in original data;
f = waitbar(0,'Realigning the phase of permutations','Name','Work in Progress');
perm_bin_OneTone_first_hit = []; perm_bin_OneTone_second_hit = [];
perm_bin_OneTone_first_miss = []; perm_bin_OneTone_second_miss = [];
perm_bin_TwoTones_first_hit = []; perm_bin_TwoTones_second_hit = []; 
perm_bin_TwoTones_first_miss = []; perm_bin_TwoTones_second_miss = []; 
perm_maxPhase_realign = []; perm_minPhase_realign = [];

for s = subjects
try
    
    for j = 1:size(perm_OneTone_first_hit(:,:,:),2)        
       
        if nansum(perm_OneTone_first_hit(s,j,:)) ~= 0       
          perm_bin_OneTone_first_hit(s,j,:) = histcounts(perm_OneTone_first_hit(s,j,:),-pi:bin:pi);
        else
          perm_bin_OneTone_first_hit(s,j,:) = zeros(1,length(theta_phase)); 
        end

        if nansum(perm_OneTone_second_hit(s,j,:)) ~= 0       
          perm_bin_OneTone_second_hit(s,j,:) = histcounts(perm_OneTone_second_hit(s,j,:),-pi:bin:pi);
        else
          perm_bin_OneTone_second_hit(s,j,:) = zeros(1,length(theta_phase)); 
        end      

        if nansum(perm_TwoTones_first_hit(s,j,:)) ~= 0       
          perm_bin_TwoTones_first_hit(s,j,:) = histcounts(perm_TwoTones_first_hit(s,j,:),-pi:bin:pi);
        else
          perm_bin_TwoTones_first_hit(s,j,:) = zeros(1,length(theta_phase)); 
        end

        if nansum(perm_TwoTones_second_hit(s,j,:)) ~= 0       
          perm_bin_TwoTones_second_hit(s,j,:) = histcounts(perm_TwoTones_second_hit(s,j,:),-pi:bin:pi);
        else
          perm_bin_TwoTones_second_hit(s,j,:) = zeros(1,length(theta_phase)); 
        end            

        if nansum(perm_OneTone_first_miss(s,j,:)) ~= 0       
          perm_bin_OneTone_first_miss(s,j,:) = histcounts(perm_OneTone_first_miss(s,j,:),-pi:bin:pi);
        else
          perm_bin_OneTone_first_miss(s,j,:) = zeros(1,length(theta_phase)); 
        end

        if nansum(perm_OneTone_second_miss(s,j,:)) ~= 0       
          perm_bin_OneTone_second_miss(s,j,:) = histcounts(perm_OneTone_second_miss(s,j,:),-pi:bin:pi);
        else
          perm_bin_OneTone_second_miss(s,j,:) = zeros(1,length(theta_phase)); 
        end      

        if nansum(perm_TwoTones_first_miss(s,j,:)) ~= 0       
          perm_bin_TwoTones_first_miss(s,j,:) = histcounts(perm_TwoTones_first_miss(s,j,:),-pi:bin:pi);
        else
          perm_bin_TwoTones_first_miss(s,j,:) = zeros(1,length(theta_phase)); 
        end

        if nansum(perm_TwoTones_second_miss(s,j,:)) ~= 0       
          perm_bin_TwoTones_second_miss(s,j,:) = histcounts(perm_TwoTones_second_miss(s,j,:),-pi:bin:pi);
        else
          perm_bin_TwoTones_second_miss(s,j,:) = zeros(1,length(theta_phase)); 
        end 

        [~,max_phase_hit] = max(perm_bin_OneTone_first_hit(s,j,:));
        max_toshift_hit = length(theta_phase)/2+1 - max_phase_hit;
        perm_maxPhase_realign.OneTone_first_hit(s,j,:) = circshift(perm_bin_OneTone_first_hit(s,j,:),max_toshift_hit);
        [~,max_phase_hit] = max(perm_bin_OneTone_second_hit(s,j,:));
        max_toshift_hit = length(theta_phase)/2+1 - max_phase_hit;
        perm_maxPhase_realign.OneTone_second_hit(s,j,:) = circshift(perm_bin_OneTone_second_hit(s,j,:),max_toshift_hit);
        [~,max_phase_hit] = max(perm_bin_TwoTones_first_hit(s,j,:));
        max_toshift_hit = length(theta_phase)/2+1 - max_phase_hit;
        perm_maxPhase_realign.TwoTones_first_hit(s,j,:) = circshift(perm_bin_TwoTones_first_hit(s,j,:),max_toshift_hit);
        [~,max_phase_hit] = max(perm_bin_TwoTones_second_hit(s,j,:));
        max_toshift_hit = length(theta_phase)/2+1 - max_phase_hit;
        perm_maxPhase_realign.TwoTones_second_hit(s,j,:) = circshift(perm_bin_TwoTones_second_hit(s,j,:),max_toshift_hit);
        [~,max_phase_miss] = max(perm_bin_OneTone_first_miss(s,j,:));
        max_toshift_miss = length(theta_phase)/2+1 - max_phase_miss;
        perm_maxPhase_realign.OneTone_first_miss(s,j,:) = circshift(perm_bin_OneTone_first_miss(s,j,:),max_toshift_miss);
        [~,max_phase_miss] = max(perm_bin_OneTone_second_miss(s,j,:));
        max_toshift_miss = length(theta_phase)/2+1 - max_phase_miss;
        perm_maxPhase_realign.OneTone_second_miss(s,j,:) = circshift(perm_bin_OneTone_second_miss(s,j,:),max_toshift_miss);
        [~,max_phase_miss] = max(perm_bin_TwoTones_first_miss(s,j,:));
        max_toshift_miss = length(theta_phase)/2+1 - max_phase_miss;
        perm_maxPhase_realign.TwoTones_first_miss(s,j,:) = circshift(perm_bin_TwoTones_first_miss(s,j,:),max_toshift_miss);
        [~,max_phase_miss] = max(perm_bin_TwoTones_second_miss(s,j,:));
        max_toshift_miss = length(theta_phase)/2+1 - max_phase_miss;
        perm_maxPhase_realign.TwoTones_second_miss(s,j,:) = circshift(perm_bin_TwoTones_second_miss(s,j,:),max_toshift_miss);
        clear max_phase_hit max_phase_miss max_toshift_hit max_toshift_miss max_phase_miss max_phase_miss max_toshift_miss max_toshift_miss 

        % compute phase bin entrainment hit/miss and size effect at first/second tones in OneTone and Two Tones conditions;
        centerbin = length(theta_phase)/2+1;
        max_phs_OneTone_first_hit = perm_maxPhase_realign.OneTone_first_hit(s,j,:);
        max_phs_OneTone_second_hit = perm_maxPhase_realign.OneTone_second_hit(s,j,:);
        max_phs_TwoTones_first_hit = perm_maxPhase_realign.TwoTones_first_hit(s,j,:);
        max_phs_TwoTones_second_hit = perm_maxPhase_realign.TwoTones_second_hit(s,j,:);
        max_phs_OneTone_first_miss = perm_maxPhase_realign.OneTone_first_miss(s,j,:);
        max_phs_OneTone_second_miss = perm_maxPhase_realign.OneTone_second_miss(s,j,:);
        max_phs_TwoTones_first_miss = perm_maxPhase_realign.TwoTones_first_miss(s,j,:);
        max_phs_TwoTones_second_miss = perm_maxPhase_realign.TwoTones_second_miss(s,j,:);  

        % hit-miss: OneTone first;
        tocomp1 = ((max_phs_OneTone_first_hit(:,centerbin-1)+ max_phs_OneTone_first_hit(:,centerbin+1))./2) - ((max_phs_OneTone_first_hit(:,2) + max_phs_OneTone_first_hit(:,end))./2);        
        tocomp2 = ((max_phs_OneTone_first_miss(:,centerbin-1)+ max_phs_OneTone_first_miss(:,centerbin+1))./2) - ((max_phs_OneTone_first_miss(:,2) + max_phs_OneTone_first_miss(:,end))./2);
        perm_diff_phs_OneTone_first(s,j,1) = tocomp1 - tocomp2;        
        % hit-miss: OneTone second;
        tocomp3 = ((max_phs_OneTone_second_hit(:,centerbin-1)+ max_phs_OneTone_second_hit(:,centerbin+1))./2) - ((max_phs_OneTone_second_hit(:,2) + max_phs_OneTone_second_hit(:,end))./2);        
        tocomp4 = ((max_phs_OneTone_second_miss(:,centerbin-1)+ max_phs_OneTone_second_miss(:,centerbin+1))./2) - ((max_phs_OneTone_second_miss(:,2) + max_phs_OneTone_second_miss(:,end))./2);
        perm_diff_phs_OneTone_second(s,j,1) = tocomp3 - tocomp4;
        % hit-miss: TwoTones first;
        tocomp5 = ((max_phs_TwoTones_first_hit(:,centerbin-1)+ max_phs_TwoTones_first_hit(:,centerbin+1))./2) - ((max_phs_TwoTones_first_hit(:,2) + max_phs_TwoTones_first_hit(:,end))./2);        
        tocomp6 = ((max_phs_TwoTones_first_miss(:,centerbin-1)+ max_phs_TwoTones_first_miss(:,centerbin+1))./2) - ((max_phs_TwoTones_first_miss(:,2) + max_phs_TwoTones_first_miss(:,end))./2);        
        perm_diff_phs_TwoTones_first(s,j,1) = tocomp5 - tocomp6;    
        % hit-miss: TwoTones second;
        tocomp7 = ((max_phs_TwoTones_second_hit(:,centerbin-1)+ max_phs_TwoTones_second_hit(:,centerbin+1))./2) - ((max_phs_TwoTones_second_hit(:,2) + max_phs_TwoTones_second_hit(:,end))./2);        
        tocomp8 = ((max_phs_TwoTones_second_miss(:,centerbin-1)+ max_phs_TwoTones_second_miss(:,centerbin+1))./2) - ((max_phs_TwoTones_second_miss(:,2) + max_phs_TwoTones_second_miss(:,end))./2);
        perm_diff_phs_TwoTones_second(s,j,1) = tocomp7 - tocomp8;    
        % effect size;
        perm_effect_size_TwoTones(s,j,1) = tocomp7 - tocomp8 - tocomp5 - tocomp6; 
        perm_effect_size_OneTone(s,j,1) = tocomp3 - tocomp4 - tocomp1 - tocomp2; 
        clear tocomp1 tocomp2 tocomp3 tocomp4 tocomp5 tocomp6 tocomp7 tocomp8 
        
    end    
catch
end  

    % Update waitbar and message;
    waitbar(s/length(subjects),f,sprintf('Realigning the phase of permutations'));

end

delete(f);

% Remove the bad subjects to exclude from statistic analyses;
bad_subjects = [1 4 25 26 27];
perm_diff_phs_OneTone_first(bad_subjects,:) = []; perm_diff_phs_OneTone_second(bad_subjects,:) = [];
perm_diff_phs_TwoTones_first(bad_subjects,:) = []; perm_diff_phs_TwoTones_second(bad_subjects,:) = [];
perm_effect_size_TwoTones(bad_subjects,:) = []; perm_effect_size_OneTone(bad_subjects,:) = [];

% run the ttest for each permutation pairs and sort pvalues ascend;
f = waitbar(0,'Processing the T-tests of all permutations','Name','Work in Progress');
perm_OneTone_p_max_adj = []; perm_TwoTones_p_max_adj = []; perm_EffectSize_p_max_adj = [];

for j = 1:size(perm_diff_phs_OneTone_first,2)
    
    %TwoTones condition;
    [~,perm_TwoTones_first_p_max_adj(j,2),~,stats1] = ttest(perm_diff_phs_TwoTones_first(:,j),0,'tail','right');
    perm_TwoTones_first_p_max_adj(j,1) = stats1.tstat;
    clear stats1    
    [~,perm_TwoTones_second_p_max_adj(j,2),~,stats2] = ttest(perm_diff_phs_TwoTones_second(:,j),0,'tail','right');
    perm_TwoTones_second_p_max_adj(j,1) = stats2.tstat;
    clear stats2    
    [~,perm_TwoTones_p_max_adj(j,2),~,stats3] = ttest2(perm_diff_phs_TwoTones_second(:,j),perm_diff_phs_TwoTones_first(:,j),'tail','right');
    perm_TwoTones_p_max_adj(j,1) = stats3.tstat;
    clear stat3    
    %OneTone condition;
    [~,perm_OneTone_first_p_max_adj(j,2),~,stats4] = ttest(perm_diff_phs_OneTone_first(:,j),0,'tail','right');
    perm_OneTone_first_p_max_adj(j,1) = stats4.tstat;
    clear stats4   
    [~,perm_OneTone_second_p_max_adj(j,2),~,stats5] = ttest(perm_diff_phs_OneTone_second(:,j),0,'tail','right');
    perm_OneTone_second_p_max_adj(j,1) = stats5.tstat;
    clear stats5       
    [~,perm_OneTone_p_max_adj(j,2),~,stats6] = ttest2(perm_diff_phs_OneTone_second(:,j),perm_diff_phs_OneTone_first(:,j),'tail','right');
    perm_OneTone_p_max_adj(j,1) = stats6.tstat;
    clear stat6
    [~,perm_EffectSize_p_max_adj(j,2),~,stats7] = ttest2(perm_effect_size_TwoTones(:,j),perm_effect_size_OneTone(:,j),'tail','right');
    perm_EffectSize_p_max_adj(j,1) = stats7.tstat;
    clear stat7
    
    %Update waitbar and message;
    waitbar(j/size(perm_diff_phs_OneTone_first,2),f,sprintf('Processing the T-tests of all permutations'));

end

delete(f);

% sort the tvalues of all the permutation ttests;
perm_OneTone_first_p_max_adj = sortrows(perm_OneTone_first_p_max_adj,1,'descend'); perm_OneTone_second_p_max_adj = sortrows(perm_OneTone_second_p_max_adj,1,'descend');
perm_TwoTones_first_p_max_adj = sortrows(perm_TwoTones_first_p_max_adj,1,'descend'); perm_TwoTones_second_p_max_adj = sortrows(perm_TwoTones_second_p_max_adj,1,'descend');
perm_OneTone_p_max_adj = sortrows(perm_OneTone_p_max_adj,1,'descend'); perm_TwoTones_p_max_adj = sortrows(perm_TwoTones_p_max_adj,1,'descend');
perm_EffectSize_p_max_adj = sortrows(perm_EffectSize_p_max_adj,1,'descend');

% calculate and display the p-values;
clc;
P_value_perm1 = (sum(perm_OneTone_first_p_max_adj(:,1) > ori_OneTone_first_p_max_adj(1,1))+1)/(length(perm_OneTone_first_p_max_adj(:,1))+1);
display(['OneTone-first: ' num2str(P_value_perm1)]); 
P_value_perm2 = (sum(perm_OneTone_second_p_max_adj(:,1) > ori_OneTone_second_p_max_adj(1,1))+1)/(length(perm_OneTone_second_p_max_adj(:,1))+1);
display(['OneTone-second: ' num2str(P_value_perm2)]);
P_value_perm5 = (sum(perm_OneTone_p_max_adj(:,1) > ori_OneTone_p_max_adj(1,1))+1)/(length(perm_OneTone_p_max_adj(:,1))+1);
display(['OneTone-Interaction: ' num2str(P_value_perm5)]); 
P_value_perm3 = (sum(perm_TwoTones_first_p_max_adj(:,1) > ori_TwoTones_first_p_max_adj(1,1))+1)/(length(perm_TwoTones_first_p_max_adj(:,1))+1);
display(['TwoTones-first: ' num2str(P_value_perm3)]); 
P_value_perm4 = (sum(perm_TwoTones_second_p_max_adj(:,1) > ori_TwoTones_second_p_max_adj(1,1))+1)/(length(perm_TwoTones_second_p_max_adj(:,1))+1);
display(['TwoTones-second: ' num2str(P_value_perm4)]);
P_value_perm6 = (sum(perm_TwoTones_p_max_adj(:,1) > ori_TwoTones_p_max_adj(1,1))+1)/(length(perm_TwoTones_p_max_adj(:,1))+1);
display(['TwoTones-Interaction: ' num2str(P_value_perm6)]); 
P_value_perm7 = (sum(perm_EffectSize_p_max_adj(:,1) > ori_EffectSize_p_max_adj(1,1))+1)/(length(perm_EffectSize_p_max_adj(:,1))+1);
display(['OnevsTwo Conditions: ' num2str(P_value_perm7)]); 

%% PLOT Hits-Miss entrainment at first and second tone in the OneTone and TwoTones conditions;

close all; figure; 
cond = [1 2];
count = 0;
for condition = cond

    count = count +1;
    subplot(2,1,count)
    
    % select condition data;
    if condition == 1
        aligned_phase_second = ori_maxPhase_realign.OneTone_second_hit./sum(ori_maxPhase_realign.OneTone_second_hit,2);
        aligned_phase_first = ori_maxPhase_realign.OneTone_first_hit./sum(ori_maxPhase_realign.OneTone_first_hit,2);
        title_plot = 'Single Tone Condition';
    elseif condition == 2
        aligned_phase_second = ori_maxPhase_realign.TwoTones_second_hit./sum(ori_maxPhase_realign.TwoTones_second_hit,2);
        aligned_phase_first = ori_maxPhase_realign.TwoTones_first_hit./sum(ori_maxPhase_realign.TwoTones_first_hit,2);
        title_plot = 'Two Tones Condition';
    end

    % Remove The bad subjects; 
    aligned_phase_second(bad_subjects,:) = [];
    aligned_phase_first(bad_subjects,:) = [];
        
    % close all;figure;
    x = 1:numel(aligned_phase_second(1,1:centerbin-1));
    x2 = [x, fliplr(x)];
    x3 = 1:numel(aligned_phase_second(1,centerbin+1:end));
    x4 = [x3, fliplr(x3)];
    phase_binned = 1:1:length(theta_phase);

    x_labels = [];
    for p1 = 1:length(phase_binned)

        max_phs = find(mean(aligned_phase_second) == max(mean(aligned_phase_second)));
        if  p1 == max_phs
            x_labels{1,p1} = 'max phase';
        elseif p1 == max_phs - 4
            x_labels{1,p1} = 'max-4'; 
        elseif p1 == max_phs - 7
            x_labels{1,p1} = 'max-7'; 
        elseif p1 == max_phs + 4
            x_labels{1,p1} = 'max+4';
        elseif p1 == max_phs + 7
            x_labels{1,p1} = 'max+7';
        else
            x_labels{1,p1} = '';
        end
    end

    % prepare plot for second tone;
    mean_phs_second = nanmean(aligned_phase_second(:,1:centerbin-1))/sum(nanmean(aligned_phase_second));
    std_phs_second = nanstd(mean_phs_second);%/sqrt(numel(aligned_phase_second(:,1)));
    curve1 = mean_phs_second + std_phs_second;
    curve2 = mean_phs_second - std_phs_second;
    inBetween_phs_second = [curve1, fliplr(curve2)];
    mean_phs2_second = nanmean(aligned_phase_second(:,centerbin+1:end))/sum(nanmean(aligned_phase_second));
    std_phs2_second = nanstd(mean_phs2_second);%/sqrt(numel(aligned_phase_second(:,1)));
    curve3 = mean_phs2_second + std_phs2_second;
    curve4 = mean_phs2_second - std_phs2_second;
    inBetween_phs2_second = [curve3, fliplr(curve4)];
    % prepare plot for first tone;
    mean_phs_first = nanmean(aligned_phase_first(:,1:centerbin-1))/sum(nanmean(aligned_phase_first));
    std_phs_first = nanstd(mean_phs_first);%/sqrt(numel(aligned_phase_first(:,1)));
    curve5 = mean_phs_first + std_phs_first;
    curve6 = mean_phs_first - std_phs_first;
    inBetween_phs_first = [curve5, fliplr(curve6)];
    mean_phs2_first = nanmean(aligned_phase_first(:,centerbin+1:end))/sum(nanmean(aligned_phase_first));
    std_phs2_first = nanstd(mean_phs2_first);%/sqrt(numel(aligned_phase_first(:,1)));
    curve7 = mean_phs2_first + std_phs2_first;
    curve8 = mean_phs2_first - std_phs2_first;
    inBetween_phs2_first = [curve7, fliplr(curve8)];


    f_phs_second = fill(x2, inBetween_phs_second, [0.8 0.8 0.8]);
    f_phs_second.EdgeColor= 'none';
    f_phs_second.FaceColor= [0.5 0 0];
    f_phs_second.FaceAlpha= 0.1; hold on
    f_phs2_second = fill(x4+centerbin, inBetween_phs2_second, [0.8 0.8 0.8]);
    f_phs2_second.EdgeColor= 'none';
    f_phs2_second.FaceColor= [0.5 0 0];
    f_phs2_second.FaceAlpha= 0.1; hold on
    f_phs_first = fill(x2, inBetween_phs_first, [0.8 0.8 0.8]);
    f_phs_first.EdgeColor= 'none';
    f_phs_first.FaceColor= [0.8 0.8 0.8];
    f_phs_first.FaceAlpha= 0.5; hold on
    f_phs2_first = fill(x4+centerbin, inBetween_phs2_first, [0.8 0.8 0.8]);
    f_phs2_first.EdgeColor= 'none';
    f_phs2_first.FaceColor= [0.8 0.8 0.8];
    f_phs2_first.FaceAlpha= 0.5; hold on

    p1 = plot(phase_binned(1:centerbin-1),mean_phs_second,'Color', [0.5 0 0], 'LineWidth', 2); hold on
    p2 = plot(phase_binned(centerbin+1:end),mean_phs2_second,'Color', [0.5 0 0], 'LineWidth', 2); hold on
    p3 = plot(phase_binned(1:centerbin-1),mean_phs_first,'Color', [0 0 0], 'LineWidth', 2); hold on
    p4 = plot(phase_binned(centerbin+1:end),mean_phs2_first,'Color', [0 0 0], 'LineWidth', 2); hold on

    scatter(12,mean_phs2_first(1),'Marker','o','SizeData',30,'MarkerEdgeColor','k','MarkerFaceColor', [0 0.6 1]); hold on
    scatter(20,mean_phs2_first(9),'Marker','o','SizeData',30,'MarkerEdgeColor','k','MarkerFaceColor','w'); hold on
    scatter(12,mean_phs2_second(1),'Marker','o','SizeData',30,'MarkerEdgeColor',[0.5 0 0],'MarkerFaceColor',[0 0.6 1]); hold on
    scatter(20,mean_phs2_second(9),'Marker','o','SizeData',30,'MarkerEdgeColor',[0.5 0 0],'MarkerFaceColor','w'); hold on
    scatter(2,mean_phs_first(2),'Marker','o','SizeData',30,'MarkerEdgeColor','k','MarkerFaceColor','w'); hold on
    scatter(10,mean_phs_first(10),'Marker','o','SizeData',30,'MarkerEdgeColor','k','MarkerFaceColor',[0 0.6 1]); hold on
    scatter(2,mean_phs_second(2),'Marker','o','SizeData',30,'MarkerEdgeColor',[0.5 0 0],'MarkerFaceColor','w'); hold on
    scatter(10,mean_phs_second(10),'Marker','o','SizeData',30,'MarkerEdgeColor',[0.5 0 0],'MarkerFaceColor',[0 0.6 1]);
    xlim([0 length(theta_phase)+1]);
    ylim([0.02 0.06]);
    xticks(phase_binned)
    xticklabels(x_labels)
    title(title_plot)
    ylabel('p(hit)')
    legend([p1 p3],'second tone','first tone','Location','SouthEast');

end

%suptitle('REALIGNED PHASE')
set(gcf,'color',[0.95 0.95 0.95])
