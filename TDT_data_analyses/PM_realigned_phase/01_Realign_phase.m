%% First, we realign the visual phase on preferred bin for each participant;
clear;clc;
addpath C:\toolbox\CircHist-master; addpath C:\toolbox\Circular_Statistics_Toolbox\;  
cd C:\ANALYSES_DATA\analyses_realign_phase\;
load ONE_TONE_ANALYSES;load TWO_TONES_ANALYSES;

% Determine the number of bins in theta phase;
bin = pi/10;
theta_phase = -pi:bin:pi-bin;
centerbin = length(theta_phase)/2+1;

% subjects IDs;
subjects = 1:29;

% I- Compile the hits from real data of all participants together;
ori_OneTone_first_hit = []; ori_OneTone_second_hit = [];
ori_OneTone_first_miss = []; ori_OneTone_second_miss = [];
ori_TwoTones_first_hit = []; ori_TwoTones_second_hit = [];
ori_TwoTones_first_miss = []; ori_TwoTones_second_miss = [];

for s = subjects
    
    OneTone_tmp = ONE_TONE_ANALYSES(ismember(ONE_TONE_ANALYSES.PARTICIPANT,s),:); 
    TwoTones_tmp = TWO_TONES_ANALYSES(ismember(TWO_TONES_ANALYSES.PARTICIPANT,s),:); 

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

    clear max_phase_hit max_phase_miss max_toshift_hit max_toshift_miss max_phase_miss max_phase_miss max_toshift_miss max_toshift_miss 

    % compute phase bin entrainment hit/miss and size effect at first/second tones in OneTone and Two Tones conditions;
    ori_maxPhase_realign.OneTone_first_hit(s,:) = ori_maxPhase_realign.OneTone_first_hit(s,:)/(sum(ori_maxPhase_realign.OneTone_first_hit(s,:),2) + sum(ori_maxPhase_realign.OneTone_first_miss(s,:),2));
    ori_maxPhase_realign.OneTone_second_hit(s,:) = ori_maxPhase_realign.OneTone_second_hit(s,:)/(sum(ori_maxPhase_realign.OneTone_second_hit(s,:),2) + sum(ori_maxPhase_realign.OneTone_second_miss(s,:),2));
    ori_maxPhase_realign.TwoTones_first_hit(s,:) = ori_maxPhase_realign.TwoTones_first_hit(s,:)/(sum(ori_maxPhase_realign.TwoTones_first_hit(s,:),2) + sum(ori_maxPhase_realign.TwoTones_first_miss(s,:),2));
    ori_maxPhase_realign.TwoTones_second_hit(s,:) = ori_maxPhase_realign.TwoTones_second_hit(s,:)/(sum(ori_maxPhase_realign.TwoTones_second_hit(s,:),2) + sum(ori_maxPhase_realign.TwoTones_second_miss(s,:),2));
    ori_maxPhase_realign.OneTone_first_miss(s,:) = ori_maxPhase_realign.OneTone_first_miss(s,:)/(sum(ori_maxPhase_realign.OneTone_first_hit(s,:),2) + sum(ori_maxPhase_realign.OneTone_first_miss(s,:),2));
    ori_maxPhase_realign.OneTone_second_miss(s,:) = ori_maxPhase_realign.OneTone_second_miss(s,:)/(sum(ori_maxPhase_realign.OneTone_second_hit(s,:),2) + sum(ori_maxPhase_realign.OneTone_second_miss(s,:),2));
    ori_maxPhase_realign.TwoTones_first_miss(s,:) = ori_maxPhase_realign.TwoTones_first_miss(s,:)/(sum(ori_maxPhase_realign.TwoTones_first_hit(s,:),2) + sum(ori_maxPhase_realign.TwoTones_first_miss(s,:),2));
    ori_maxPhase_realign.TwoTones_second_miss(s,:) = ori_maxPhase_realign.TwoTones_second_miss(s,:)/(sum(ori_maxPhase_realign.TwoTones_second_hit(s,:),2) + sum(ori_maxPhase_realign.TwoTones_second_miss(s,:),2)); 
    
    centerbin = length(theta_phase)/2+1;
    max_phs_OneTone_first_hit = ori_maxPhase_realign.OneTone_first_hit(s,:);
    max_phs_OneTone_second_hit = ori_maxPhase_realign.OneTone_second_hit(s,:);
    max_phs_TwoTones_first_hit = ori_maxPhase_realign.TwoTones_first_hit(s,:);
    max_phs_TwoTones_second_hit = ori_maxPhase_realign.TwoTones_second_hit(s,:);
    max_phs_OneTone_first_miss = ori_maxPhase_realign.OneTone_first_miss(s,:);
    max_phs_OneTone_second_miss = ori_maxPhase_realign.OneTone_second_miss(s,:);
    max_phs_TwoTones_first_miss = ori_maxPhase_realign.TwoTones_first_miss(s,:);
    max_phs_TwoTones_second_miss = ori_maxPhase_realign.TwoTones_second_miss(s,:);  

    %hit-miss: OneTone first;
    tocomp1(s,1) = ((max_phs_OneTone_first_hit(:,centerbin-1)+ max_phs_OneTone_first_hit(:,centerbin+1))./2) - ((max_phs_OneTone_first_hit(:,2) + max_phs_OneTone_first_hit(:,end))./2);        
    tocomp2(s,1) = ((max_phs_OneTone_first_miss(:,centerbin-1)+ max_phs_OneTone_first_miss(:,centerbin+1))./2) - ((max_phs_OneTone_first_miss(:,2) + max_phs_OneTone_first_miss(:,end))./2);
    ori_diff_phs_OneTone_first(s,1) = tocomp1(s,1) - tocomp2(s,1);
    %hit-miss: OneTone second;
    tocomp3(s,1) = ((max_phs_OneTone_second_hit(:,centerbin-1)+ max_phs_OneTone_second_hit(:,centerbin+1))./2) - ((max_phs_OneTone_second_hit(:,2) + max_phs_OneTone_second_hit(:,end))./2);        
    tocomp4(s,1) = ((max_phs_OneTone_second_miss(:,centerbin-1)+ max_phs_OneTone_second_miss(:,centerbin+1))./2) - ((max_phs_OneTone_second_miss(:,2) + max_phs_OneTone_second_miss(:,end))./2);
    ori_diff_phs_OneTone_second(s,1) = tocomp3(s,1) - tocomp4(s,1);
    %hit-miss: TwoTones first;
    tocomp5(s,1) = ((max_phs_TwoTones_first_hit(:,centerbin-1)+ max_phs_TwoTones_first_hit(:,centerbin+1))./2) - ((max_phs_TwoTones_first_hit(:,2) + max_phs_TwoTones_first_hit(:,end))./2);        
    tocomp6(s,1) = ((max_phs_TwoTones_first_miss(:,centerbin-1)+ max_phs_TwoTones_first_miss(:,centerbin+1))./2) - ((max_phs_TwoTones_first_miss(:,2) + max_phs_TwoTones_first_miss(:,end))./2);        
    ori_diff_phs_TwoTones_first(s,1) = tocomp5(s,1) - tocomp6(s,1);    
    %hit-miss: TwoTones second;
    tocomp7(s,1) = ((max_phs_TwoTones_second_hit(:,centerbin-1)+ max_phs_TwoTones_second_hit(:,centerbin+1))./2) - ((max_phs_TwoTones_second_hit(:,2) + max_phs_TwoTones_second_hit(:,end))./2);        
    tocomp8(s,1) = ((max_phs_TwoTones_second_miss(:,centerbin-1)+ max_phs_TwoTones_second_miss(:,centerbin+1))./2) - ((max_phs_TwoTones_second_miss(:,2) + max_phs_TwoTones_second_miss(:,end))./2);
    ori_diff_phs_TwoTones_second(s,1) = tocomp7(s,1) - tocomp8(s,1);
    %effect size;
    ori_effect_size_TwoTones(s,1) = tocomp7(s,1) - tocomp8(s,1) - tocomp5(s,1) - tocomp6(s,1); 
    ori_effect_size_OneTone(s,1) = tocomp3(s,1) - tocomp4(s,1) - tocomp1(s,1) - tocomp2(s,1); 

    % clear tocomp1 tocomp2 tocomp3 tocomp4 tocomp5 tocomp6 tocomp7 tocomp8 

end

% Remove the bad subjects to exclude from statistic analyses;
bad_subjects = [1 4 25 26 27];
tocomp1(bad_subjects,:) = []; tocomp2(bad_subjects,:) = []; tocomp3(bad_subjects,:) = []; tocomp4(bad_subjects,:) = [];
tocomp5(bad_subjects,:) = []; tocomp6(bad_subjects,:) = []; tocomp7(bad_subjects,:) = []; tocomp8(bad_subjects,:) = [];
ori_diff_phs_OneTone_first(bad_subjects,:) = []; ori_diff_phs_OneTone_second(bad_subjects,:) = [];
ori_diff_phs_TwoTones_first(bad_subjects,:) = []; ori_diff_phs_TwoTones_second(bad_subjects,:) = [];
ori_effect_size_TwoTones(bad_subjects,:) = []; ori_effect_size_OneTone(bad_subjects,:) = [];
tocomp1(tocomp1==0,:) = 0.00001; tocomp2(tocomp2==0,:) = 0.00001; tocomp3(tocomp3==0,:) = 0.00001; tocomp4(tocomp4==0,:) = 0.00001;
tocomp5(tocomp5==0,:) = 0.00001; tocomp6(tocomp6==0,:) = 0.00001; tocomp7(tocomp7==0,:) = 0.00001; tocomp8(tocomp8==0,:) = 0.00001;
ori_diff_phs_OneTone_first(ori_diff_phs_OneTone_first== 0,:) = 0.00001; ori_diff_phs_OneTone_second(ori_diff_phs_OneTone_second== 0,:) = 0.00001;
ori_diff_phs_TwoTones_first(ori_diff_phs_TwoTones_first== 0,:) = 0.00001; ori_diff_phs_TwoTones_second(ori_diff_phs_TwoTones_second== 0,:) = 0.00001;
ori_effect_size_TwoTones(ori_effect_size_TwoTones== 0,:) = 0.00001; ori_effect_size_OneTone(ori_effect_size_OneTone== 0,:) = 0.00001;

% Store data for ANOVAs;
one_way_ANOVA_data = [];
one_way_ANOVA_data(:,1) = [ones(size(tocomp1,1),1);2*ones(size(tocomp3,1),1);3*ones(size(tocomp5,1),1); 4*ones(size(tocomp7,1),1)];
one_way_ANOVA_data(:,2) = [tocomp1;tocomp3;tocomp5;tocomp7];
two_way_ANOVA_data = [];
two_way_ANOVA_data(:,1) = tocomp1;
two_way_ANOVA_data(:,2) = tocomp3;
two_way_ANOVA_data(:,3) = tocomp5;
two_way_ANOVA_data(:,4) = tocomp7;

% T-tests (before correcting multiple comparisons);
Ttests_PM = [];
[~,Ttests_PM.OneTone_first(1,1),~,stat1] = ttest(tocomp1,0,'tail','right');
[~,Ttests_PM.OneTone_second(1,1),~,stat2] = ttest(tocomp3,0,'tail','right');
[~,Ttests_PM.TwoTones_first(1,1),~,stat3] = ttest(tocomp5,0,'tail','right');
[~,Ttests_PM.TwoTones_second(1,1),~,stat4] = ttest(tocomp7,0,'tail','right');
Ttests_PM.OneTone_first(1,2) = stat1.tstat;
Ttests_PM.OneTone_second(1,2) = stat2.tstat;
Ttests_PM.TwoTones_first(1,2) = stat3.tstat; 
Ttests_PM.TwoTones_second(1,2) = stat4.tstat;
Ttests_PM.Cohen_d = mean([tocomp1 tocomp3 tocomp5 tocomp7])./std([tocomp1 tocomp3 tocomp5 tocomp7]); 
clear stat1 stat2 stat3 stat4 tocomp1 tocomp2 tocomp3 tocomp4

%% PLOT Hits first vs second tone in the OneTone and TwoTones conditions;
close all; 

% Remove The bad subjects; 
if size(ori_maxPhase_realign.OneTone_first_hit,1) >  size(ori_diff_phs_OneTone_first,1) 
    ori_maxPhase_realign.OneTone_first_hit(bad_subjects,:) = [];
    ori_maxPhase_realign.OneTone_second_hit(bad_subjects,:) = [];
    ori_maxPhase_realign.TwoTones_first_hit(bad_subjects,:) = [];
    ori_maxPhase_realign.TwoTones_second_hit(bad_subjects,:) = [];
    ori_maxPhase_realign.OneTone_first_miss(bad_subjects,:) = [];
    ori_maxPhase_realign.OneTone_second_miss(bad_subjects,:) = [];
    ori_maxPhase_realign.TwoTones_first_miss(bad_subjects,:) = [];
    ori_maxPhase_realign.TwoTones_second_miss(bad_subjects,:) = [];
end

% condition (single tone:1; two tones:2);
cond = [1 2];

for cond_tone = cond

    subplot(2,1,cond_tone)
    
    if cond_tone == 1    
        % select condition data;
        aligned_phase_hit = mean(ori_maxPhase_realign.OneTone_first_hit);
        aligned_phase_hit2 = mean(ori_maxPhase_realign.OneTone_second_hit);
        title_plot = 'Single Tone';    
    elseif cond_tone == 2 
        % select condition data;
        aligned_phase_hit = mean(ori_maxPhase_realign.TwoTones_first_hit);
        aligned_phase_hit2 = mean(ori_maxPhase_realign.TwoTones_second_hit);
        title_plot = 'Two Tones';    
    end

    %close all;figure;
    x = 1:numel(aligned_phase_hit(1,1:centerbin-1));
    x2 = [x, fliplr(x)];
    x3 = 1:numel(aligned_phase_hit(1,centerbin+1:end));
    x4 = [x3, fliplr(x3)];
    phase_binned = 1:1:length(theta_phase);

    x_labels = [];
    for p1 = 1:length(phase_binned)

        max_phs = 11;
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

    % prepare plot for first tone;
    mean_phs_first = aligned_phase_hit(:,1:centerbin-1);
    std_phs_first = nanstd(mean_phs_first);%/sqrt(size(aligned_phase_hit(:,1),1));
    curve1 = mean_phs_first + std_phs_first;
    curve2 = mean_phs_first - std_phs_first;
    inBetween_phs_first = [curve1, fliplr(curve2)];
    mean_phs2_first = aligned_phase_hit(:,centerbin+1:end);
    std_phs2_first = nanstd(mean_phs2_first);%/sqrt(size(aligned_phase_hit(:,1),1));
    curve3 = mean_phs2_first + std_phs2_first;
    curve4 = mean_phs2_first - std_phs2_first;
    inBetween_phs2_first = [curve3, fliplr(curve4)];
    % prepare plot for second tone;
    mean_phs_second = aligned_phase_hit2(:,1:centerbin-1);
    std_phs_second = nanstd(mean_phs_second);%/sqrt(size(aligned_phase_hit(:,1),1));
    curve5 = mean_phs_second + std_phs_second;
    curve6 = mean_phs_second - std_phs_second;
    inBetween_phs_second = [curve5, fliplr(curve6)];
    mean_phs2_second = aligned_phase_hit2(:,centerbin+1:end);
    std_phs2_second = nanstd(mean_phs2_second);%/sqrt(size(aligned_phase_hit(:,1),1));
    curve7 = mean_phs2_second + std_phs2_second;
    curve8 = mean_phs2_second - std_phs2_second;
    inBetween_phs2_second = [curve7, fliplr(curve8)];

    % configure the std patches;
    f_phs_second = fill(x2, inBetween_phs_first, [0.8 0.8 0.8]);
    f_phs_second.EdgeColor= 'none';
    f_phs_second.FaceColor= [1 0 0];
    f_phs_second.FaceAlpha= 0.1; hold on
    f_phs2_second = fill(x4+centerbin, inBetween_phs2_first, [0.8 0.8 0.8]);
    f_phs2_second.EdgeColor= 'none';
    f_phs2_second.FaceColor= [1 0 0];
    f_phs2_second.FaceAlpha= 0.1; hold on
    f_phs_second2 = fill(x2, inBetween_phs_second, [0.8 0.8 0.8]);
    f_phs_second2.EdgeColor= 'none';
    f_phs_second2.FaceColor= [0 0.8 0];
    f_phs_second2.FaceAlpha= 0.1; hold on
    f_phs2_second2 = fill(x4+centerbin, inBetween_phs2_second, [0.8 0.8 0.8]);
    f_phs2_second2.EdgeColor= 'none';
    f_phs2_second2.FaceColor= [0 0.8 0];
    f_phs2_second2.FaceAlpha= 0.1; hold on

    % plot the phase modulations (before and aftre central bin);
    p1 = plot(phase_binned(1:centerbin-1),mean_phs_first,'Color', [1 0 0], 'LineWidth', 2); hold on
    p2 = plot(phase_binned(centerbin+1:end),mean_phs2_first,'Color', [1 0 0], 'LineWidth', 2); hold on
    p3 = plot(phase_binned(1:centerbin-1),mean_phs_second,'Color', [0 0.8 0], 'LineWidth', 2); hold on
    p4 = plot(phase_binned(centerbin+1:end),mean_phs2_second,'Color', [0 0.8 0], 'LineWidth', 2);hold on

    % mark important dots for visual;
    scatter(12,mean_phs2_first(1),'Marker','o','SizeData',30,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[0 0.6 1]); hold on
    scatter(20,mean_phs2_first(9),'Marker','o','SizeData',30,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor','w'); hold on
    scatter(2,mean_phs_first(2),'Marker','o','SizeData',30,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor','w'); hold on
    scatter(10,mean_phs_first(10),'Marker','o','SizeData',30,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[0 0.6 1]); hold on
    scatter(12,mean_phs2_second(1),'Marker','o','SizeData',30,'MarkerEdgeColor',[0 0.8 0],'MarkerFaceColor',[0 0.6 1]); hold on
    scatter(20,mean_phs2_second(9),'Marker','o','SizeData',30,'MarkerEdgeColor',[0 0.8 0],'MarkerFaceColor','w'); hold on
    scatter(2,mean_phs_second(2),'Marker','o','SizeData',30,'MarkerEdgeColor',[0 0.8 0],'MarkerFaceColor','w'); hold on
    scatter(10,mean_phs_second(10),'Marker','o','SizeData',30,'MarkerEdgeColor',[0 0.8 0],'MarkerFaceColor',[0 0.6 1]);

    xlim([0 length(theta_phase)+1]);
    ylim([0.01 0.06]);
    xticks(phase_binned)
    xticklabels(x_labels)
    title(title_plot)
    ylabel('p(hit)')

end

set(gcf,'color',[0.95 0.95 0.95],'Position',[777.8 133.8 416.8 547.2]);

%%