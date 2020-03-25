%%Compile the hits/miss and their reaction times;
clear;clc;
cd XXX\;
load allsubj_FA_sorted ; %FA information already computed;

%subjects ID;
subjects = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29];


%-------------------------------------------------------------------------------------%
%       Compute the performances for T1/T2 in the One Tone condition;
%-------------------------------------------------------------------------------------%

clear ONE_TONES_ANALYSES
count = 0;
for s = subjects
    
    cd(['XXX\subj' num2str(s)]);
    load(['subj' num2str(s) '-Phase_analyses']) 
    
    one_tone_temp = PHASE_ANALYSES(PHASE_ANALYSES.CONDITION ==1,:);
    
    for i = 1:length(one_tone_temp.TONE_POSITION)
        
            count = count+1;
            ONE_TONES_ANALYSES(count,1) = one_tone_temp.PARTICIPANT(i);
            ONE_TONES_ANALYSES(count,2) = one_tone_temp.ONSET(i);

                if ONE_TONES_ANALYSES(count,2)/44100 <= 2.5
                    ONE_TONES_ANALYSES(count,3) = 1;
                elseif ONE_TONES_ANALYSES(count,2)/44100 > 2.5
                    ONE_TONES_ANALYSES(count,3) = 2;
                end

            ONE_TONES_ANALYSES(count,4) = one_tone_temp.THETA_PHASE{i};
            ONE_TONES_ANALYSES(count,5) = one_tone_temp.FREQPEAK(i);
            ONE_TONES_ANALYSES(count,6) = one_tone_temp.HIT(i);
            ONE_TONES_ANALYSES(count,7) = one_tone_temp.REACTION_TIME(i);           
    
    end
    
    clear one_tone_temp    

end

ONE_TONES_ANALYSES = array2table(ONE_TONES_ANALYSES,'Variable',{'PARTICIPANT','OneTone_Onset','OneTone_position','OneTone_phase','Freqpeak','Hit', 'Rt'});
    
%compute the CRs and RTs for separate T1/T2 equivalent tones;
clear hit_sum_oneTone
s = 0;

for i = subjects
    
    s = s +1;
    hit_sum_oneTone(s,1) = sum(ONE_TONES_ANALYSES.Hit(ONE_TONES_ANALYSES.PARTICIPANT == i & ONE_TONES_ANALYSES.OneTone_position == 1))./sum(ONE_TONES_ANALYSES.PARTICIPANT == i & ONE_TONES_ANALYSES.OneTone_position == 1);
    hit_sum_oneTone(s,2) = sum(ONE_TONES_ANALYSES.Hit(ONE_TONES_ANALYSES.PARTICIPANT == i & ONE_TONES_ANALYSES.OneTone_position == 2))./sum(ONE_TONES_ANALYSES.PARTICIPANT == i & ONE_TONES_ANALYSES.OneTone_position == 2);       
    hit_sum_oneTone(s,3) = FA_rate.FA_T1_window(i);
    hit_sum_oneTone(s,4) = FA_rate.FA_T2_window(i);  
    hit_sum_oneTone(s,7) = mean(ONE_TONES_ANALYSES.Rt(ONE_TONES_ANALYSES.PARTICIPANT == i & ONE_TONES_ANALYSES.OneTone_position == 1 & ONE_TONES_ANALYSES.Hit ==1));
    hit_sum_oneTone(s,8) = mean(ONE_TONES_ANALYSES.Rt(ONE_TONES_ANALYSES.PARTICIPANT == i & ONE_TONES_ANALYSES.OneTone_position == 2 & ONE_TONES_ANALYSES.Hit ==1));    
           
end

%compute the dprimes;
hit_sum_oneTone(hit_sum_oneTone(:,:)==0) = 0.0001;
hit_sum_oneTone(hit_sum_oneTone(:,:)==1) = 0.9999;
hit_sum_oneTone(:,5) = dprime_simple(hit_sum_oneTone(:,1),hit_sum_oneTone(:,3));
hit_sum_oneTone(:,6) = dprime_simple(hit_sum_oneTone(:,2),hit_sum_oneTone(:,4));

OneTone_T1vsT2_performances = cell(length(subjects), 9);  
OneTone_T1vsT2_performances = array2table(OneTone_T1vsT2_performances,'Variable',{'PARTICIPANT','hit_T1','hit_T2','fa_T1','fa_T2','dp_T1','dp_T2','rt_T1','rt_T2'});
OneTone_T1vsT2_performances.PARTICIPANT = subjects';
OneTone_T1vsT2_performances.hit_T1 = hit_sum_oneTone(:,1);  
OneTone_T1vsT2_performances.hit_T2 = hit_sum_oneTone(:,2); 
OneTone_T1vsT2_performances.fa_T1 = hit_sum_oneTone(:,3);  
OneTone_T1vsT2_performances.fa_T2 = hit_sum_oneTone(:,4);  
OneTone_T1vsT2_performances.dp_T1 = hit_sum_oneTone(:,5);
OneTone_T1vsT2_performances.dp_T2 = hit_sum_oneTone(:,6);
OneTone_T1vsT2_performances.rt_T1 = hit_sum_oneTone(:,7);
OneTone_T1vsT2_performances.rt_T2 = hit_sum_oneTone(:,8);


%-------------------------------------------------------------------------------------%
%       Compute the performances for T1/T2 in the Two Tones condition;
%-------------------------------------------------------------------------------------%

clear TWO_TONES_ANALYSES
count = 0;
for s = subjects
    
    cd(['D:\Birmingham_2018\BEHAVIOR_project\ANALYSES_DATA\subj' num2str(s)]);
    load(['subj' num2str(s) '-Phase_analyses']) 
    
    two_tones_temp = PHASE_ANALYSES(PHASE_ANALYSES.CONDITION ==2,:);
    
    for i = 1:length(two_tones_temp.TONE_POSITION)
        
            count = count+1;
            TWO_TONES_ANALYSES(count,1) = two_tones_temp.PARTICIPANT(i);
            TWO_TONES_ANALYSES(count,2) = two_tones_temp.ONSET(i);
            TWO_TONES_ANALYSES(count,3) = two_tones_temp.TONE_POSITION(i);
            TWO_TONES_ANALYSES(count,4) = two_tones_temp.THETA_PHASE{i};
            TWO_TONES_ANALYSES(count,5) = two_tones_temp.FREQPEAK(i);
            TWO_TONES_ANALYSES(count,6) = two_tones_temp.HIT(i);
            TWO_TONES_ANALYSES(count,7) = two_tones_temp.REACTION_TIME(i);           
    
    end
    
    clear two_tones_temp 

end

TWO_TONES_ANALYSES = array2table(TWO_TONES_ANALYSES,'Variable',{'PARTICIPANT','Tone_Onset','Tone_position','Tone_phase','Freqpeak','Hit', 'Rt'});
    
%compute the CRs and RTs for separate T1/T2 tones;
clear hit_sum_twoTones
s = 0;

for i = subjects
    
    s = s +1;
    hit_sum_twoTones(s,1) = sum(TWO_TONES_ANALYSES.Hit(TWO_TONES_ANALYSES.PARTICIPANT == i & TWO_TONES_ANALYSES.Tone_position == 1))./sum(TWO_TONES_ANALYSES.PARTICIPANT == i & TWO_TONES_ANALYSES.Tone_position == 1);
    hit_sum_twoTones(s,2) = sum(TWO_TONES_ANALYSES.Hit(TWO_TONES_ANALYSES.PARTICIPANT == i & TWO_TONES_ANALYSES.Tone_position == 2))./sum(TWO_TONES_ANALYSES.PARTICIPANT == i & TWO_TONES_ANALYSES.Tone_position == 2);       
    hit_sum_twoTones(s,3) = FA_rate.FA_T1_window(i);
    hit_sum_twoTones(s,4) = FA_rate.FA_T2_window(i); 
    hit_sum_twoTones(s,7) = mean(TWO_TONES_ANALYSES.Rt(TWO_TONES_ANALYSES.PARTICIPANT == i & TWO_TONES_ANALYSES.Tone_position == 1 & TWO_TONES_ANALYSES.Hit ==1));
    hit_sum_twoTones(s,8) = mean(TWO_TONES_ANALYSES.Rt(TWO_TONES_ANALYSES.PARTICIPANT == i & TWO_TONES_ANALYSES.Tone_position == 2 & TWO_TONES_ANALYSES.Hit ==1));    
           
end

%compute the dprimes;
hit_sum_twoTones(hit_sum_twoTones(:,:)==0) = 0.0001;
hit_sum_twoTones(hit_sum_twoTones(:,:)==1) = 0.9999;
hit_sum_twoTones(:,5) = dprime_simple(hit_sum_twoTones(:,1),hit_sum_twoTones(:,3));
hit_sum_twoTones(:,6) = dprime_simple(hit_sum_twoTones(:,2),hit_sum_twoTones(:,4));

TwoTones_T1vsT2_performances = cell(length(subjects), 9);  
TwoTones_T1vsT2_performances = array2table(TwoTones_T1vsT2_performances,'Variable',{'PARTICIPANT','hit_T1','hit_T2','fa_T1','fa_T2','dp_T1','dp_T2','rt_T1','rt_T2'});
TwoTones_T1vsT2_performances.PARTICIPANT = subjects';
TwoTones_T1vsT2_performances.hit_T1 = hit_sum_twoTones(:,1);  
TwoTones_T1vsT2_performances.hit_T2 = hit_sum_twoTones(:,2); 
TwoTones_T1vsT2_performances.fa_T1 = hit_sum_twoTones(:,3);  
TwoTones_T1vsT2_performances.fa_T2 = hit_sum_twoTones(:,4);  
TwoTones_T1vsT2_performances.dp_T1 = hit_sum_twoTones(:,5);
TwoTones_T1vsT2_performances.dp_T2 = hit_sum_twoTones(:,6);
TwoTones_T1vsT2_performances.rt_T1 = hit_sum_twoTones(:,7);
TwoTones_T1vsT2_performances.rt_T2 = hit_sum_twoTones(:,8);

% %save all participants performances;
% save allsubj_T1vsT2_performances OneTone_T1vsT2_performances TwoTones_T1vsT2_performances 

%% Exclude outliers based on Grand averaged performances (i.e. blind to any condition);  

%design final matrix;
MEAN_HIT = cell(length(subjects), 3);
MEAN_HIT = array2table(MEAN_HIT, 'Variable',{'PARTICIPANT','GA_HITS','GA_RT'});


for s = subjects
    
    cd(['XXX\subj' num2str(s)]);
    load(['subj' num2str(s) '-Phase_analyses']) 

    MEAN_HIT.PARTICIPANT{s} = s;
    MEAN_HIT.GA_HITS{s} = mean(PHASE_ANALYSES.HIT);
    MEAN_HIT.GA_RT{s} = mean(PHASE_ANALYSES.REACTION_TIME(PHASE_ANALYSES.HIT ==1));
 
    clear alltrials participant_parameters PHASE_ANALYSES

end


for i = 1: length(MEAN_HIT.PARTICIPANT)
    
    %exclude based on absurde GA Hits rates (<50% chance level or 100% perfect);
    if MEAN_HIT.GA_HITS{i} < 0.5 || MEAN_HIT.GA_HITS{i} == 1   
        MEAN_HIT.Var4{i} = 0;
    else
        MEAN_HIT.Var4{i} = 1;
    end
    
    %exclude based on RT faster or slower than GA Rts +/- 2 standard deviations;
    if MEAN_HIT.GA_RT{i} < mean(cell2mat(MEAN_HIT.GA_RT))- 2*std(cell2mat(MEAN_HIT.GA_RT)) || MEAN_HIT.GA_RT{i} > mean(cell2mat(MEAN_HIT.GA_RT))+ 2*std(cell2mat(MEAN_HIT.GA_RT))   
        MEAN_HIT.Var5{i} = 0;
    else
        MEAN_HIT.Var5{i} = 1;
    end

    if sum(MEAN_HIT.Var4{i} + MEAN_HIT.Var5{i}) == 2 
        MEAN_HIT.Var6{i} = 1;
    else
        MEAN_HIT.Var6{i} = 0;    
    end
    
end

%keep the bad subject inf0 (#1 was an experimenter so put it too);
bad_subjects = [1 cell2mat(MEAN_HIT.PARTICIPANT(cell2mat(MEAN_HIT.Var6) == 0))'];

%% PLOT dprimes;

%remove bad subjects;
if length(OneTone_T1vsT2_performances.PARTICIPANT) == 29 && length(TwoTones_T1vsT2_performances.PARTICIPANT) == 29     
    OneTone_T1vsT2_performances(bad_subjects,:) = [];
    TwoTones_T1vsT2_performances(bad_subjects,:) = [];
end

%select the columns of dprimes;
dat = [table2array(OneTone_T1vsT2_performances(:,6:7)) table2array(TwoTones_T1vsT2_performances(:,6:7))];          
col = [0 0 0];                        
dotsize = 15;                              
smooth = 200;    

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
labels = {'\fontsize{12}T1','\fontsize{17}  \fontsize{14}\bfSingle Tone','\fontsize{12}T2','\fontsize{12}T1','\fontsize{17}  \fontsize{14}\bfTwo Tones','\fontsize{12}T2'};
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
x1 = linspace(1.1,1.7);
y1 = linspace(6.3,6.3);
plot(x1,y1, 'Color', [0 0 0],'LineWidth',0.5);
text(1.3,6.4,'*','FontSize', 20,'FontWeight','light');
set(gcf,'color','w','Positio',[727 375 432 517]);

%% Plot Reaction Times;

%select the columns of dprimes;
dat = [table2array(OneTone_T1vsT2_performances(:,8:9)) table2array(TwoTones_T1vsT2_performances(:,8:9))].*1000;         
col = [0 0 0];                        
dotsize = 15;                              
smooth = 200;    

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
labels = {'\fontsize{12}T1','\fontsize{17}  \fontsize{14}\bfSingle Tone','\fontsize{12}T2','\fontsize{12}T1','\fontsize{17}  \fontsize{14}\bfTwo Tones','\fontsize{12}T2'};
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