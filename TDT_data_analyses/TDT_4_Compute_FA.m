%%Compile the hits/miss and their reaction times;
clearvars;clc;
cd XXX\;
fs = 44100;

%subjects ID;
subjects = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29];


for s = subjects
    
    cd(['XXX\subj' num2str(s)]);
    load(['subj' num2str(s) '-Phase_analyses']);
    
    %compile the first tone onsets;
    TEMP_first = PHASE_ANALYSES;
    TEMP_first(TEMP_first.CONDITION ~= 2 | TEMP_first.TONE_POSITION ~= 1,:) = [];
    TEMP_first(:,[2 3 5 7:9]) = [];
    
    for j = 1:length(TEMP_first{:,1})
        MEAN_ONSETS.two_tones_first{(s-1)*100+j,1} = table2array(TEMP_first(j,1));
        MEAN_ONSETS.two_tones_first{(s-1)*100+j,2} = table2array(TEMP_first(j,2));
        MEAN_ONSETS.two_tones_first{(s-1)*100+j,3} = table2array(TEMP_first(j,3));
        MEAN_ONSETS.two_tones_first{(s-1)*100+j,4} = cell2mat(table2array(TEMP_first(j,4)));
    end
        
    %compile the second tone onsets;
    TEMP_second = PHASE_ANALYSES;
    TEMP_second(TEMP_second.CONDITION ~= 2 | TEMP_second.TONE_POSITION ~= 2,:) = [];
    TEMP_second(:,[2 3 5 7:9]) = [];
    
    for j = 1:length(TEMP_second{:,1})
        MEAN_ONSETS.two_tones_second{(s-1)*100+j,1} = table2array(TEMP_second(j,1));
        MEAN_ONSETS.two_tones_second{(s-1)*100+j,2} = table2array(TEMP_second(j,2));
        MEAN_ONSETS.two_tones_second{(s-1)*100+j,3} = table2array(TEMP_second(j,3));
        MEAN_ONSETS.two_tones_second{(s-1)*100+j,4} = cell2mat(table2array(TEMP_second(j,4)));
    end
   
    %compile the one tone onsets;
    TEMP_one = PHASE_ANALYSES;
    TEMP_one(TEMP_one.CONDITION ~= 1,:)= [];
    TEMP_one(:,[2 3 5 7:9]) = [];
    
    for j = 1:length(TEMP_one{:,1})
        MEAN_ONSETS.one_tone{(s-1)*100+j,1} = table2array(TEMP_one(j,1));
        MEAN_ONSETS.one_tone{(s-1)*100+j,2} = table2array(TEMP_one(j,2));
        MEAN_ONSETS.one_tone{(s-1)*100+j,3} = table2array(TEMP_one(j,3));
        MEAN_ONSETS.one_tone{(s-1)*100+j,4} = cell2mat(table2array(TEMP_one(j,4)));           
    end     
   
    %compile the GA onsets;
    TEMP_GA = PHASE_ANALYSES;
    TEMP_GA(TEMP_GA.CONDITION == 0,:)= [];
    TEMP_GA(:,[2 3 5 7:9]) = [];
    
    for j = 1:length(TEMP_GA{:,1})
        MEAN_ONSETS.GA{(s-1)*300+j,1} = table2array(TEMP_GA(j,1));
        MEAN_ONSETS.GA{(s-1)*300+j,2} = table2array(TEMP_GA(j,2));
        MEAN_ONSETS.GA{(s-1)*300+j,3} = table2array(TEMP_GA(j,3));
        MEAN_ONSETS.GA{(s-1)*100+j,4} = cell2mat(table2array(TEMP_GA(j,4)));             
    end  
        
        
    clear alltrials participant_parameters PHASE_ANALYSES TEMP_first TEMP_second TEMP_one TEMP_GA j
    
end

%remove missing rows from participant 3;
MEAN_ONSETS.GA(791:900,:) = [];
MEAN_ONSETS.one_tone(261:300,:) = [];
MEAN_ONSETS.two_tones_first(266:300,:) = [];
MEAN_ONSETS.two_tones_second(266:300,:) = [];

%remove bad subjects;
bad_subjects = [4 12 25 26 27];
idx1 = ismember(cell2mat(MEAN_ONSETS.GA(:,1)),bad_subjects);
idx2 = ismember(cell2mat(MEAN_ONSETS.one_tone(:,1)),bad_subjects);
idx3 = ismember(cell2mat(MEAN_ONSETS.two_tones_first(:,1)),bad_subjects);
idx4 = ismember(cell2mat(MEAN_ONSETS.two_tones_second(:,1)),bad_subjects);
MEAN_ONSETS.GA(idx1,:) = []; 
MEAN_ONSETS.one_tone(idx2,:) = []; 
MEAN_ONSETS.two_tones_first(idx3,:) = []; 
MEAN_ONSETS.two_tones_second(idx4,:) = []; 

clear idx1 idx2 idx3 idx4 

%Calculate relevant infos about tones onsets in both conditions;
MEAN_ONSETS.min_two_tones_first = min(cell2mat(MEAN_ONSETS.two_tones_first(:,3)))/fs;
MEAN_ONSETS.max_two_tones_first = max(cell2mat(MEAN_ONSETS.two_tones_first(:,3)))/fs;
MEAN_ONSETS.min_two_tones_second = min(cell2mat(MEAN_ONSETS.two_tones_second(:,3)))/fs;
MEAN_ONSETS.max_two_tones_second = max(cell2mat(MEAN_ONSETS.two_tones_second(:,3)))/fs;
MEAN_ONSETS.min_one_tone = min(cell2mat(MEAN_ONSETS.one_tone(:,3)))/fs;
MEAN_ONSETS.max_one_tone = max(cell2mat(MEAN_ONSETS.one_tone(:,3)))/fs;
MEAN_ONSETS.min_GA = min(cell2mat(MEAN_ONSETS.GA(:,3)))/fs;
MEAN_ONSETS.max_GA = max(cell2mat(MEAN_ONSETS.GA(:,3)))/fs;
MEAN_ONSETS.mean_two_tones_first = mean(cell2mat(MEAN_ONSETS.two_tones_first(:,3)))/fs;
MEAN_ONSETS.std_two_tones_first = std(cell2mat(MEAN_ONSETS.two_tones_first(:,3)))/fs;
MEAN_ONSETS.mean_two_tones_second = mean(cell2mat(MEAN_ONSETS.two_tones_second(:,3)))/fs;
MEAN_ONSETS.std_two_tones_second = std(cell2mat(MEAN_ONSETS.two_tones_second(:,3)))/fs;
MEAN_ONSETS.mean_one_tone = mean(cell2mat(MEAN_ONSETS.one_tone(:,3)))/fs;
MEAN_ONSETS.std_one_tone = std(cell2mat(MEAN_ONSETS.one_tone(:,3)))/fs;
MEAN_ONSETS.mean_GA = mean(cell2mat(MEAN_ONSETS.GA(:,3)))/fs;
MEAN_ONSETS.std_GA = std(cell2mat(MEAN_ONSETS.GA(:,3)))/fs;
MEAN_ONSETS.tone_windows = MEAN_ONSETS.mean_two_tones_second - (MEAN_ONSETS.mean_two_tones_second - MEAN_ONSETS.mean_two_tones_first)/2;


%--------------------------------------------------------------------------%
%                        Calculate False Alarms;
%--------------------------------------------------------------------------%

FALSE_ALARMS = cell(length(subjects)*100, 4);
FALSE_ALARMS = array2table(FALSE_ALARMS, 'Variable',{'PARTICIPANT','FA','FA_Reaction_Time','Trial_numb'});

for s = subjects

    cd(['XXX\subj' num2str(s)]);
    load(['subj' num2str(s) '-Tone_detectionTask_clean']); 
    
    %remove the first 4 trials of practice;
    respMat_ANALYSES(1:4,:) = [];
    FA_TEMP = respMat_ANALYSES;
    FA_TEMP(cell2mat(FA_TEMP.TONE)==1 | cell2mat(FA_TEMP.TONE)==2, :) = [];

    for j = 1:length(FA_TEMP.PARTICIPANT)
        FALSE_ALARMS.PARTICIPANT{(s-1)*length(FA_TEMP.PARTICIPANT)+j} = s;
        FALSE_ALARMS.FA{(s-1)*length(FA_TEMP.PARTICIPANT) +j} = 1;
        FALSE_ALARMS.FA_Reaction_Time{(s-1)*length(FA_TEMP.PARTICIPANT)+j} = FA_TEMP.REACTION_TIME1{j};
        FALSE_ALARMS.Trial_numb{(s-1)*length(FA_TEMP.PARTICIPANT)+j} = FA_TEMP.TRIAL_NUMB{j};
    end
    
    clear respMat_ANALYSES FA_TEMP alltrials participant_parameters
    
end

FALSE_ALARMS(214:300,:)= [];
FALSE_ALARMS(isnan(cell2mat(FALSE_ALARMS.FA_Reaction_Time)),:) = [];
FALSE_ALARMS.FA_Reaction_Time = cell2mat(FALSE_ALARMS.FA_Reaction_Time)*1000;

temp_fa_window = zeros(1,1);    
    
for i = 1:length(FALSE_ALARMS.PARTICIPANT)
    
    if FALSE_ALARMS.FA_Reaction_Time(i)/1000 <  MEAN_ONSETS.tone_windows
        temp_fa_window(i,1) = 1;
    else
        temp_fa_window(i,1) = 2;
    end
end

WINDOW = temp_fa_window;
FALSE_ALARMS = addvars(FALSE_ALARMS,WINDOW,'After','FA_Reaction_Time');
clear i temp_fa_window WINDOW 

%Compile Correct responses and False alarms together;
count_FA = zeros(length(subjects),7);
temp_count = zeros(1,1);

for i =1:length(FALSE_ALARMS.WINDOW)
    
    temp_count(i,1) = cell2mat(FALSE_ALARMS.PARTICIPANT(i));  
    
    if FALSE_ALARMS.WINDOW(i)== 1
        
            temp_count(i,2) = 1;
            temp_count(i,3) = 0;
            temp_count(i,4) = FALSE_ALARMS.FA_Reaction_Time(i); 
            temp_count(i,5) = nan;            
            temp_count(i,6) = cell2mat(FALSE_ALARMS.Trial_numb(i));
            
    elseif FALSE_ALARMS.WINDOW(i)== 2 
        
            temp_count(i,2) = 0;
            temp_count(i,3) = 1;
            temp_count(i,4) = nan;
            temp_count(i,5) = FALSE_ALARMS.FA_Reaction_Time(i);  
            temp_count(i,6) = cell2mat(FALSE_ALARMS.Trial_numb(i));
    end
end

for i= 1:length(temp_count(:,4))
    
    if temp_count(i,6) <= 150
        temp_count(i,7) = 1;
    elseif temp_count(i,6) > 150
        temp_count(i,7) = 2;
    end
    
end

FA_sorted = array2table(temp_count,'Variable',{'PARTICIPANT','FA_T1_window','FA_T2_window','FA_T1_rt','FA_T2_rt','Trial_numb','Exp_half'});

%now calculate the FA rates for each participants for later dprimes;
clear fa_rate_temp
count = 0;
for i = 1:length(subjects)
    
    count = count +1;
    fa_rate_temp(i,1) = i;
    fa_rate_temp(i,2) = sum(FA_sorted.FA_T1_window(FA_sorted.PARTICIPANT ==i))/100;
    fa_rate_temp(i,3) = sum(FA_sorted.FA_T2_window(FA_sorted.PARTICIPANT ==i))/100;

end

FA_rate = array2table(fa_rate_temp,'Variable',{'PARTICIPANT','FA_T1_window','FA_T2_window'});

%save the FA info for all the subjects;
cd XXX\;
save allsubj_FA_sorted FA_sorted FA_rate

%%