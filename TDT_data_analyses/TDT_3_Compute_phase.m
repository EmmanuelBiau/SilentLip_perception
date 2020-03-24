%% Now, we will calculate the phase information for every tones in single subjects (Hits and Misses);
clear;clc;
cd XXX\Signal_Infos\;
addpath XXX\Circular_Statistics_Toolbox\;
load('Phase_mina');

subjects = 1:29;

for s = subjects

cd(['XXX\subj' num2str(s)]);
load(['subj' num2str(s) '-Tone_detectionTask_clean']); 
fs = 44100;

%remove the first 4 trials of practice;
respMat_ANALYSES(1:4,:) = [];

%Create a new table with needed information for phase analyses; extract need info;
condition = cell2mat([respMat_ANALYSES.TONE(cell2mat(respMat_ANALYSES.TONE)==1);respMat_ANALYSES.TONE(cell2mat(respMat_ANALYSES.TONE)==2);respMat_ANALYSES.TONE(cell2mat(respMat_ANALYSES.TONE)==2)]);
position = [ones(length(respMat_ANALYSES.HITS_TONE1(cell2mat(respMat_ANALYSES.TONE)==1)),1);ones(length(respMat_ANALYSES.HITS_TONE1(cell2mat(respMat_ANALYSES.TONE)==2)),1);2*ones(length(respMat_ANALYSES.HITS_TONE2(cell2mat(respMat_ANALYSES.TONE)==2)),1)];
hits = cell2mat([respMat_ANALYSES.HITS_TONE1(cell2mat(respMat_ANALYSES.TONE)==1);respMat_ANALYSES.HITS_TONE1(cell2mat(respMat_ANALYSES.TONE)==2);respMat_ANALYSES.HITS_TONE2(cell2mat(respMat_ANALYSES.TONE)==2)]);
rt_1 = cell2mat(respMat_ANALYSES.REACTION_TIME1(cell2mat(respMat_ANALYSES.TONE)==1));
rt_1(rt_1 > 4.95) = 0;
reaction_time = [rt_1 - cell2mat(respMat_ANALYSES.TONE_ONSET1(cell2mat(respMat_ANALYSES.TONE)==1))/fs; cell2mat(respMat_ANALYSES.REACTION_TIME1(cell2mat(respMat_ANALYSES.TONE)==2)) - cell2mat(respMat_ANALYSES.TONE_ONSET1(cell2mat(respMat_ANALYSES.TONE)==2))/fs; cell2mat(respMat_ANALYSES.REACTION_TIME2(cell2mat(respMat_ANALYSES.TONE)==2)) - cell2mat(respMat_ANALYSES.TONE_ONSET2(cell2mat(respMat_ANALYSES.TONE) == 2))/fs];
reaction_time(reaction_time < 0 | reaction_time > 4.95) = 5;
onsets = cell2mat([respMat_ANALYSES.TONE_ONSET1(cell2mat(respMat_ANALYSES.TONE)==1);respMat_ANALYSES.TONE_ONSET1(cell2mat(respMat_ANALYSES.TONE)==2);respMat_ANALYSES.TONE_ONSET2(cell2mat(respMat_ANALYSES.TONE)==2)]);
video = [respMat_ANALYSES.VIDEO(cell2mat(respMat_ANALYSES.TONE)==1);respMat_ANALYSES.VIDEO(cell2mat(respMat_ANALYSES.TONE)==2);respMat_ANALYSES.VIDEO(cell2mat(respMat_ANALYSES.TONE)==2)];
freqpeak = cell2mat([respMat_ANALYSES.FREQPEAK(cell2mat(respMat_ANALYSES.TONE)==1);respMat_ANALYSES.FREQPEAK(cell2mat(respMat_ANALYSES.TONE)==2);respMat_ANALYSES.FREQPEAK(cell2mat(respMat_ANALYSES.TONE)==2)]);
lag = cell2mat([respMat_ANALYSES.VA_LAG(cell2mat(respMat_ANALYSES.TONE)==1);respMat_ANALYSES.VA_LAG(cell2mat(respMat_ANALYSES.TONE)==2);respMat_ANALYSES.VA_LAG(cell2mat(respMat_ANALYSES.TONE)==2)]);
participant = s*ones(length(hits),1);
trial_numb = cell2mat([respMat_ANALYSES.TRIAL_NUMB(cell2mat(respMat_ANALYSES.TONE)==1);respMat_ANALYSES.TRIAL_NUMB(cell2mat(respMat_ANALYSES.TONE)==2);respMat_ANALYSES.TRIAL_NUMB(cell2mat(respMat_ANALYSES.TONE)==2)]);

for tn = 1:length(trial_numb)
    
    if trial_numb(tn) <=150
        exp_half(tn,1) = 1;
    elseif trial_numb(tn) > 150 
        exp_half(tn,1) = 2;
    end
    
end

%feed the table;
PHASE_ANALYSES = cell(length(hits), 12);
PHASE_ANALYSES = array2table(PHASE_ANALYSES, 'Variable',{'PARTICIPANT','CONDITION','TONE_POSITION','HIT','REACTION_TIME','ONSET','VIDEO','LAG', 'FREQPEAK','THETA_PHASE','TRIAL_NUMB','EXP_HALF'});
PHASE_ANALYSES.PARTICIPANT = participant;
PHASE_ANALYSES.CONDITION = condition ;
PHASE_ANALYSES.TONE_POSITION = position;
PHASE_ANALYSES.HIT = hits;
PHASE_ANALYSES.REACTION_TIME = reaction_time;
PHASE_ANALYSES.ONSET = onsets;
PHASE_ANALYSES.VIDEO = video;
PHASE_ANALYSES.FREQPEAK = freqpeak;
PHASE_ANALYSES.LAG = lag;
PHASE_ANALYSES.TRIAL_NUMB = trial_numb-4;
PHASE_ANALYSES.EXP_HALF = exp_half;

%calculate the phase angle for each tone;
for i = 1:length(PHASE_ANALYSES.VIDEO)
    
    phase_video = cell2mat(Phase_mina.trial(PHASE_ANALYSES.VIDEO(i))); % phase of the video;
        
        %Audio onset lags: add the lag to the tone onset to get real phase time-point in visual information;
        if PHASE_ANALYSES.LAG(i) >= 0
            PHASE_ANALYSES.THETA_PHASE{i} = phase_video(round((PHASE_ANALYSES.ONSET(i)/fs - PHASE_ANALYSES.LAG(i))*250)+1);
        
        %Audio onset leads: remove the lag to the tone onset to get real phase time-point in visual information;
        elseif PHASE_ANALYSES.LAG(i) < 0
            PHASE_ANALYSES.THETA_PHASE{i} = phase_video(round((PHASE_ANALYSES.ONSET(i)/fs + PHASE_ANALYSES.LAG(i))*250)+1); 
        end
        
end

clear rt_1 reaction_time onsets video freqpeak participant hits position condition lag 

cd(['XXX\subj' num2str(s)]);
save (['subj' num2str(s) '-Phase_analyses'],'PHASE_ANALYSES','alltrials','participant_parameters');

clearvars -except Phase_mina

end

%%