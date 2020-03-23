close all; clear all; clc; cd C:\experiments\Manu\BEHAVIOR_Experiment\;
AssertOpenGL;GetSecs;WaitSecs(0.1);
%Initialise LabJack
lab_init; 
%Setup PTB with some default values
PsychDefaultSetup(2);
InitializePsychSound(1);

%----------------------------------------------------------------------------------------------------------
%   Set up the audio
%----------------------------------------------------------------------------------------------------------
reqlatencyclass = 4;
freq = 44100;
buffersize = 0;
latbias = [];
pahandle = PsychPortAudio('Open', [], [], reqlatencyclass, freq, 2);
prelat = PsychPortAudio('LatencyBias', pahandle, latbias);
postlat = PsychPortAudio('LatencyBias', pahandle);
%load white noise;
load('whitenoise'); %(wnseq = our white noise);
wnseq = (wnseq - min(wnseq))*(1-(-1))/(max(wnseq) - min(wnseq)) + (-1);
wnseq(2,:) = wnseq(1,:); 

%Send the trigger;
lab_put_code(L,1);

%play white noise;
PsychPortAudio('FillBuffer', pahandle, wnseq);
PsychPortAudio('Start', pahandle, 1, 0, 1);

%Close LabJack2
lab_close; 


%%




