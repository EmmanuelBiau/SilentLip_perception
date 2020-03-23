%Adapt paths - Psychotoolbox - triggersender to your machine;
close all; clear all; clc; 
cd XXX\; addpath XXX\;
addpath('XXX\Psychtoolbox\MatlabWindowsFilesR2007a\');
addpath(genpath('XXX\Psychtoolbox\'));
AssertOpenGL;
%Force GetSecs and WaitSecs into memory to avoid latency later on
GetSecs;
WaitSecs(0.1);


% %Initialise LabJack
% lab_init; 

%Setup PTB with some default values
PsychDefaultSetup(2);
InitializePsychSound(1);

%----------------------------------------------------------------------------------------------------------
%   Set up the audio
%----------------------------------------------------------------------------------------------------------
%Request latency mode 4: request the most aggressive settings for the
%given device fail if device can't meet the strictest requirements.
reqlatencyclass = 4;
%Requested output frequency. Must set this. 96000, 48000, or 44100. Ours are 44100
freq = 44100;
%Pointless to set this. Auto-selected to be optimal.
buffersize = 0;
%Best left alone, only here as manual override in case all the auto-tuning cleverness fails
%suggestedLatencySecs = 0.002;
%Set based on work with oscilliscope
latbias = [];
%ASIO4ALL v2. May need to change if new hardware is installed sca
%deviceid = 14; % PC

%Open audio device for low-latency output
%pahandle = PsychPortAudio('Open', deviceid, [], reqlatencyclass, freq, 2, buffersize, suggestedLatencySecs);
pahandle = PsychPortAudio('Open', [], [], reqlatencyclass, freq, 2);

%Tell driver about hardware's inherent latency, determined via calibration once
prelat = PsychPortAudio('LatencyBias', pahandle, latbias);
postlat = PsychPortAudio('LatencyBias', pahandle);

%Generate a beep at 1000Hz, 0.1 sec, 50% amplitude, to avoid latency later on
mynoise(1,:) = 0.5*MakeBeep(1000, 0.1, freq);
mynoise(2,:) = mynoise(1,:);
%Fill buffer with audio data
PsychPortAudio('FillBuffer', pahandle, mynoise);
%Perform one audio warmup trial, to get the sound hardware fully up and running, performing whatever lazy
%initialization only happens at real first use. This useless warmup will allow for lower latency for start
%of playback during actual use of the audio driver in the real trials
PsychPortAudio('Start', pahandle, 1, 0, 1);
PsychPortAudio('Stop', pahandle, 1);
%----------------------------------------------------------------------------------------------------------
%   Ok so now the audio hardware is fully initialised and our driver is on hot-standby, ready to start playback
%   of any sound with minimal latency
%----------------------------------------------------------------------------------------------------------

%Create dialog box to gather participant number
participant_num = inputdlg('Participant Number', 'Participant Number');
participant_num = str2double(participant_num);

%%

try


%HideCursor;
    
%----------------------------------------------------------------------------------------------------------
%   Set up the display and timing
%----------------------------------------------------------------------------------------------------------
    % select external screen
    screenid = max(Screen('Screens'));
    % create a window
    white = WhiteIndex(screenid);
    black = BlackIndex(screenid);
    [win, winRect] = PsychImaging('OpenWindow', screenid, 0.5);
    % define size of x and y based on screen size
    [screenXpixels, screenYpixels] = Screen('WindowSize', win);
    % define screen center based on screen size
    [xCenter, yCenter] = RectCenter(winRect);
    % set test font size
    Screen('TextSize', win, 20);
    % set text font
    Screen('TextFont', win, 'Helvetica');
    
    % get the flip interval, based on screen refresh rate
    ifi = Screen('GetFlipInterval', win);
    % Set based on work with oscilliscope
    waitframes = 1;
    % Set based on work with oscilliscope
    isiTimeSecs = 1;
    % the above converted into frames
    isiTimeFrames = round(isiTimeSecs / ifi);
    % set Matlab priority to maximum
    topPriorityLevel = MaxPriority(win);
    
    %key
    escapeKey = KbName('ESCAPE');
    oneKey = KbName('1');
    expeKey = KbName('k');
 
    rng('shuffle');

%---------------------------------------------------------------------------------------%
    
    %Create our 2s tone and prepare the whitenoise parameters;
    sig_freq = 1000; 
    fs = 44100;
    t = 0:1/fs:(2-1/fs);
    sig = sin(2*pi*t*sig_freq);
    sig = (sig - min(sig))*(1-(-1))/(max(sig) - min(sig)) + (-1);
    load('whitenoise'); %(wnseq = our white noise);
    wnseq = wnseq(1:88200);
    wnseq = (wnseq - min(wnseq))*(1-(-1))/(max(wnseq) - min(wnseq)) + (-1);
    sigPower = 1;
    SNR =1;
 
%---------------------------------------------------------------------------------------%
    %Prepare the trial matrix with counterbalanced conditions and random target onsets;  
%---------------------------------------------------------------------------------------%    
    
    %loop of n iterations to make sure we give enough chances to generate at least one satisfying sequence meeting constraints;
    for iteration = 1:10000

    %Continue even if the loop gives an error;   
    try

        %Previous element: 0 at the beginning and define the first trial;
        permuts = zeros(100,1);   
        pf = 0; 
        t = randi([0,1],1,1); 

        %This takes care of the constraint in 1st column;
        while(t == pf)
          t = randi([0,1],99,1);      
        end

        %Add the first trial condition of the generated sequence; 
        permuts(1,1)= t;
        pf = t;
        cons = 2;

        %This should avoid that more than 3 consecutive elements appear in each row; 
        for c = 4:100
           t = randi([0,1],1,1);     
           cons = cons+(t == permuts(c-1,1));

               if(cons == 3)
                 t = randi([0,1],1,1);
                 while(t == permuts(c-2,1))  
                      t = randi([0,1],1,1);
                 end
                 cons = 1;
               end
           permuts(c,1) = t;

        end

        %Continue the loop even if this one is sequence does not meet the constraint (which should give an error);
        catch
            continue

    end


        %Now we have a sequence that works, but we need to make sure that the 3 conditiosn are counterbalanced (100 trials each);
        condition_0 = sum(permuts(:) == 0);
        condition_1 = sum(permuts(:) == 1);

        %If ok, we have a good pseudorandom sequence and we keep it, so we can scape and save it;
        if condition_0 == condition_1 
            break
        end

    end
    
    clear t pf i cons c condition_0 condition_1 condition_2 iteration 
 %---------------------------------------------------------------------------------------%
 
    %Now feed the trial matrix;
    alltrials(:,2) = permuts; 
    clear permuts
             
    
    for i = 1:length(alltrials)
             
        if alltrials(i,2) == 0
            alltrials(i,3) = 0;
        elseif alltrials(i,2) == 1
            alltrials(i,3) = randsample(0.3*fs:1.4*fs,1);
        end
    end
    
    
    %add 4 initials trials of examples;
    alltrials = [zeros(4,3);alltrials];
    alltrials(:,1) = 1:length(alltrials(:,2));
    alltrials(1:4,2) = [1 1 0 1]; %The first 4 trials servent as examples;
    alltrials(1:4,3) = [44100 44100 0 44100] ; %for the first 4 trials, the tone onset is fixed at 1s if present;
    alltrials(:,4) = [ones(24,1);2*ones(20,1);3*ones(20,1);4*ones(20,1);5*ones(20,1)];
    
    %Prepare the output matrix of the participant;   
    respMat = cell(length(alltrials), 9);
    respMat = array2table(respMat, 'Variable',{'PARTICIPANT','BLOCK','TRIAL_NUMB','TONE','TONE_ONSET','SNR','RESP_PART','HITS', 'REACTION_TIME'});
  

%---------------------------------------------------------------------------------------%     
                            %Now the experiment starts;
%---------------------------------------------------------------------------------------%  

for block =1:max(alltrials(:,4))

trial_block = find(alltrials(:,4) == block)'; 

for trial = trial_block  
    
    %Adapt the SNR after 2 consecutive hits or 1 miss (reqSNR = 0: only whitenoise; reqSNR = 1: tone = whitenoise); 
    
    if trial <= 4
        
        SNR_trial = 0.01; 
        hit(trial,2)= SNR_trial;
    
    elseif trial == 5
        
        SNR_trial = 0.01; 
        hit(trial,2)= SNR_trial;
    
    elseif trial >5
        
        if hit(trial-1,1)== 1 && hit(trial-2,1)== 1 && alltrials(trial-1,2)== 1 && alltrials(trial-2,2)== 1 || hit(trial-1,1)== 1 && hit(trial-2,1)== 1 && alltrials(trial-1,2)== 1 && alltrials(trial-2,2)== 0  
            SNR_trial = hit(trial-1,2)- 0.01*0.02;
                %floor SNR_trial at 5% minimum;
                if SNR_trial <= 0.00001
                    SNR_trial = 0.00001;
                end  
            hit(trial,2)= SNR_trial;
                        
        elseif hit(trial-1,1)== 1 && hit(trial-2,1)== 1 && alltrials(trial-1,2)== 0 && alltrials(trial-2,2)== 1 || hit(trial-1,1)== 1 && hit(trial-2,1)== 1 && alltrials(trial-1,2)== 0 && alltrials(trial-2,2)== 0 || hit(trial-1,1) == 1 && hit(trial-2,1) == 0
            SNR_trial = hit(trial-1,2);
            hit(trial,2)= SNR_trial;                       
            
        elseif hit(trial-1,1) == 0 
            SNR_trial = hit(trial-1,2)+ 0.01*0.02;
                %floor SNR_trial at 95% maximum;
                if  SNR_trial >= 0.01   
                    SNR_trial = 0.01;
                end 
            hit(trial,2)= SNR_trial;
        end
        
    end
    
    %Generate the signal + whitenoise at each trial;
    noisePower = sigPower/SNR_trial;
    y = sqrt(noisePower)*wnseq;
    tone_temp = sig + y; 
    tone_temp = (tone_temp - min(tone_temp))*(1-(-1))/(max(tone_temp) - min(tone_temp)) + (-1);
    
    %Now create the new stimulus at each trial;
    if alltrials(trial,2) == 0      
        %No tone in the trial;
        stimulus = [wnseq'];
        stimulus(:,2) = stimulus(:,1);
      
    elseif alltrials(trial,2) == 1
        
        %Tone with the appropriate onset;
        n1 = wnseq(1:alltrials(trial,3)-1);
        tone_temp2 = tone_temp(alltrials(trial,3):alltrials(trial,3)+4410);
        n2 = wnseq(alltrials(trial,3)+4410+1:end);       
        stimulus = [n1';tone_temp2';n2'];
        stimulus(:,2) = stimulus(:,1);

    end
    
        %Buffer the stimulus, ready for the trial;
        PsychPortAudio('FillBuffer', pahandle, stimulus');
        

%---------------------------------------------------------------------------------------%    
    
    % If this is the first trial we present the instruction screen and wait for a key-press
    if block == 1 && trial == 1
        DrawFormattedText(win, 'Thank you for your participation: In this first task, your will have to detect the auditory targets (tones).\n\n\n\n\n\nEvery time you detect a tone, press 1 as fast as possible.\n\n\n\nIf you dont hear any tone, wait for the next trial without responding.\n\n\n\nPlease focus on the fixation cross during the experiment.\n\n\n\n\n\nWhen ready, press space bar to start with some examples',...
            'center', 'center', black);
        Screen('Flip', win);
        pause(1)
        respToBeMade = true; % set it so that a response is required
        while respToBeMade==true
            [keyIsDown,secs, keyCode] = KbCheck;
            if keyCode(escapeKey) % if ESC is pressed, terminate the experiment
                ShowCursor;
                sca;
                return
            elseif keyIsDown && ~keyCode(escapeKey)
                respToBeMade = false;
            end
        end
    end     

       
%-------------------------------------------------------------------------------        
    %Trial structure;

    %draw the fixation cross for 500 - 1000ms;
    rframes = randi([100,200]);
    for frame = 1:rframes
    %Draw the fixation point
    DrawFormattedText(win, '+', 'center', 'center', black);
    Screen('Flip', win);
    end   

%     %Send the trigger;
%     lab_put_code(L,1);

    %draw a red fixation cross during the stimulus playing(2s);
    DrawFormattedText(win, '+', 'center', 'center', [255 0 0]); %[255,0,0]
    Screen('Flip', win);          
    %Play the stimulus;
    PsychPortAudio('Start', pahandle, 1, 0, 1);

    %Start the timer for calculating the response times:onset of the tone;
    tStart = GetSecs;
    
    %%%%RESPONSE PARTICIPANT;
    %Check the keyboard for responses
    respToBeMade = true; % set it so that a response is required
    while respToBeMade==true
        [keyIsDown,secs, keyCode] = KbCheck;
        if keyCode(escapeKey) % if ESC is pressed, terminate the experiment
            ShowCursor;
            sca;
            return
        elseif keyCode(oneKey)
            response_trial = 1;               
            respToBeMade = false;
            rt_trial = GetSecs - tStart; 

                if response_trial == 1 && alltrials(trial,2) == 1 && rt_trial > alltrials(trial,3)/fs
                    hit(trial,1)= 1;
                elseif response_trial == 1 && alltrials(trial,2) == 1 && rt_trial < alltrials(trial,3)/fs
                    hit(trial,1)= 0;                        
                elseif response_trial == 1 && alltrials(trial,2) == 0
                    hit(trial,1)= 0;
                end

        elseif GetSecs - tStart > 2 %if no response is given after 2 s
            response_trial = 0;
            respToBeMade = false;
            rt_trial = GetSecs - tStart;

                if response_trial == 0 && alltrials(trial,2) == 0
                    hit(trial,1)= 1;
                elseif response_trial == 0 && alltrials(trial,2) == 1
                    hit(trial,1)= 0;
                end
        end
    end
    
    %Stop the sound and close the trial;
    WaitSecs(0.1);
    PsychPortAudio('Stop', pahandle, 1); 
    Screen('Close');

    %Feed participants data after each trial;      
    respMat.PARTICIPANT{trial,1} = participant_num;
    respMat.BLOCK{trial,1} = block;
    respMat.TRIAL_NUMB{trial,1} = trial;
    respMat.TONE{trial,1} = alltrials(trial,2);
    respMat.TONE_ONSET{trial,1} = alltrials(trial,3);
    respMat.SNR{trial,1}= hit(trial,2);
    respMat.RESP_PART{trial,1} = response_trial;
    respMat.HITS{trial,1} = hit(trial,1);
    respMat.REACTION_TIME{trial,1} = rt_trial;

    %End of examples slide;
    if trial == 4 && block ==1

        DrawFormattedText(win, 'End of the examples.\n\nIf you have no further questions, press space bar to continue the experiment.',...
            'center', 'center', black);
        Screen('Flip', win);
        pause(1)
        respToBeMade = true; % set it so that a response is required
        while respToBeMade==true
            [keyIsDown,secs, keyCode] = KbCheck;
            if keyCode(escapeKey) % if ESC is pressed, terminate the experiment
                ShowCursor;
                sca;
                return
            elseif keyIsDown && ~keyCode(escapeKey)
                respToBeMade = false;
            end
        end
    end

    %Break slide;  
    if trial == max(trial_block) && block < max(alltrials(:,4))

        DrawFormattedText(win, 'Take a break now.\n\n\n\nPress space bar to continue the experiment when you are ready.',...
            'center', 'center', black);
        Screen('Flip', win);
        pause(1)
        respToBeMade = true; % set it so that a response is required
        while respToBeMade==true
            [keyIsDown,secs, keyCode] = KbCheck;
            if keyCode(escapeKey) % if ESC is pressed, terminate the experiment
                ShowCursor;
                sca;
                return
            elseif keyIsDown && ~keyCode(escapeKey)
                respToBeMade = false;
            end
        end

    %End the experiment;
    elseif trial == max(trial_block) && block == max(alltrials(:,4))

        DrawFormattedText(win, 'Well done! The first task is finished.\n\n\n\nWait for the experimenter to start the second experiment',...
            'center', 'center', black);
        Screen('Flip', win);
        WaitSecs(1);
        %Key 'K' has to be pressed to continue --- CAREFULL !!!!----
        respToBeMade = true;
        while respToBeMade==true
             [keyIsDown,secs, keyCode] = KbCheck;
             if keyCode(escapeKey) % if ESC is pressed, terminate the experiment;
                 ShowCursor;
                 sca;
                 return
             elseif keyCode(expeKey) % if 'K' is pressed, terminate the experiment;
                 respToBeMade = false;
             else
                 respToBeMade = true;
             end
        end
    end

end


end


    %Participant outputs to configure the following task parameter;participant_parameters = [];
    participant_parameters(1,1) = mean(cell2mat(respMat.REACTION_TIME(cell2mat(respMat.HITS)==1 & cell2mat(respMat.TONE)==1)) - cell2mat(respMat.TONE_ONSET(cell2mat(respMat.HITS)==1 & cell2mat(respMat.TONE)==1))/fs); 
    participant_parameters(1,2) = std(cell2mat(respMat.REACTION_TIME(cell2mat(respMat.HITS)==1 & cell2mat(respMat.TONE)==1)) - cell2mat(respMat.TONE_ONSET(cell2mat(respMat.HITS)==1 & cell2mat(respMat.TONE)==1))/fs); 
    participant_parameters(2,1) = mean(cell2mat(respMat.SNR(50:104)));
    participant_parameters(2,2) = std(cell2mat(respMat.SNR(50:104)));
    save(['subj', num2str(participant_num), '_parameters'], 'participant_parameters');

    %Export the participant logfile;
    cd XXX\; 
    mkdir(['subj' num2str(participant_num)]);
    cd(['XXX\subj' num2str(participant_num)]);
    save(['subj' num2str(participant_num) '-Calibration_task'], 'respMat','participant_parameters','alltrials'); 
     

catch
end
   
%     %Close LabJack2
%     lab_close; 

    %Same close sequence and psychlasterror
    cd XXX\;
    sca
    ShowCursor;
    Priority(0);

        
%% 