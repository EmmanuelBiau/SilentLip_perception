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
%ASIO4ALL v2. May need to change if new hardware is installed
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


% HideCursor;
    
%----------------------------------------------------------------------------------------------------------%
%   Set up the display and timing
%----------------------------------------------------------------------------------------------------------%
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
    
    %-----------------------------------------------------%
    %Define rectangle for frame counting with photodiode;
    baseRect = [0 0 50 50];
    %Center the rectangle on the centre of the screen using fractional pixel values
    Rect = CenterRectOnPointd(baseRect, xCenter-615, yCenter);
    %Set the color of our square to white.
    rectColor = [1 1 1];
    %-----------------------------------------------------%

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
    
    %Key
    escapeKey = KbName('ESCAPE');
    oneKey = KbName('1');
    expeKey = KbName('k');

    rng('shuffle');
    
%---------------------------------------------------------------------------------------%
    
    %Create our 2s tone and prepare the whitenoise parameters;
    sig_freq = 1000; 
    fs = 44100;
    t = 0:1/fs:(5-1/fs);
    sig = sin(2*pi*t*sig_freq);
    sig = (sig - min(sig))*(1-(-1))/(max(sig) - min(sig)) + (-1);%normalize(sig, 'range',[-1 1]);
    load('whitenoise'); %(wnseq = our white noise);
    wnseq = (wnseq - min(wnseq))*(1-(-1))/(max(wnseq) - min(wnseq)) + (-1);
    sigPower = 1;
    
    %load participant's and the SNR determined from calibration task;
    load(['subj', num2str(participant_num), '_parameters']);

%---------------------------------------------------------------------------------------%
    %Prepare the trial matrix with counterbalanced conditions and random target onsets;  
%---------------------------------------------------------------------------------------%    
    
    %loop of n iterations to make sure we give enough chances to generate at least one satisfying sequence meeting constraints;
    for iteration = 1:10000

    %Continue even if the loop gives an error;   
    try

        %Previous element: 0 at the beginning and define the first trial;
        permuts = zeros(300,1);   
        pf = 0; 
        t = randi([0,2],1,1); 

        %This takes care of the constraint in 1st column;
        while(t == pf)
          t = randi([0,2],299,1);      
        end

        %Add the first trial condition of the generated sequence; 
        permuts(1,1)= t;
        pf = t;
        cons = 2;

        %This should avoid that more than 3 consecutive elements appear in each row; 
        for c = 4:300
           t = randi([0,2],1,1);     
           cons = cons+(t == permuts(c-1,1));

               if(cons == 3)
                 t = randi([0,2],1,1);
                 while(t == permuts(c-2,1))  
                      t = randi([0,2],1,1);
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
        condition_2 = sum(permuts(:) == 2);

        %If ok, we have a good pseudorandom sequence and we keep it, so we can scape and save it;
        if condition_0 == condition_1 && condition_0 == condition_2
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
           alltrials(i,4) = 0; 
           
        elseif alltrials(i,2) == 1
            alltrials(i,3) = randsample(0.3*fs:4.3*fs,1);
            alltrials(i,4) = 0;
            
        elseif alltrials(i,2) == 2
            alltrials(i,3) = randsample(0.3*fs:3*fs,1);
            alltrials(i,4) = randsample((alltrials(i,3)+1*fs):4.5*fs,1);
        end
    end  
    
    
    %add 4 initials trials of examples;
    alltrials = [zeros(4,4);alltrials];
    alltrials(:,1) = 1:length(alltrials(:,2));
    alltrials(1:4,2) = [1 2 0 1]; %The first 4 trials servent as examples;
    alltrials(1:4,3) = [44100 44100 0 2.5*44100] ; %for the first 4 trials, the tone onset is fixed at 1s if present;
    alltrials(1:4,4) = [0 3*44100 0 0] ; 
    alltrials(:,5) = [ones(54,1);2*ones(50,1);3*ones(50,1);4*ones(50,1);5*ones(50,1);6*ones(50,1)];
    
    %load video list and randomize order;
    video_list_temp = dataset('File', 'List_movies.csv', 'Delimiter' ,',');
    video_list_temp2 = video_list_temp(5:end,1:2);
    video_list = [video_list_temp(1:4,:);video_list_temp2(randperm(size(video_list_temp2, 1)),:);video_list_temp2(randperm(size(video_list_temp2, 1)),:);video_list_temp2(randperm(size(video_list_temp2, 1)),:);video_list_temp2(randperm(size(video_list_temp2, 1)),:);video_list_temp2(randperm(size(video_list_temp2, 1)),:)];
    clear video_list_temp video_list_temp2
    
    %Prepare the output matrix of the participant;   
    respMat = cell(length(alltrials), 13);
    respMat = array2table(respMat, 'Variable',{'PARTICIPANT','BLOCK','TRIAL_NUMB','VIDEO','FREQPEAK','SNR','TONE','TONE_ONSET1','TONE_ONSET2','KEYPRESS', 'REACTION_TIMES', 'SCREEN_ONSET', 'SCREEN_OFFSET'});
  

%---------------------------------------------------------------------------------------%     
                            %Now the experiment starts;
%---------------------------------------------------------------------------------------%  


for block =1:max(alltrials(:,5))

trial_block = find(alltrials(:,5) == block)'; 

for trial = trial_block  
    
    %Load and prepare the video of the trial;
    cd XXX\; 
    avi_object= VideoReader(video_list.moviename{trial,1});
    % read in all frames of avi_object. ***** Display resolution: 1280 x1024 with 75 Hz refreshing rate***** 
    avi_frames2= read(avi_object);
    avi_frames=avi_frames2(2:719,2:1279,:,:); % resize the frames from the movie to avoid black edges;
    % get number of frames in avi_object
    num_frames=get(avi_object, 'NumberofFrames');
    % get framerate of avi_object
    framerate=get(avi_object, 'FrameRate');
    
  
    %Prepare the auditory stimulus of the trial;
    %Set the SNR a bit louder for the first 4 example trials and then, at the SNR determined by the calibration task; 
    if trial <= 4 && block ==1
        SNR_trial = 0.01;         
    
    elseif trial > 4
        SNR_trial = participant_parameters(2,1);
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
         
    elseif alltrials(trial,2) == 2
        
        %Tone with the appropriate onset;
        n1 = wnseq(1:alltrials(trial,3)-1);
        tone_temp2 = tone_temp(alltrials(trial,3):alltrials(trial,3)+4410);
        n2 = wnseq(alltrials(trial,3)+4410+1:alltrials(trial,4)-1);
        tone_temp3 = tone_temp(alltrials(trial,4):alltrials(trial,4)+4410);
        n3 = wnseq(alltrials(trial,4)+4410+1:end);
        
        stimulus = [n1';tone_temp2';n2';tone_temp3';n3'];
        stimulus(:,2) = stimulus(:,1);      
        
    end
    
        %Buffer the stimulus, ready for the trial;
        PsychPortAudio('FillBuffer', pahandle, stimulus');
        

%---------------------------------------------------------------------------------------%    
    
    % If this is the first trial we present the instruction screen and wait for a key-press
    if block == 1 && trial == 1
        DrawFormattedText(win, 'SECOND TASK: You will have to detect the same audio tones.\n\nBut now, there will be silent videos displayed at the same time.\n\n\n\n\n\nFor each trial, there can be only 0, 1 or 2 Tones.\n\n\n\nEvery time you detect a tone, press 1 as fast as possible.\n\nIf you do not hear any tone, wait for the next trial without responding.\n\n\n\n\nPlease focus on the fixation cross during the experiment.\n\n\n\n\n\n\n\n(press space bar when you are ready to start with some examples)',...
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
      
%---------------------------------------------------------------------------------------%        
    
    %draw the fixation cross for 500 - 1250ms;
    rframes = randi([100,250]);
    for frame = 1:rframes
        DrawFormattedText(win, '+', 'center', 'center', [250 0 0]);
        Screen('Flip', win);
    end
    
    %Make the video by grabbing the frames;
    for frame = 1:num_frames
        frameToShow{frame} = Screen('MakeTexture', win, avi_frames(:,:,:,frame));
    end
    
    %Screen with 75 hz refresh rate and videos at 25 fps: We need 3 displays of each
    %video frame for a 3-second video presentation. 
    flipsPerFrame = 3;
    PsychPortAudio('Start', pahandle, 1, inf, 0);
    WaitSecs(0.1);   
    
    %Initiate the Key Response functions for every trial;
    KbQueueCreate;
    KbQueueStart; 
    
    %sync to the vertical retrace and pause execution of the script until the Flip has happened;
    [refreshTimestamp, StimulusOnsetTime] = Screen('Flip', win);

    %play the sound after 4 ifi after stimulus onset; 
    PsychPortAudio('RescheduleStart', pahandle, StimulusOnsetTime+(ifi*4));

%-------------------------------------------------------------------------------------------------------------%
    % Here little trick to maximize the control about the exact timing between video and sound:
    % (1) minimize buffering effect when building the video from the frames;
    % (2) display a flicker on the left of the screen to measure the onset of the video with the photodiode;
%--------------------------------------------------------------------------------------------------------------%    

    %draw the white square and the first frame into buffer;
    ifrm = 1; 
     Screen('DrawTexture', win, frameToShow{1}, [], []);
     Screen('FillRect', win, rectColor, Rect);
   
%     %STIMULI STARTED TRIGGER;
%     lab_put_code(L,2);      

    %Flip to show our first frame and the white square;
    [refreshTimestamp] = Screen('Flip', win, refreshTimestamp + (4-0.5)*ifi); 

    %Prepare the response detection config;
    count = 0;
    key_press_count = zeros(1,1); %store the number of key pressed every trials; 
    key_press_rt = zeros(1,1); %store the onsets of key pressed every trials; 
    
    %Get the onset of the stimulus;
    tStart = GetSecs;
    
    %Draw the movie;
    %ifrm = 0; % the current frame (starting at 0);
    for frame = 2:num_frames % for each frame in the movie;
        ifrm = ifrm+1; % the current frame;
        for iflip = 1:flipsPerFrame
            Screen('DrawTexture', win, frameToShow{frame}, [], []);
            Screen('FillRect', win, rectColor, Rect);
            DrawFormattedText(win, '+', 'center', 'center', black);
            % ('Flip', 'window', 'when'). When = 0.5 * ifi sec after last
            % refresh timestamp. ie, flip at the first possible
            % retrace after the last flip timestamp plus the
            % duration in sec of half of an inter-frame interval             
            [refreshTimestamp] = Screen('Flip', win, refreshTimestamp + 0.5*ifi);
            [ pressed, firstPress]= KbQueueCheck;
            if firstPress(oneKey)
                count = count + 1;
                key_press_count(count,1) = 1;
                key_press_rt(count,1) = firstPress(oneKey)- tStart;
            elseif frame == num_frames && iflip == flipsPerFrame && count == 0
                key_press_count(1,1) = 0;
                key_press_rt(1,1) = GetSecs - tStart;
            end
            
            %clear the firsPress for the next participant's key press;
            KbQueueFlush;
            
        end
    end
    
    %Stops the key information storage;
    KbQueueRelease;
    PsychPortAudio('Stop', pahandle, 1); 
    
    %Get the offset of the stimulus;
    tEnd = GetSecs-tStart;
    
    %Clear the textures created by reading the frames of the movie
    for frame = 1:num_frames
        Screen('Close', frameToShow{frame});
    end

                   
    %Feed participants data after each trial;      
    respMat.PARTICIPANT{trial,1} = participant_num;
    respMat.BLOCK{trial,1} = block;
    respMat.TRIAL_NUMB{trial,1} = trial;
    respMat.VIDEO{trial,1} = video_list.moviename{trial,1};
    respMat.FREQPEAK{trial,1} = video_list.freqpeak(trial,1);
    respMat.TONE{trial,1} = alltrials(trial,2);
    respMat.TONE_ONSET1{trial,1} = alltrials(trial,3);
    respMat.TONE_ONSET2{trial,1} = alltrials(trial,4);
    respMat.SNR{trial,1} = SNR_trial;
    respMat.KEYPRESS{trial,1} = key_press_count;
    respMat.REACTION_TIMES{trial} = key_press_rt;
    respMat.SCREEN_ONSET{trial,1} = tStart;
    respMat.SCREEN_OFFSET{trial,1} = tEnd;
    
    clear key_press_count key_press_rt

%---------------------------------------------------------------------------------------%    
    %End of examples slide;
    if trial == 4 && block ==1

        DrawFormattedText(win, 'End of the examples.\n\nIf you have no further question, press space bar to continue the experiment.',...
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
        
    %Break slide;
    elseif trial == max(trial_block) && block < max(alltrials(:,5))

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
    elseif trial == max(trial_block) && block == max(alltrials(:,5))

        DrawFormattedText(win, 'Well done! This task is finished.\n\n\n\nWait for the experimenter to start the second experiment',...
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

    %Export the participant logfile;
    cd(['XXX\subj' num2str(participant_num)]);
    save(['subj' num2str(participant_num) '-Tone_detectionTask'],'respMat','participant_parameters','alltrials');

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
