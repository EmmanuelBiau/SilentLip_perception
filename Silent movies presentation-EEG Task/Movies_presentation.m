%Adapt paths - Psychotoolbox - triggersender to your machine;
close all; clear all; clc; 
cd XXX\; addpath XXX\;
addpath('XXX\Psychtoolbox\MatlabWindowsFilesR2007a\');
addpath(genpath('XXX\Psychtoolbox\'));
AssertOpenGL;

% Initialise LabJack;
lab_init; 

% Setup PTB with some default values;
PsychDefaultSetup(2);
InitializePsychSound(1);

%----------------------------------------------------------------------------------------------------------
%                                       Set up the audio/ports
%----------------------------------------------------------------------------------------------------------
% Request latency mode 2
reqlatencyclass = 2;
% Requested output frequency. Must set this. 96000, 48000, or 44100. Ours are 44100
freq = 44100;
buffersize = 512; %this is not used later?
latbias = []; 
% Open audio device for low-latency output
pahandle = PsychPortAudio('Open', [], [], reqlatencyclass, freq, 1);

% Tell driver about hardware's inherent latency, determined via calibration once
prelat = PsychPortAudio('LatencyBias', pahandle, latbias);
postlat = PsychPortAudio('LatencyBias', pahandle); %not sure if a latbias should be added 29/11/2017

% Generate a beep at 1000Hz, 0.1 sec, 50% amplitude, to avoid latency later on;
mynoise(1,:) = 0.5*MakeBeep(1000, 0.1, freq);
% Fill buffer with audio data;
PsychPortAudio('FillBuffer', pahandle, mynoise);
PsychPortAudio('Start', pahandle, 1, 0, 1);
PsychPortAudio('Stop', pahandle, 1);

%----------------------------------------------------------------------------------------------------------%
%                                            START THE EXPERIMENT NOW 
%----------------------------------------------------------------------------------------------------------%

% Create dialog box to gather participant number: 
participant_num = inputdlg('Participant Number', 'Participant Number');
participant_num = str2double(participant_num); 


try
    % Comment this out for testing, uncomment before real participants;
    HideCursor;
    
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
    
    %----------------------------------------------------------------------------------------------------------
    %   Display and timing setup complete
    %----------------------------------------------------------------------------------------------------------
    
    %----------------------------------------------------------------------------------------------------------
    %   Keyboard information: Define the keyboard keys that are listened
    %   for. We will be using the QWERTY keys 1-5 as response keys and the
    %   escape key as an exit/reset key
    %----------------------------------------------------------------------------------------------------------
    escapeKey = KbName('ESCAPE');
    oneKey = KbName('1!');
    twoKey = KbName('2@');
    threeKey = KbName('3#');
    fourKey = KbName('4$');
    fiveKey = KbName('5%');

    % Shuffle trial order;
    rng('shuffle');
    

%%   
    %----------------------------------------------------------------------------------------------------------
    %   3.0 MOVIE ONLY BLOCK;
    %----------------------------------------------------------------------------------------------------------
    
    % load Movie only files;
    cd XXX\;
    Movies = dataset('File', 'Movies_list.csv', 'Delimiter' ,',');
    
    numTrials = 60;

    % Define the response matrix for the sync test trials.
    respMat = cell(numTrials, 9);
    respMat = array2table(respMat, 'Variable',{'PARTICIPANT','BLOCK','TRIAL','STIM_NUMB','MOVIENAME','FREQPEAK','PARTICIP_RESP','REACTION_TIME', 'RATING_ONSET'});
    
    % Shuffle trial order
    shuffler = Shuffle(1:numTrials);
    
    %----------------------------------------------------------------------------------------------------------
    %   2.0 EXPERIMENTAL LOOP - MOVIE TEST TRIALS
    %-----------------------------------------------------------------------------------------------------------
    
    % Begin the trial loop for the SyncTest task.
    for trial = 1:numTrials
        
        %------------------------------------------------------------------------------------------------------
        %   2.0 Read in and prepare video - see comments in the encoding section
        %------------------------------------------------------------------------------------------------------
        cd XXX\; 
        avi_object=VideoReader(Movies.MovieFile{shuffler(trial),1});
        avi_frames2= read(avi_object);
        avi_frames=avi_frames2(2:719,2:1279,:,:); % resize the frames from the movie to avoid black edges;
        num_frames=get(avi_object, 'NumberofFrames');
        framerate=get(avi_object, 'FrameRate');
        
        %------------------------------------------------------------------------------------------------------
        %   2.0 BEGIN MOVIE TRIALS
        %------------------------------------------------------------------------------------------------------
        
        % If this is the first trial we present the instruction screen and wait for a key-press;
        if trial == 1

            % Instructions for the sound only block;
            DrawFormattedText(win, 'Finally, you will be presented only with video clips. \n\nFor each video, your task is to rate how emotional you think the video is. \n\nPress Spacebar To Begin.',...
            'center', 'center', 0);
            Screen('Flip', win);
            pause(1)
            KbStrokeWait;
            
        end
        
        %Now we present the isi interval with fixation point;
        rframes = randi([100 250]);
        for frame = 1:isiTimeFrames+rframes
            % Draw the fixation point;
            DrawFormattedText(win, '+', 'center', 'center', 0);
            Screen('Flip', win);
        end
        
        %------------------------------------------------------------------------------------------------------
        %   2.0 START STIMULUS DISPLAY LOOP FOR THE CURRENT BLOCK
        %------------------------------------------------------------------------------------------------------
        
        % Make the video by grabbing the frames;
        for frame = 1:num_frames
            frameToShow{frame} = Screen('MakeTexture', win, avi_frames(:,:,:,frame));
        end

        % Use a screen with 75 hz refresh rate. The videos have a 25 framerate. So we need 3 displays of each
        % video frame for a 3-second video presentation; 
        flipsPerFrame = 3;
        
        % VIDEO TEST STIMULI ON TRIGGER
        lab_put_code(L,4); 
 
        % Draw the movie;
        [refreshTimestamp] = Screen('Flip', win); 
        ifrm = 0; % the current frame (starting at 0)
        for frame = 1:num_frames % for each frame in the movie
            ifrm = ifrm+1; % the current frame
            for iflip = 1:flipsPerFrame
                Screen('DrawTexture', win, frameToShow{frame}, [], []);
                DrawFormattedText(win, '+', 'center', 'center', 0);
                [refreshTimestamp] = Screen('Flip', win, refreshTimestamp + 0.5*ifi);
            end
        end
        
        % Clear the textures created by reading the frames of the movie;
        for frame = 1:num_frames
            Screen('Close', frameToShow{frame});
        end
        
        % Draw the rating screen;
        DrawFormattedText(win, 'Using the numbers 1 through 5, rate how emotional you think the sound was.\n\n 1 means that the sound was completely neutral, \n\n and 5 means that the sound was very emotional',...
        'center', 'center', 0);
        Screen('Flip', win);
        
        % Start the timer for calculating the response times;
        tStart = GetSecs;
        
        % Check the keyboard for responses;
        respToBeMade = true; % set it so that a response is required;
        while respToBeMade==true
            [keyIsDown,secs, keyCode] = KbCheck;
            if keyCode(escapeKey) % if ESC is pressed, terminate the experiment;
                ShowCursor;
                sca;
                return
            elseif keyCode(oneKey)
                response = 1;
                respToBeMade = false;
            elseif keyCode(twoKey)
                response = 2;
                respToBeMade = false;
            elseif keyCode(threeKey)
                response = 3;
                respToBeMade = false;
            elseif keyCode(fourKey)
                response = 4;
                respToBeMade = false;
            elseif keyCode(fiveKey)
                response = 5;
                respToBeMade = false;
            end
        end
        
        % Get the time when the response is made, and calculate the RT;
        tEnd = GetSecs;
        responsert = tEnd - tStart;
       
        % Record trial data to the data matrix;
        respMat.PARTICIPANT{trial,1} = participant_num;
        respMat.BLOCK{trial,1} = 'video'; % block name;
        respMat.TRIAL{trial,1} = 300+trial; 
        respMat.STIM_NUMB{trial,1} = Movies.stim_number(shuffler(trial),1); % the original stimulus number;
        respMat.MOVIENAME{trial,1} = Movies.MovieFile{shuffler(trial),1}; % video name;
        respMat.FREQPEAK{trial,1} = Movies.freqpeak(shuffler(trial),1); % frequency of peak of stimulus;
        respMat.PARTICIP_RESP{trial,1} = response; % the participant response key;
        respMat.REACTION_TIME{trial,1} = responsert; % the participant RT;
        respMat.RATING_ONSET{trial,1} = tStart; % the time of the presentation of the rating screen;
               
    end

    %Save the Ratings data to a file;
    cd XXX\; 
    mkdir(['subj' num2str(participant_num)]);
    cd(['XXX\subj_' num2str(participant_num)]);
    writetable(respMat,['subj' num2str(participant_num) '-Movie.csv'],'Delimiter', ',','QuoteStrings',true);
        
    %------------------------------------------------------------------------------------------------------
    %   2.0 END STIMULUS DISPLAY LOOP FOR THE CURRENT BLOCK
    %   THIS ENDS 2.0
    %------------------------------------------------------------------------------------------------------
  
    DrawFormattedText(win, 'You are now finished! \n\nThank you for your participation. \n\nPress Spacebar to end the experiment.',...
    'center', 'center', 0);
    Screen('Flip', win);
    pause(1)
    KbStrokeWait;
    sca
    
catch  
    ShowCursor;
    Priority(0);

    %Close LabJack2
    lab_close

end
sca

%%