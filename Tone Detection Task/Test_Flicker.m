close all; clear all; clc; cd C:\experiments\Manu\BEHAVIOR_Experiment\;
AssertOpenGL; GetSecs; WaitSecs(0.1);
%Setup PTB with some default values
PsychDefaultSetup(2);

%Initialise LabJack
lab_init; 

%Key
escapeKey = KbName('ESCAPE');

%-----------------------------------------------------%
%               Background Properties
%-----------------------------------------------------%

screenid = max(Screen('Screens'));
%create a window
white = WhiteIndex(screenid);
black = BlackIndex(screenid);
[win, winRect] = PsychImaging('OpenWindow', screenid, 0.5);
%define size of x and y based on screen size
[screenXpixels, screenYpixels] = Screen('WindowSize', win);
%define screen center based on screen size
[xCenter, yCenter] = RectCenter(winRect);
%set test font size
Screen('TextSize', win, 20);
%set text font
Screen('TextFont', win, 'Helvetica');

%-----------------------------------------------------%
%               Flicker Properties
%-----------------------------------------------------%

%Define rectangle for frame counting with photodiode;
baseRect = [0 0 50 50];
%Center the rectangle on the centre of the screen using fractional pixel values
Rect = CenterRectOnPointd(baseRect, xCenter-615, yCenter);
%Set the color of our square to white.
rectColor = [1 1 1];

%-----------------------------------------------------%

%Send the trigger;
lab_put_code(L,2);

%Now Draw the grey background with the flicker on the left;
DrawFormattedText(win, '+', 'center', 'center', black); 
Screen('FillRect', win, rectColor, Rect);
Screen('Flip', win);
pause(1)
respToBeMade = true;
while respToBeMade==true
    [keyIsDown,secs, keyCode] = KbCheck;
    if keyCode(escapeKey)
        ShowCursor;
        sca;
        return
    end
end

%Close LabJack2
lab_close; 
 clear all;

%%
