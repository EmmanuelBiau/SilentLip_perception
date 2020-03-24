%Clen the logfile: Import the participant logfile and detect Hits,Miss, Reaction times and trial infos;
clearvars;clc;
fs = 44100;

subjects = 1:29;

for s = subjects

%load the raw logfile after the experiment;
cd(['XXX\subj' num2str(s)]);
load(['subj' num2str(s) '-Tone_detectionTask_lags']);

%Prepare the output logfile of the participant;   
Logfile = cell(length(alltrials), 18);
Logfile = array2table(Logfile, 'Variable',{'PARTICIPANT','BLOCK','TRIAL_NUMB','VIDEO','FREQPEAK','SNR','TONE','TONE_ONSET1','TONE_ONSET2','RESP_PART1','RESP_PART2','HITS_TONE1','HITS_TONE2','REACTION_TIME1','REACTION_TIME2', 'SCREEN_ONSET', 'SCREEN_OFFSET','VA_LAG'});


%---------------------- HITS/MISS STATEMENTS ------------------------------%

%1- DETECT THE KEY RESPONSES (RTs) IN THE CORRECT ZONE FOR TONE1 and TONE2; 

    for trial = 1:length(respMat.TRIAL_NUMB)

        count = sum(respMat.KEYPRESS{trial});
        key_press_rt = respMat.REACTION_TIMES{trial};
        key_press_count = respMat.KEYPRESS{trial};

        switch respMat.TONE{trial}

            case 0

                    if count == 0
                       rt_tone1_hit = 1;
                       rt_tone2_hit = 1;

                    elseif count ~= 0
                       rt_tone1_miss = key_press_rt(find(key_press_rt > 0.1 & key_press_rt < 4.95 ,1, 'first'));
                       rt_tone2_miss = 0;
                       clear rt_tone1_hit 
                       clear rt_tone2_hit 
                    end

            case 1   

                    %Hit 1 Tone with only one keypress;
                    if length(key_press_rt) == 1 && key_press_rt(1,1) > alltrials(trial,3)/fs + 0.1 && key_press_rt(1,1) < 4.95
                        rt_tone1_hit = key_press_rt(1,1);

                    %Hit 1 Tone more than one keypress;
                    elseif length(key_press_rt) > 1 && length(key_press_rt) <= 2
                        rt_tone1_hit = key_press_rt(find(key_press_rt > alltrials(trial,3)/fs + 0.1 & key_press_rt < 4.95 , 1, 'first'));

                             if isempty(rt_tone1_hit)
                                rt_tone1_miss = key_press_rt(find(key_press_rt > alltrials(trial,3)/fs + 0.1 & key_press_rt < 4.95, 1, 'first'));
                                clear rt_tone1_hit   

                                    if isempty(rt_tone1_miss)
                                        rt_tone1_miss = 0;
                                        clear rt_tone1_hit
                                    end
                             end   

                    %miss 1 Tone because keypress too slow;    
                    elseif length(key_press_rt) <= 2 && key_press_rt(1,1) > 4.95
                        rt_tone1_miss = key_press_rt(find(key_press_rt > 4.95 , 1, 'first'));
                        clear rt_tone1_hit

                    %miss 1 Tone because too many keypress or no keypress at all;               
                    elseif length(key_press_rt) > 2 || count == 0
                        rt_tone1_miss = 0;
                        clear rt_tone1_hit                                  
                    end


             case 2     
                    
                    %Hit 2 Tones with 2 or 3 keypress;
                    if length(key_press_rt) >= 1 && length(key_press_rt) <= 3
                        rt_tone1_hit = key_press_rt(find(key_press_rt > alltrials(trial,3)/fs + 0.1 & key_press_rt < alltrials(trial,4)/fs, 1, 'first'));
                        rt_tone2_hit = key_press_rt(find(key_press_rt > alltrials(trial,4)/fs + 0.1 & key_press_rt < 4.95 , 1, 'first'));

                            if isempty(rt_tone1_hit)
                                rt_tone1_miss = key_press_rt(find(key_press_rt > alltrials(trial,3)/fs + 0.1 & key_press_rt < alltrials(trial,4)/fs, 1, 'first'));
                                clear rt_tone1_hit

                                    if isempty(rt_tone1_miss)
                                        rt_tone1_miss = 0;
                                        clear rt_tone1_hit
                                    end
                            end

                            if isempty(rt_tone2_hit)   
                                rt_tone2_miss = key_press_rt(find(key_press_rt > alltrials(trial,4)/fs + 0.1 & key_press_rt < 4.95 ,1, 'first'));
                                clear rt_tone2_hit

                                    if isempty(rt_tone2_miss)
                                        rt_tone2_miss = 0;
                                        clear rt_tone2_hit
                                    end
                            end
                            
                    %Too many keypresses;
                    elseif length(key_press_rt) > 3
                        rt_tone1_miss = 0;
                        rt_tone2_miss = 0;
                        clear rt_tone1_hit 
                        clear rt_tone2_hit                   
                    end                   
        end

        
        %2- FEED PARTICIPANT's respMAT OF THE CORRESPONDING TRIAL: IF the RTs are ok for their corresponding Tone onsets, then they are hits (1);

        switch alltrials(trial,2)

            case 0
                    %if hit for no Tone trials;
                    if exist('rt_tone1_hit','var') && exist('rt_tone2_hit','var')
                        Logfile.RESP_PART1{trial,1} = NaN;
                        Logfile.RESP_PART2{trial,1} = NaN;
                        Logfile.HITS_TONE1{trial,1} = 1;
                        Logfile.HITS_TONE2{trial,1} = 1;
                        Logfile.REACTION_TIME1{trial,1} = NaN;
                        Logfile.REACTION_TIME2{trial,1} = NaN;   

                    %if FA for no Tone trials;
                    elseif ~exist('rt_tone1_hit','var') && ~exist('rt_tone2_hit','var')
                        Logfile.RESP_PART1{trial,1} = 1;
                        Logfile.RESP_PART2{trial,1} = NaN;
                        Logfile.HITS_TONE1{trial,1} = 0;
                        Logfile.HITS_TONE2{trial,1} = 0;
                        Logfile.REACTION_TIME1{trial,1} = rt_tone1_miss;
                        Logfile.REACTION_TIME2{trial,1} = NaN; 
                    end

            case 1
                    %if hit for sinle tone;
                    if exist('rt_tone1_hit','var')
                        Logfile.RESP_PART1{trial,1} = key_press_count(find(key_press_rt == rt_tone1_hit ,1, 'first'));
                        Logfile.RESP_PART2{trial,1} = NaN;
                        Logfile.HITS_TONE1{trial,1} = 1;
                        Logfile.HITS_TONE2{trial,1} = NaN;
                        Logfile.REACTION_TIME1{trial,1} = rt_tone1_hit;
                        Logfile.REACTION_TIME2{trial,1} = NaN;

                    %if miss for single tone but keypress;
                    elseif ~exist('rt_tone1_hit','var') && rt_tone1_miss ~=0
                        Logfile.RESP_PART1{trial,1} = key_press_count(find(key_press_rt == rt_tone1_miss ,1, 'first'));
                        Logfile.RESP_PART2{trial,1} = NaN;
                        Logfile.HITS_TONE1{trial,1} = 0;
                        Logfile.HITS_TONE2{trial,1} = NaN;
                        Logfile.REACTION_TIME1{trial,1} = rt_tone1_miss;
                        Logfile.REACTION_TIME2{trial,1} = NaN;

                    %if miss for single tone and no keypress;
                    elseif ~exist('rt_tone1_hit','var') && rt_tone1_miss ==0
                        Logfile.RESP_PART1{trial,1} = 0;
                        Logfile.RESP_PART2{trial,1} = NaN;
                        Logfile.HITS_TONE1{trial,1} = 0;
                        Logfile.HITS_TONE2{trial,1} = NaN;
                        Logfile.REACTION_TIME1{trial,1} = rt_tone1_miss;
                        Logfile.REACTION_TIME2{trial,1} = NaN;
                    end

            case 2
                    %if hit for Tone1 and Tone2;
                    if exist('rt_tone1_hit','var') && exist('rt_tone2_hit','var')
                        Logfile.RESP_PART1{trial,1} = key_press_count(find(key_press_rt == rt_tone1_hit ,1, 'first'));
                        Logfile.RESP_PART2{trial,1} = key_press_count(find(key_press_rt == rt_tone2_hit ,1, 'first'));
                        Logfile.HITS_TONE1{trial,1} = 1;
                        Logfile.HITS_TONE2{trial,1} = 1;
                        Logfile.REACTION_TIME1{trial,1} = rt_tone1_hit;
                        Logfile.REACTION_TIME2{trial,1} = rt_tone2_hit;

                    %if miss for Tone1 (with keypress) and Tone2 (with keypress) ;
                    elseif ~exist('rt_tone1_hit','var')&& rt_tone1_miss ~= 0 && ~exist('rt_tone2_hit','var') && rt_tone2_miss ~= 0
                        Logfile.RESP_PART1{trial,1} = key_press_count(find(key_press_rt == rt_tone1_miss ,1, 'first'));
                        Logfile.RESP_PART2{trial,1} = key_press_count(find(key_press_rt == rt_tone2_miss ,1, 'first'));
                        Logfile.HITS_TONE1{trial,1} = 0;
                        Logfile.HITS_TONE2{trial,1} = 0;
                        Logfile.REACTION_TIME1{trial,1} = rt_tone1_miss;
                        Logfile.REACTION_TIME2{trial,1} = rt_tone2_miss;

                    %if miss for Tone1 (with keypress) and Tone2 (with no keypress) ;
                    elseif ~exist('rt_tone1_hit','var')&& rt_tone1_miss ~= 0 && ~exist('rt_tone2_hit','var') && rt_tone2_miss == 0
                        Logfile.RESP_PART1{trial,1} = key_press_count(find(key_press_rt == rt_tone1_miss ,1, 'first'));
                        Logfile.RESP_PART2{trial,1} = 0;
                        Logfile.HITS_TONE1{trial,1} = 0;
                        Logfile.HITS_TONE2{trial,1} = 0;
                        Logfile.REACTION_TIME1{trial,1} = rt_tone1_miss;
                        Logfile.REACTION_TIME2{trial,1} = 0;  

                    %if miss for Tone1 (with no keypress) and Tone2 (with keypress) ;
                    elseif ~exist('rt_tone1_hit','var')&& rt_tone1_miss == 0 && ~exist('rt_tone2_hit','var') && rt_tone2_miss ~= 0
                        Logfile.RESP_PART1{trial,1} = 0;
                        Logfile.RESP_PART2{trial,1} = key_press_count(find(key_press_rt == rt_tone2_miss ,1, 'first'));
                        Logfile.HITS_TONE1{trial,1} = 0;
                        Logfile.HITS_TONE2{trial,1} = 0;
                        Logfile.REACTION_TIME1{trial,1} = 0;
                        Logfile.REACTION_TIME2{trial,1} = rt_tone2_miss;

                    %if miss for Tone1 (with no keypress) and Tone2 (with no keypress) ;
                    elseif ~exist('rt_tone1_hit','var')&& rt_tone1_miss == 0 && ~exist('rt_tone2_hit','var') && rt_tone2_miss == 0
                        Logfile.RESP_PART1{trial,1} = 0;
                        Logfile.RESP_PART2{trial,1} = 0;
                        Logfile.HITS_TONE1{trial,1} = 0;
                        Logfile.HITS_TONE2{trial,1} = 0;
                        Logfile.REACTION_TIME1{trial,1} = 0;
                        Logfile.REACTION_TIME2{trial,1} = 0;

                    %if hit for Tone1 and miss Tone2 (with key press);
                    elseif exist('rt_tone1_hit','var') && ~exist('rt_tone2_hit','var') && rt_tone2_miss ~= 0 
                        Logfile.RESP_PART1{trial,1} = key_press_count(find(key_press_rt == rt_tone1_hit ,1, 'first'));
                        Logfile.RESP_PART2{trial,1} = key_press_count(find(key_press_rt == rt_tone2_miss ,1, 'first'));                      
                        Logfile.HITS_TONE1{trial,1} = 1;
                        Logfile.HITS_TONE2{trial,1} = 0;
                        Logfile.REACTION_TIME1{trial,1} = rt_tone1_hit;
                        Logfile.REACTION_TIME2{trial,1} = rt_tone2_miss;

                    %if hit for Tone1 and miss Tone2 (with no key press);
                    elseif exist('rt_tone1_hit','var') && ~exist('rt_tone2_hit','var') && rt_tone2_miss == 0 
                        Logfile.RESP_PART1{trial,1} = key_press_count(find(key_press_rt == rt_tone1_hit ,1, 'first'));
                        Logfile.RESP_PART2{trial,1} = 0;                      
                        Logfile.HITS_TONE1{trial,1} = 1;
                        Logfile.HITS_TONE2{trial,1} = 0;
                        Logfile.REACTION_TIME1{trial,1} = rt_tone1_hit;
                        Logfile.REACTION_TIME2{trial,1} = rt_tone2_miss;

                    %if miss for Tone1 and hit Tone2 (with key press);    
                    elseif ~exist('rt_tone1_hit','var') && rt_tone1_miss ~= 0 && exist('rt_tone2_hit','var') 
                        Logfile.RESP_PART1{trial,1} = key_press_count(find(key_press_rt == rt_tone1_miss ,1, 'first'));
                        Logfile.RESP_PART2{trial,1} = key_press_count(find(key_press_rt == rt_tone2_hit ,1, 'first'));                        
                        Logfile.HITS_TONE1{trial,1} = 0;
                        Logfile.HITS_TONE2{trial,1} = 1;
                        Logfile.REACTION_TIME1{trial,1} = rt_tone1_miss;
                        Logfile.REACTION_TIME2{trial,1} = rt_tone2_hit;

                    %if miss for Tone1 (with no key press) and hit Tone2;    
                    elseif ~exist('rt_tone1_hit','var') && rt_tone1_miss == 0 && exist('rt_tone2_hit','var') 
                        Logfile.RESP_PART1{trial,1} = 0;
                        Logfile.RESP_PART2{trial,1} = key_press_count(find(key_press_rt == rt_tone2_hit ,1, 'first'));                        
                        Logfile.HITS_TONE1{trial,1} = 0;
                        Logfile.HITS_TONE2{trial,1} = 1;
                        Logfile.REACTION_TIME1{trial,1} = rt_tone1_miss;
                        Logfile.REACTION_TIME2{trial,1} = rt_tone2_hit;    
                    end
        end


        %Now fill up the participant's information;
        Logfile.PARTICIPANT{trial,1} = respMat.PARTICIPANT{trial,1};
        Logfile.BLOCK{trial,1} = respMat.BLOCK{trial,1};
        Logfile.TRIAL_NUMB{trial,1} = respMat.TRIAL_NUMB{trial,1};
        Logfile.FREQPEAK{trial,1} = respMat.FREQPEAK{trial,1};
        Logfile.VIDEO{trial,1} = respMat.VIDEO{trial,1};
        Logfile.SNR{trial,1} = respMat.SNR{trial,1};
        Logfile.TONE{trial,1} = respMat.TONE{trial,1};
        Logfile.TONE_ONSET1{trial,1} = respMat.TONE_ONSET1{trial,1}; 
        Logfile.TONE_ONSET2{trial,1} = respMat.TONE_ONSET2{trial,1};
        Logfile.SCREEN_ONSET{trial,1} = respMat.SCREEN_ONSET{trial,1}; 
        Logfile.SCREEN_OFFSET{trial,1} = respMat.SCREEN_OFFSET{trial,1};
        Logfile.VA_LAG{trial,1} = respMat.VA_LAG(trial,1);

        clear count key_press_rt key_press_count rt_tone1_hit rt_tone2_hit rt_tone1_miss rt_tone2_miss 
        
        

    end
    
%Export the participant logfile;
respMat_ANALYSES = Logfile;

%change the video number in num format for further analyses;
for i =1:length(respMat_ANALYSES.VIDEO)
    temp_respMat_ANALYSES{i,1} = respMat_ANALYSES.VIDEO{i}(7:end);
end

respMat_ANALYSES.VIDEO = str2double(strtok(temp_respMat_ANALYSES,'_'));

%save and export participant's cleaned logfile;
cd(['XXX\subj' num2str(s)]);
save(['subj' num2str(s) '-Tone_detectionTask_clean'],'respMat_ANALYSES','alltrials','participant_parameters');

keep fs subjects
        
end


%%