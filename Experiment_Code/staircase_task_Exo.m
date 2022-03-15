% -------------------------------------------------------------------------
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%% >>>>>>>>>>>>>>>>>>>    N-Up-1-Down Staircasing    <<<<<<<<<<<<<<<<<<<<<<
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
% -------------------------------------------------------------------------
% **** See García-Pérez (1998) for more information about staircasing **** 
% -------------------------------------------------------------------------
function UpDn = staircase_task_Exo
% This function is meant to be used with OrientationWheelTask_v3
% There is a lot of code that is not needed for this task. This code
% orginally comes from tCAS_VisEntrain_Central_v5 and then the task was
% modified to better match the orientation wheel task. Therefore, there is
% a lot of unused code that is kept because it was part of one of those two
% original codes. Most of that unused code will be found in the local
% functions towards the end of this code (and will be saved in the prefs 
% structure). 
% 
% For changing the actual staircasing of each subject, look to the
% parameters found in the UpDn structure.
% 
% 
% 
% Information about typical stimulus color RGB values:
% (use uisetcolor to see them in a pop-up window)
% white = 255
% black = 0
% gray/background color = 127.5
% target color = 77.5
% mask color = 0
% entrainer color = 127.5 (no entrainers in this task)
% 
% /////////////////////////////////////////////////////////////////////////
%
% ==== TRIGGERS ====
%
% -- Targets Present --
% Trial start (fixation): 1
% Entrainers: 61-68
% Cue: 
%     left:  2
%     right: 3
% Target:
%     left:  30
%     right: 20
% Mask: 90
% Response screen Y/N: 25  
%   Response Yes: 180 
%   Response No: 185 
%   Response invalid: 150
% 
%
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

try

prepareEnvironment;

% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------    
% Input participant's number, name, date of testing, preferences
part_num = input('Participant Number:','s');
Filename = ['M:\Experiments\OrientWheel_Exo\Orient_Data\Staircasing\' part_num '_stair_exo.mat'];

% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////

window = openWindow(); %get info for psychtoolbox
prefs = getPreferences(); %get task variable info

% -------------------------------------------------------------------------   
% Get presentation timing information
refresh = Screen('GetFlipInterval', window.onScreen); % Get flip refresh rate
slack = refresh/2; % Divide by 2 to get slack

% -------------------------------------------------------------------------   
% Get rects for each item.
rects = cell(1, max(prefs.setSizes));
for i = 1:max(prefs.setSizes)
    rects{i} = circularArrayRects([0, 0, prefs.squareSize, prefs.squareSize], ...
        i, prefs.radius, window.centerX, window.centerY)';
end

% -------------------------------------------------------------------------
% Set default colors of stimuli
prefs.mask_gray = 0; %default is black
prefs.targ_gray = 0; %default is black

% -------------------------------------------------------------------------
% Center the target oval & cues on the centre of the screen
centeredRectresp = CenterRectOnPointd(prefs.baseCircleresp,window.centerX,window.centerY);
centeredRectrespR = CenterRectOnPointd(prefs.baseCircleresp,window.centerXR,window.centerY); 
centeredRectrespL = CenterRectOnPointd(prefs.baseCircleresp,window.centerXL,window.centerY);

% -------------------------------------------------------------------------
% Create offscreen window with the texture of the mask
maskwin = createMasktex(window,prefs); 

% ------------------------------------------------------------------------- 
% Get coordinates of the place holder
placeRect = createPlaceHolder(window,prefs);

% ------------------------------------------------------------------------- 
% Create offscreen window with the cues
CueL = createCueL(window); %create offscreen window with cue
CueR = createCueR(window); %create offscreen window with cue
  
% -------------------------------------------------------------------------  
% Put up instructions for staircasing task and wait for keypress.
instruct_staircase(window,prefs);

% -------------------------------------------------------------------------    
% Location of orientations on screen
orientWheelLocations = orientwheelLocations(window,prefs);

% -------------------------------------------------------------------------  

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% Set up the random mix of lags and target position and cue for each block

ntotaltrial = prefs.nTrials;

% Set up the random mix of lags for each block
all_lags = [1:prefs.n_lags];
if prefs.lags_per_block > 1
    for i_lag = 2:prefs.lags_per_block
        all_lags = [all_lags 1:prefs.n_lags];
    end
end
i_lags = randperm(prefs.lags_per_block * (prefs.n_lags));
all_lags = all_lags(i_lags);
clear i_lag


%set up the right/left target location
position = 1;
for i_pos = 2:ntotaltrial
    if i_pos <= round(ntotaltrial/2)
        position = [position 1]; % 1 is right
    else
        position = [position 0]; % 0 is left
    end
end
position = position(randperm(ntotaltrial));


%set up the valid/invalid cues (1=informative)
valid = 1;
validtrials = ntotaltrial*prefs.p_validity;
for i_valid = 2:ntotaltrial
    if i_valid > round(validtrials)
        valid = [valid 0]; 
    else
        valid = [valid 1];
    end
end
clear i_valid
valid = valid(randperm(ntotaltrial));
    

% Pre-allocate some variables
allCoords_targ = cell(1,ntotaltrial);
trialOrient = cell(1,ntotaltrial);  


% -------------------------------------------------------------------------   
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% +++++++++++++     Set-Up Structure UpDn for Staircasing     +++++++++++++
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% -------------------------------------------------------------------------   
% Note. DECREASE in stimulus level = INCREASE in ability to see target 
%       (staircase procedures usually assume this is a positive relationship) 
%       >> up step = decrease difficulty = closer to black (RGB=0)
%       >> down step = increase difficulty = closer to background (RGB=128)

% -------------------------------------------------------------------------
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% Modifiable paramaters of the staircasing
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% -------------------------------------------------------------------------

UpDn.xStartLevel = 77.5; %Initial target color: ~2/3 between background and black to make it easy
   
UpDn.xCurrent = UpDn.xStartLevel; %Stimulus level to be used on current trial

UpDn.sUp   = 1; % Number of consecutive undetected responses after which stimulus difficulty decreases 
%    <<< after n trials with target undetected: target color - sSizeUp = darker target color (= easier to see)
UpDn.sDown = 2; % Number of consecutive detected responses after which stimulus difficulty increases
%    <<< after n trials with target detected: target color + sSizeDown = lighter target color (= closer to background = harder to see)


% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% ``` Step Sizes ```
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Size of decrease in stimulus difficulty (-sSizeUp to reduce RGB value = darker target color).
UpDn.sSizeUp = [3 4 6]; %aiming for .55-.60: [2 6 4]
% Size of increase in stimulus difficulty (+sSizeDown to increase RGB value = lighter target color).
UpDn.sSizeDown   = [3 5 5]; %aiming for .55-.60: [4 11 9]
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% ---see García-Pérez (1998) about the ratio of sSizeDown/sSizeUp when using fixed step sizes--- 
% Use formula to determine step sizes to get targets proportion correct 
%  ii=1;
%  targetP = (UpDn.sSizeUp(ii)./(UpDn.sSizeUp(ii)+UpDn.sSizeDown(ii))).^(UpDn.sUp./UpDn.sDown)
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% *** Desired Proportion of Targets Detected ***
% UpDn.propDetect = .50;
UpDn.propDetect = .70;
UpDn.currProp = 0; % starting value for i_sRun = 1


% *** Define Rules for Termination of Staircasing ***
UpDn.stopType = 'trials'; % Can be either ‘trials’ or ‘revels’
        % When set to ‘trials’, staircase will terminate after the number of trials in stopRule.
        % When set to ‘revels’, staircase will terminate after the number of reversals in stopRule.
            % note. code is not currently set-up for reversals
            
% stopRule = Number of trials or reversals before run terminates
    % set-up to have the total # of staircase trials equal to two experiment blocks                  
UpDn.nStairsRun = length(UpDn.sSizeUp); % # of blocks task will run 
% want the total number of trials to be = 2*trials_per_block 
UpDn.stopRule = round((2*prefs.trials_per_block) ./ UpDn.nStairsRun);
       

% Max stimulus level to be used in staircase.
    %if set to empty ([]) no max is applied.
UpDn.xStimMax = window.gray - 2;
% Min stimulus intensity to be used in staircase.
    %if set to empty ([]) no min is applied. 
UpDn.xStimMin = 0; %black

trialIndex = 1; %trial counter (for compatability)

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%% @@@@@@@@@@@@@@@@ Loop for multiple runs of staircasing @@@@@@@@@@@@@@@@@
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
for i_sRun = 1:UpDn.nStairsRun
    
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%         UpDn storage fields of each trials staircasing result
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::    
    
    i_trial = 1; % restart at trial 1 
    
     % Contains stimulus level to be used on current trial    
    UpDn.xCurrent = UpDn.xStartLevel;

    % Stores stimulus level on trials, and passes them to to UpDn.x 
    UpDn.x(i_sRun,i_trial) = UpDn.xCurrent;
    
    % Stores stimulus level on trials
    UpDn.xStairs(i_sRun,i_trial) = UpDn.x(i_sRun,i_trial);

    UpDn.reverse(i_sRun,i_trial) = 0; % Stores for each trial whether a reversal occurred. 
            % It contains a 0 for trials where no reversal occurred.
            % It contains the count of the reversal for trials where a reversal did occur. 
            
    % Counts the # of consecutive incorrect responses until Up rule is met (equals UpDn.sUp)
    UpDn.u = 0;   % Resets to 0 on the next trial 
    % Counts the # of consecutive correct responses until Down rule is met (equals UpDn.sDown) 
    UpDn.d = 0;   % Resets to 0 on the next trial 

    % *** Used as a Termination Flag *** 
    UpDn.stop = 0; %reset to 0
        % While stop criterion has not been reached, UpDn.stop = 0 
        % When stop criterion is reached, UpDn.stop = 1  
        
     % Get proportion detected before each run after the first run  
     if i_sRun > 1 
        [mR,nR] = size(UpDn.stairResp);
        prev_trials = (mR.*nR); %total trials from previous
     end   

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
% /////////////////////////////////////////////////////////////////////////   
%% === Run staircasing procedure until meet stop criteria ===
% /////////////////////////////////////////////////////////////////////////
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



%     for i_trial = 1:UpDn.stopRule
    while ~UpDn.stop  %loop until stop criteria is met     
        
% ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
        % update record of proportion targets seen if not first block
        if i_sRun > 1
            [mR,nR] = size(UpDn.stairResp);
            UpDn.currProp = sum(sum(UpDn.stairResp,2))./(prev_trials + (i_trial-1));  
            UpDn.detectP(i_sRun-1,i_trial) = UpDn.currProp; %to be saved
        end
% ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
        
        % Get lag id # for triggersw
        lag = prefs.lags(all_lags(i_trial))/5;
        
        % Determine how many items there are on this trial and the duration (this is left over code from original task)
        nItems = 1; %should always be 1 because there will only be 1 target
        retentionInterval = prefs.retentionIntervals;
        
        % Pick the orientation for this trial.
        trialOrient{trialIndex} = prefs.degslocs(prefs.selectorientTrial(i_trial)); 
        
% ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 

% ========================================================================= 

%       $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%       ++++++ Run Staircasing Task ++++++
%       $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

% ========================================================================= 
        %% Present Fixation
        pause(window); 
        Screen('FillRect', window.onScreen, window.gray);
        Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255);
        Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect); %placeholder
        Screen('FillRect',window.onScreen, Vpixx2Vamp(1), prefs.trigger_size);
        t_fixate_onset = Screen('Flip', window.onScreen);  
% ========================================================================= 
        %% Interval (refresh screen)
        Screen('FillOval', window.onScreen, window.gray, rects{nItems});
        Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255); %keep fixation till target
        Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect); %placeholder
        Screen('FillRect',window.onScreen,Vpixx2Vamp(0),prefs.trigger_size);
        Screen('Flip', window.onScreen);
% =========================================================================         
        %% Cue
        if position(trialIndex) == 0 
            % present the Left cue
            Screen('DrawTexture', window.onScreen, CueL, [], [], [], 0);
            Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255);%draw fixation
            Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect); %placeholder
            Screen('FillRect',window.onScreen,Vpixx2Vamp(2),prefs.trigger_size);
            tcue_onset = Screen('Flip',window.onScreen,t_fixate_onset + prefs.fixation_length*refresh - slack);
            prefs.fixation_time(trialIndex) = tcue_onset - t_fixate_onset;
        else
            % present the Right cue
            Screen('DrawTexture', window.onScreen, CueR, [], [], [], 0);
            Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255);%draw fixation
            Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect); %placeholder
            Screen('FillRect',window.onScreen,Vpixx2Vamp(3),prefs.trigger_size);
            tcue_onset = Screen('Flip',window.onScreen,t_fixate_onset + prefs.fixation_length*refresh - slack);
            prefs.fixation_time(trialIndex) = tcue_onset - t_fixate_onset;
        end
        
        % Post-cue interval
        Screen('FillOval', window.onScreen, window.gray, rects{nItems});
        Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255); %keep fixation till target
        Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect); %placeholder
        Screen('FillRect',window.onScreen,Vpixx2Vamp(0),prefs.trigger_size);
        tblank_onset = Screen('Flip',window.onScreen,tcue_onset + prefs.cue_length*refresh - slack);
% =========================================================================       
        %% Entrainers
        if prefs.n_entrs > 0
            Screen('FillOval', window.onScreen, window.gray, rects{nItems});
            Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255); %keep fixation till target
            Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect); %placeholder
            Screen('FillRect',window.onScreen,Vpixx2Vamp(61),prefs.trigger_size);
            tentr_onset = Screen(window.onScreen, 'Flip', tblank_onset + prefs.preblank_length*refresh - slack);
            if prefs.n_entrs > 1
                for i_entr = 2:prefs.n_entrs
                    Screen('FillOval', window.onScreen, window.gray, rects{nItems});
                    Screen('FillRect',window.onScreen,Vpixx2Vamp(0),prefs.trigger_size);
                    tblank_onset = Screen(window.onScreen, 'Flip', tentr_onset + prefs.entr_length*refresh - slack);
                    Screen('FillOval', window.onScreen,  (prefs.entr_grey+window.gray),...
                        [rects{nItems}(1)-prefs.maskwidth, rects{nItems}(2)-prefs.maskwidth,rects{nItems}(3)+prefs.maskwidth,...
                        rects{nItems}(4)+prefs.maskwidth]);
                    Screen('FillOval', window.onScreen, window.gray, rects{nItems});
                    Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect); %placeholder
                    Screen('FillRect',window.onScreen,Vpixx2Vamp(60 + i_entr),prefs.trigger_size);
                    tentr_onset = Screen(window.onScreen, 'Flip', tblank_onset + prefs.entr_gap_length*refresh - slack);
                end
            end
        end       
% =========================================================================         
            %% Lag
            Screen('FillOval', window.onScreen, window.gray, rects{nItems});
            Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255); %keep fixation till target
            Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect); %placeholder
            Screen('FillRect',window.onScreen,Vpixx2Vamp(0),prefs.trigger_size);
            tlag_onset = Screen(window.onScreen, 'Flip', tentr_onset + prefs.entr_length*refresh - slack);
% =========================================================================          
            %% Draw the target stimulus
            yadj = prefs.orientwheel(2,trialOrient{trialIndex}); %adjust lines y position 
            xadj = prefs.orientwheel(1,trialOrient{trialIndex}); %adjust lines x position 
            yadj2 = prefs.orientwheel2(2,trialOrient{trialIndex}); %adjust lines y position 
            xadj2 = prefs.orientwheel2(1,trialOrient{trialIndex}); %adjust lines x position 
            % Set the coordinates (these are all relative to zero we will let
            % the drawing routine center it for us)
            xCoords = [xadj -xadj2 xadj xadj2 xadj 0];
            yCoords = [yadj yadj2 yadj -yadj2 yadj 0];
            allCoords_targ{trialIndex} = [xCoords; yCoords];
% ========================================================================= 
            %% Right Target
            if position(trialIndex) == 1 %right target
                % Draw the orientation lines, set it to the center of our screen and set good quality antialiasing
                Screen('DrawLines', window.onScreen, allCoords_targ{trialIndex}, prefs.lineWidthPix, UpDn.xCurrent, [window.centerXR window.centerY], 2);
                % Draw the central circle to the screen
                Screen('FillOval', window.onScreen, UpDn.xCurrent, centeredRectrespR');
                % Draw fixation
                Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255);
                Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect); %placeholder
                % Draw trigger
                Screen('FillRect',window.onScreen,Vpixx2Vamp(20),prefs.trigger_size);

% =========================================================================                 
            else
            %% Left Target
                % Draw the orientation lines, set it to the center of our screen and set good quality antialiasing
                Screen('DrawLines', window.onScreen, allCoords_targ{trialIndex}, prefs.lineWidthPix, UpDn.xCurrent, [window.centerXL window.centerY], 2);
                % Draw the central circle to the screen
                Screen('FillOval', window.onScreen, UpDn.xCurrent, centeredRectrespL');
                % Draw fixation
                Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255);
                Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect); %placeholder
                % Draw trigger
                Screen('FillRect',window.onScreen,Vpixx2Vamp(30),prefs.trigger_size);
            
            end       
% =========================================================================              
            %% Post the stimulus
            ttarget_onset = Screen('Flip', window.onScreen, tlag_onset + prefs.lagISI(all_lags(trialIndex))*refresh - slack);
% =========================================================================               
            %% Interval
            Screen('FillOval', window.onScreen, window.gray, rects{nItems});
            %draw fixation
            Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255);
            Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect); %placeholder
            Screen('FillRect',window.onScreen,Vpixx2Vamp(0),prefs.trigger_size);
            tISI_onset = Screen('Flip', window.onScreen, ttarget_onset + prefs.targ_length*refresh - slack);        
% =========================================================================               
            %% Mask
            Screen('DrawTexture', window.onScreen, maskwin, [], [], [], 0);
            %draw fixation
            Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255);
            Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect); %placeholder
            Screen('FillRect',window.onScreen,Vpixx2Vamp(90),prefs.trigger_size);
            t_maskonset = Screen('Flip', window.onScreen, tISI_onset + prefs.maskISI*refresh - slack);
% =========================================================================              
            %% Retention Interval
            Screen('FillOval', window.onScreen, window.gray, rects{nItems});
            Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255); %Draw fixation
            Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect); %placeholder
            Screen('FillRect',window.onScreen,Vpixx2Vamp(0),prefs.trigger_size);
            Screen('Flip', window.onScreen, t_maskonset + prefs.mask_length*refresh - slack);

            % Retention interval
            WaitSecs(retentionInterval);
% ========================================================================= 
           %% Response            
           UpDn.presentedOrientDegrees(i_sRun,i_trial) = trialOrient{trialIndex}; 
           
           Screen('FillRect',window.onScreen, Vpixx2Vamp(25), prefs.trigger_size);
           Screen('Flip', window.onScreen)
           %one more refresh to clear triggers
           Screen('FillRect',window.onScreen,window.black,prefs.trigger_size);
           Screen('Flip', window.onScreen);
           
           t1 = GetSecs; %record for RT
           % Get key press
           keyIsDown = 0;
           while  ~keyIsDown
               [keyIsDown, secs, keyCode] = KbCheck;
           end
           
          % Keep a log of the subject's responses
          response = find(keyCode>0);   %1 is 49, 5 is 53, left arrow is 37, right is 39     
          % detected (left arrow)
            if response == 37
                % Stores responses for all staircasing trials
                UpDn.stairResp(i_sRun,i_trial) = 1; %detected
                Screen('FillRect',window.onScreen, Vpixx2Vamp(180), prefs.trigger_size);
                res_trig = Screen('Flip',window.onScreen,[],0);
          % undetected (right arrow)
            elseif response == 39
                % Stores responses for all staircasing trials
                UpDn.stairResp(i_sRun,i_trial) = 0; %undetected
                Screen('FillRect',window.onScreen, Vpixx2Vamp(185), prefs.trigger_size);
                res_trig = Screen('Flip',window.onScreen,[],0);
            else
                % Stores responses for all staircasing trials
                UpDn.stairResp(i_sRun,i_trial) = 9; %NULL
                Screen('FillRect',window.onScreen, Vpixx2Vamp(150), prefs.trigger_size);
                res_trig = Screen('Flip',window.onScreen,[],0);
            end
           
         % Blank
%          Screen('CopyWindow',blank ,window,rect,rect);
         Screen('FillRect',window.onScreen,Vpixx2Vamp(0),prefs.trigger_size);
         Screen(window.onScreen,'Flip',res_trig + refresh - slack);
         
            
         
% ========================================================================= 

         % Compute and store response time
         UpDn.subject_rt(i_sRun,i_trial) = secs-t1; 
            
            

        
% /////////////////////////////////////////////////////////////////////////
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%% <<<<<<<<<<<<<<<<    Update UpDn Structure    >>>>>>>>>>>>>>>>>>>>>>>>>>>
% /////////////////////////////////////////////////////////////////////////
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

 %% Target Detected (step down = increase difficulty = darken target)
        if i_trial == 1
            if UpDn.stairResp(i_sRun,i_trial) == 1
               % If a reversal happened on the current trial, contains the direction of that reversal 
                UpDn.direction = -1;       
            else
                UpDn.direction = 1;
            end
        end

        %% IF DETECTED...
        if UpDn.stairResp(i_sRun,i_trial) == 1
            UpDn.d = UpDn.d + 1; %detections counter
            
            % NO CHANGE IN TARGET COLOR:
            % If reach desired proportion and step rule has not been met
%             if ((UpDn.currProp == UpDn.propDetect) && i_sRun > 1) && UpDn.d ~= UpDn.sDown %No step down change
            if UpDn.d ~= UpDn.sDown %No step down change if step rule has not been met
                UpDn.xStairs(i_sRun,i_trial+1) = UpDn.xStairs(i_sRun,i_trial); %for next trial
                
            % CHANGE IN TARGET COLOR:
            % ~~TRUE if d counter == step down rule ||or|| reversal has not yet occurred 
            % ||or|| after 1st run, proportion detected > desired proportion detected    
%             elseif UpDn.d == UpDn.sDown || max(UpDn.reverse(i_sRun),[],2) < 1 ||...
%                     ((UpDn.currProp > UpDn.propDetect) && i_sRun > 1)
            elseif UpDn.d == UpDn.sDown || max(UpDn.reverse(i_sRun),[],2) < 1    
                %% Step down change = target color + step down size (increasing color #)
                % New target color
                UpDn.xStairs(i_sRun,i_trial+1) = UpDn.xStairs(i_sRun,i_trial) + UpDn.sSizeDown(i_sRun);
                if UpDn.xStairs(i_sRun,i_trial+1) > UpDn.xStimMax %target can't be darker than background
                   UpDn.xStairs(i_sRun,i_trial+1) = UpDn.xStimMax; 
                end
                UpDn.u = 0; %reset when there is a reversal
                UpDn.d = 0; %reset when there is a reversal
                UpDn.reverse(i_sRun,i_trial) = 0;
                if UpDn.direction == 1
                    UpDn.reverse(i_sRun,i_trial) = sum(UpDn.reverse(i_sRun,i_trial)~=0) + 1;
                else
                    UpDn.reverse(i_sRun,i_trial) = 0;
                end
                UpDn.direction = -1;
      
            else %No step down change
                UpDn.xStairs(i_sRun,i_trial+1) = UpDn.xStairs(i_sRun,i_trial);
            end 

        else %% IF UNDETECTED...
            UpDn.u = UpDn.u + 1; %undetected counter
            
            % NO CHANGE IN TARGET COLOR:
            % If reach desired proportion or step rule has not been met
%             if ((UpDn.currProp == UpDn.propDetect) && i_sRun > 1) && UpDn.u ~= UpDn.sUp %No step down change
             if UpDn.u ~= UpDn.sUp %No step down change if step rule has not been met  
                UpDn.xStairs(i_sRun,i_trial+1) = UpDn.xStairs(i_sRun,i_trial);
                
            % CHANGE IN TARGET COLOR:    
            % ~~TRUE if d counter == step up rule ||or|| reversal has not yet occurred
%             elseif UpDn.u == UpDn.sUp ||  max(UpDn.reverse(i_sRun),[],2) < 1 ||...
%                     ((UpDn.currProp < UpDn.propDetect) && i_sRun > 1)
            elseif UpDn.u == UpDn.sUp ||  max(UpDn.reverse(i_sRun),[],2) < 1   
                %% Step up change = target color - step up size (decrease color #)
                % New target color
                UpDn.xStairs(i_sRun,i_trial+1) = UpDn.xStairs(i_sRun,i_trial)- UpDn.sSizeUp(i_sRun);
                if UpDn.xStairs(i_sRun,i_trial+1) < UpDn.xStimMin %target's color limit is white
                    UpDn.xStairs(i_sRun,i_trial+1) = UpDn.xStimMin;
                end
                UpDn.u = 0; %reset when there is a reversal
                UpDn.d = 0; %reset when there is a reversal
                UpDn.reverse(i_sRun,i_trial) = 0;
                if UpDn.direction == -1
                    UpDn.reverse(i_sRun,i_trial) = sum(UpDn.reverse(i_sRun,i_trial)~=0) + 1;
                else
                    UpDn.reverse(i_sRun,i_trial) = 0;
                end
                UpDn.direction = 1;
            else %No step up change
                UpDn.xStairs(i_sRun,i_trial+1) = UpDn.xStairs(i_sRun,i_trial);
            end    
        end

        % When to stop this run
        if strncmpi(UpDn.stopType,'reversals',4) &&...
                sum(UpDn.reversal(i_sRun,i_trial)~=0) == UpDn.stopRule
            UpDn.stop = 1;
        end
        if strncmpi(UpDn.stopType,'trials',4) && i_trial == UpDn.stopRule
            UpDn.stop = 1;
        end
        
        % If we didn't stop, update x variable
        if ~UpDn.stop
            UpDn.x(i_sRun,i_trial+1) = UpDn.xStairs(i_sRun,i_trial+1);
            if UpDn.x(i_sRun,i_trial+1) >= UpDn.xStimMax %limit color to darker than backgorund
                UpDn.x(i_sRun,i_trial+1) = UpDn.xStimMax;
            elseif UpDn.x(i_sRun,i_trial+1) < UpDn.xStimMin %limit color to lighter than black
                UpDn.x(i_sRun,i_trial+1) = UpDn.xStimMin;
            end
            UpDn.xCurrent = UpDn.x(i_sRun,i_trial+1); %update
        end

        i_trial = 1 + i_trial; %UpDn trial counter
        trialIndex = 1 + trialIndex; %Orientwheel trial counter
        
        WaitSecs(0.5); %wait 1/2 second before next trial
        
        
    end %while loop for meeting stop criteria
    
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||   
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% For Reference: determine and display targeted proportion correct and stimulus intensity
%     targetP = (UpDn.sSizeUp(i_sRun)./(UpDn.sSizeUp(i_sRun)+UpDn.sSizeDown(i_sRun))).^(1./UpDn.sDown)
%     targetX = PAL_Gumbel(trueParams, targetP,'inverse'); 
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||   
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    if i_sRun < UpDn.nStairsRun %take a break
        rest(window);
    elseif i_sRun == UpDn.nStairsRun %end of task
        finish(window);
    end

     
end %Loop for running multiple staircasing blocks
clear i_trial i_sRun

% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%% Record info for later analysis
% Record that target was present
    % Total number of trials across all runs = i_trial
UpDn.subject_present(1,1:trialIndex) = 1;

% Record responses for later analysis
strial = 1;
for run = 1:UpDn.nStairsRun
    for tri = 1:UpDn.stopRule
        UpDn.subject_answer_all(1,strial) = UpDn.stairResp(run,tri);
        UpDn.subject_rt_all(1,strial) = UpDn.subject_rt(run,tri);
        strial = strial + 1;
    end
    clear tri
end
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
%% Final target color (average of all the runs)

xtarget = round(mean(UpDn.x(1:UpDn.nStairsRun,length(UpDn.x))));

% if our final proportion detected is off
% detdiff = round(mean(UpDn.detectP(1:size(UpDn.detectP),length(UpDn.x)))) - UpDn.propDetect;
% if detdiff > 0.08 %detection rate too high
%     xtarget = xtarget + 5; %5 points closer to grey
% elseif detdiff < 0.08 %detection rate too low
%     xtarget = xtarget - 5; %5 points darker than grey
% end


% To later record the target color
UpDn.targ_grey = xtarget;



% =========================================================================     
% =========================================================================  

save(Filename, 'prefs', 'UpDn', 'xtarget');
postpareEnvironment;

% =========================================================================     
% =========================================================================     
catch
    postpareEnvironment;
    psychrethrow(psychlasterror);
    
end % end try/catch

% ========================================================================= 
% ========================================================================= 
end % end whole script

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% #########################################################################
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ========================================================================= 
% -------------------------------------------------------------------------
% #########################################################################
% #########################################################################
% -------------------------------------------------------------------------
% ========================================================================= 


% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% #########################################################################
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function prepareEnvironment

ccc

HideCursor; % Comment out when debugging

commandwindow; % Select the command window to avoid typing in open scripts

% Seed the random number generator.
% RandStream.setDefaultStream(RandStream('mt19937ar', 'seed', sum(100*clock)));
RandStream.setGlobalStream(RandStream('mt19937ar', 'seed', sum(100*clock)));

% ListenChar(2); % Don't print to MATLAB command window
% Screen('Preference', 'SkipSyncTests', 1); %choose 1 for test on an LCD

% /////////////////////////////////////////////////////////////////////////
%% Set up parallel port (for when triggers are sent with Vpixx2Vamp)
%initialize the inpoutx64 low-level I/O driver
config_io;
%optional step: verify that the inpoutx64 driver was successfully installed
global cogent;
if( cogent.io.status ~= 0 )
    error('inp/outp installation failed');
end
%write a value to the default LPT1 printer output port (at 0x378)
address_eeg = hex2dec('B010');
outp(address_eeg,0);  %set pins to zero
% /////////////////////////////////////////////////////////////////////////
end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------
%% End of experiment
function postpareEnvironment
ShowCursor;
%   ListenChar(0);
Screen('CloseAll');
end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------
%% Gives the subject a set of instructions about the staircasing task
function instruct_staircase(window,prefs)

%Screen 1
% prefs = getPreferences();
Screen('TextSize', window.onScreen, window.fontsize);
Screen('DrawText', window.onScreen, 'You are going to be asked to perform a visual task.',(window.centerX-400),(window.centerY+20), 255);
Screen('DrawText', window.onScreen, 'The next few screens will briefly explain the first task.',...
    (window.centerX-400),(window.centerY+50), 255); %line every +30
Screen('DrawText', window.onScreen, 'On each screen press any key to continue.',(window.centerX-400),(window.centerY+80), 255); %line every +25
Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
Screen('Flip', window.onScreen);
% GetClicks(window.onScreen);
WaitSecs(1);
KbWait %wait for subject to press button


%Screen 2
Screen('TextSize', window.onScreen, window.fontsize);
Screen('DrawText', window.onScreen, 'A fixation dot will appear. Please remain focused on the white central dot during the ENTIRE task.',...
    (window.centerX-600),(window.centerY+50), 255);
Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255); %255 is text color (white)
Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
%place holder
placeRect = createPlaceHolder(window,prefs);
Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect);
Screen('Flip', window.onScreen);
WaitSecs(1);
KbWait %wait for subject to press button


%Screen 2.2
Screen('TextSize', window.onScreen, window.fontsize); 
% CueL = createCueL(window,prefs); %create offscreen window with cue
CueR = createCueR(window); %create offscreen window with cue
Screen('DrawTexture', window.onScreen, CueR, [], [], [], 0);
Screen('DrawText', window.onScreen, 'An arrow will then briefly appear.',...
    (window.centerX-400),(window.centerY+50), 255);
Screen('DrawText', window.onScreen, 'The arrow may be pointing right...',...
    (window.centerX-400),(window.centerY+80), 255); %line every +30
%place holder
placeRect = createPlaceHolder(window,prefs);
Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect);
%put stimulus and text on screen
Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
Screen('Flip', window.onScreen);
WaitSecs(1);
KbWait %wait for subject to press button


%Screen 2.3
Screen('TextSize', window.onScreen, window.fontsize); 
CueL = createCueL(window); %create offscreen window with cue
% CueR = createCueR(window,prefs); %create offscreen window with cue
Screen('DrawTexture', window.onScreen, CueL, [], [], [], 0);
Screen('DrawText', window.onScreen, 'Or, the arrow may be pointing left',...
    (window.centerX-400),(window.centerY+50), 255);
Screen('DrawText', window.onScreen, 'The arrow indicates the side a target stimulus will most likely appear.',...
    (window.centerX-400),(window.centerY+80), 255); %line every +30
Screen('DrawText', window.onScreen, 'While keeping your eyes focused on the center,',...
    (window.centerX-400),(window.centerY+110), 255); %line every +30
Screen('DrawText', window.onScreen, 'pay attention COVERTLY (without moving your eyes) to the circle on the side the arrow is pointing.',...
    (window.centerX-400),(window.centerY+140), 255); %line every +30
%place holder
placeRect = createPlaceHolder(window,prefs);
Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect);
%put stimulus and text on screen
Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
Screen('Flip', window.onScreen);
WaitSecs(1);
KbWait %wait for subject to press button



%Screen 3
Screen('TextSize', window.onScreen, window.fontsize);
Screen('DrawText', window.onScreen, 'The arrow will disappear and a target pointing in a certain direction will briefly appear inside the circle on the left OR right.',...
    (window.centerX-600),(window.centerY+50), 255);
Screen('DrawText', window.onScreen, 'Your task is to try to detect which direction the target is pointing.',...
    (window.centerX-400),(window.centerY+80), 255); %line every +30
Screen('DrawText', window.onScreen, 'It might be difficult to see the target, but still try your best.',...
    (window.centerX-400),(window.centerY+110), 255); %line every +30
%---draw target stimulus---
xCoords = [15 1.73 15 -1.73 15 0];
yCoords = [-25.98 1 -25.98 -1 -25.98 0];
allCoords_targ = [xCoords; yCoords];
%right target
Screen('DrawLines', window.onScreen, allCoords_targ, prefs.lineWidthPix, prefs.targ_gray, [window.centerXR window.centerY], 2);
centeredRectresp = CenterRectOnPointd([0,0,10,10], window.centerXR, window.centerY); %location of center for oval
Screen('FillOval', window.onScreen, prefs.targ_gray, centeredRectresp');
%left target
Screen('DrawLines', window.onScreen, allCoords_targ, prefs.lineWidthPix, prefs.targ_gray, [window.centerXL window.centerY], 2);
centeredRectresp = CenterRectOnPointd([0,0,10,10], window.centerXL, window.centerY); %location of center for oval
Screen('FillOval', window.onScreen, prefs.targ_gray, centeredRectresp');
%place holder
placeRect = createPlaceHolder(window,prefs);
Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect);
%put stimulus and text on screen
Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
Screen('Flip', window.onScreen);
WaitSecs(1);
KbWait %wait for subject to press button


%Screen 4
Screen('TextSize', window.onScreen, window.fontsize);
maskwin = createMasktex(window,prefs); %create offscreen window with mask
Screen('DrawTexture', window.onScreen, maskwin, [], [], [], 0);
Screen('DrawText', window.onScreen, 'A star-like shape will then quickly appear and then disappear.',...
    (window.centerX-400),(window.centerY+50), 255);
Screen('DrawText', window.onScreen, 'Remember that you want to detect the pointing target, not the star.',...
    (window.centerX-400),(window.centerY+80), 255); %line every +30
%place holder
placeRect = createPlaceHolder(window,prefs);
Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect);
Screen('FillRect',window.onScreen,Vpixx2Vamp(0),prefs.trigger_size);
Screen('Flip', window.onScreen);
WaitSecs(1);
KbWait %wait for subject to press button


%Screen 5
%draw target stimulus
xCoords = [-21.21 -1.41 -21.21 1.41 -21.21 0];
yCoords = [21.21 -1.41 21.21 1.41 21.21 0];
allCoords_targ = [xCoords; yCoords];
Screen('DrawLines', window.onScreen, allCoords_targ, prefs.lineWidthPix, prefs.targ_gray, [window.centerX window.centerY], 2);
centeredRectresp = CenterRectOnPointd([0,0,10,10], window.centerX, window.centerY); 
Screen('FillOval', window.onScreen, prefs.targ_gray, centeredRectresp');
Screen('TextSize', window.onScreen, window.fontsize);
Screen('DrawText', window.onScreen, 'Press the LEFT arrow if you can see the target.',(window.centerX-880),(window.centerY-50), 255);
Screen('DrawText', window.onScreen, 'Press the RIGHT arrow if you only see the star',(window.centerX-880),(window.centerY), 255);
Screen('DrawText', window.onScreen, 'If you are unsure whether the pointing target was there,',(window.centerX-880),(window.centerY+60), 255);
Screen('DrawText', window.onScreen, 'please provide your best guess.',(window.centerX-880),(window.centerY+90), 255);
Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
Screen('Flip', window.onScreen);
WaitSecs(1);
KbWait %wait for subject to press button


%Screen 6
Screen('TextSize', window.onScreen, window.fontsize);
Screen('DrawText', window.onScreen, 'Remember...Press the LEFT arrow if you see any indication of the pointing target being present.',(window.centerX-400),(window.centerY+20), 255);
Screen('DrawText', window.onScreen, 'Press the RIGHT arrow if you do NOT see the pointing target.',(window.centerX-400),(window.centerY+50), 255);
Screen('DrawText', window.onScreen, 'Press any key when you are ready to start.',(window.centerX-400),(window.centerY+80), 255);
Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
Screen('Flip', window.onScreen);
WaitSecs(1);
KbWait %wait for subject to press button


end %end function instruct_staircase

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------

function pause(window)
prefs = getPreferences();
Screen('FillRect', window.onScreen, window.gray);
Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
Screen('Flip', window.onScreen);
WaitSecs(0.5);
end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------

function rest(window)
prefs = getPreferences();

%Screen 1
Screen('TextSize', window.onScreen, window.fontsize);
Screen('DrawText', window.onScreen, 'Feel free to take a break at this time.',...
    (window.centerX-400),(window.centerY+20), 255);
Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
Screen('Flip', window.onScreen);
% GetClicks(window.onScreen);
WaitSecs(1);
KbWait %wait for subject to press button

%Screen 2
Screen('TextSize', window.onScreen, window.fontsize);
Screen('DrawText', window.onScreen, 'Remember...Press the LEFT arrow if you see any indication of the pointing target being present.',(window.centerX-400),(window.centerY+20), 255);
Screen('DrawText', window.onScreen, 'Press the RIGHT arrow if you do NOT see the pointing target.',(window.centerX-400),(window.centerY+50), 255);
Screen('DrawText', window.onScreen, 'Press any key when you are ready to continue.',(window.centerX-400),(window.centerY+80), 255);
Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
Screen('Flip', window.onScreen);
% GetClicks(window.onScreen);
WaitSecs(1);
KbWait %wait for subject to press button

WaitSecs(1);

end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------

function finish(window)
prefs = getPreferences();
Screen('TextSize', window.onScreen, window.fontsize);
Screen('DrawText', window.onScreen, 'You are now finished with the task.',(window.centerX-400),(window.centerY+20), 255);
Screen('DrawText', window.onScreen, 'Please let the experimenter know by using the call box.',(window.centerX-400),(window.centerY+50), 255);
Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
Screen('Flip', window.onScreen);
% GetClicks(window.onScreen);
KbWait; %wait for subject to press button
end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------

function offsets = circularArrayOffsets(n, centerX, centerY, radius, rotation)
degreeStep = 360/n;
offsets = [sind(0:degreeStep:(360-degreeStep) + rotation)'.* radius, ...
    cosd(0:degreeStep:(360-degreeStep) + rotation)'.* radius];
end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------

function rects = circularArrayRects(rect, nItems, radius, centerX, centerY)
coor = circularArrayOffsets(nItems, centerX, centerY, radius, 0) + repmat([centerX centerY], nItems, 1);
rects = [coor(:, 1)-rect(3)/2 , coor(:, 2)-rect(3)/2, coor(:, 1)+rect(3)/2, coor(:, 2)+rect(3)/2];
end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------
%% ~~~ Open the main window and get dimensions ~~~
function window = openWindow()

window.screenNumber = max(Screen('Screens'));
window.onScreen = Screen('OpenWindow', window.screenNumber, [127.5 127.5 127.5]);
Screen('BlendFunction', window.onScreen, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
[window.screenX, window.screenY] = Screen('WindowSize', window.onScreen); % check resolution
window.screenRect  = [0, 0, window.screenX, window.screenY]; % screen rect
window.centerX = window.screenX * 0.5; % center of screen in X direction
window.centerY = window.screenY * 0.5; % center of screen in Y direction
window.centerXL = floor(mean([0, window.centerX])); % center of left half of screen in X direction
window.centerXR = floor(mean([window.centerX, window.screenX])); % center of right half of screen in X direction
% targets need to be closer to the center (~6.3 degrees of arc)
window.centerXL = floor(mean([window.centerXL, window.centerX])); % center of the center left half of screen in X direction
window.centerXR = floor(mean([window.centerXR, window.centerX])); % center of the center right half of screen in X direction

% Basic drawing and screen variables.
window.black      = BlackIndex(window.onScreen); % RGB value = 0
window.white      = WhiteIndex(window.onScreen); % RGB value = 255
window.gray       = mean([window.black window.white]); % RGB value = 127.5
window.fontsize   = 26; % size of instruction text
window.bgcolor    = window.gray; % set background color

end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------
function L = orientwheelLocations(window,prefs)
L = [sind(1:360).*prefs.orientWheelRadius + window.centerX;...
    cosd(1:360).*prefs.orientWheelRadius + window.centerY];
end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------
%% *****Set-up off screen window with mask*****
function maskwin = createMasktex(window,prefs)

% Open a new window off screen
maskwin = Screen('OpenOffscreenwindow', window.onScreen, window.gray); 

% Set up alpha-blending for smooth (anti-aliased) lines in mask window
Screen('BlendFunction', maskwin, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% Draws each target stimulus to create mask
degslocs = [15:15:360,360:-15:15]; %deg of orientation for drawing mask
for ii = 1:length(degslocs)

    yadj = prefs.orientwheel(2,degslocs(ii)); %adjust lines y position 
    xadj = prefs.orientwheel(1,degslocs(ii)); %adjust lines x position 
    yadj2 = prefs.orientwheel2(2,degslocs(ii)); %for edge lines 
    xadj2 = prefs.orientwheel2(1,degslocs(ii)); %for edge lines

    % Set the coordinates (these are all relative to zero we will let
    % the drawing routine center to our monitor for us)
    % Set the new coordinates
    xCoords = [xadj -xadj2 xadj xadj2 xadj 0];
    yCoords = [yadj yadj2 yadj -yadj2 yadj 0];
    allCoords = [xCoords; yCoords];

    % Draw the orientation lines, set it to the center of our screen and
    % set good quality antialiasing
    Screen('DrawLines', maskwin, allCoords, prefs.lineWidthPix, prefs.mask_gray,...
        [window.centerXR window.centerY], 2);
    Screen('DrawLines', maskwin, allCoords, prefs.lineWidthPix, prefs.mask_gray,...
        [window.centerXL window.centerY], 2);
    
    clear xCoords yCoords allCoords
    
end
clear yadj xadj yadj2 xadj2

% Draw the circle to the screen
centeredRect = CenterRectOnPointd([0,0,14,14], window.centerXR, window.centerY); 
Screen('FillOval', maskwin, prefs.mask_gray, centeredRect');
centeredRect = CenterRectOnPointd([0,0,14,14], window.centerXL, window.centerY); 
Screen('FillOval', maskwin, prefs.mask_gray, centeredRect');

end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------
%% *****Set-up locations of taret place holders*****
function placeRect = createPlaceHolder(window,prefs)

% Make a base Rect of diameter size
baseRect = [0 0 prefs.placeRadius*2 prefs.placeRadius*2];

% Center the location of rec
placeRect(:,1) = CenterRectOnPointd(baseRect, window.centerXR, window.centerY);
placeRect(:,2) = CenterRectOnPointd(baseRect, window.centerXL, window.centerY);

end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------
%% *****Set-up off screen window with RIGHT cue*****
function CueR = createCueR(window)

% Open a new window off screen
CueR = Screen('OpenOffscreenwindow', window.onScreen, window.gray); 

% Set up alpha-blending for smooth (anti-aliased) lines in mask window
Screen('BlendFunction', CueR, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% Set the color of the triangle
rectColor = window.black;

% Cue to tell PTB that the polygon is convex (concave polygons require much
% more processing)
isConvex = 1;

right_arrow = [window.centerX+10 window.centerY-30; window.centerX-10 window.centerY-40; window.centerX-10 window.centerY-20];

% Draw the rect to the screen
Screen('FillPoly', CueR, rectColor, right_arrow, isConvex);

end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------
%% *****Set-up off screen window with LEFT cue*****
function CueL = createCueL(window)

% Open a new window off screen
CueL = Screen('OpenOffscreenwindow', window.onScreen, window.gray); 

% Set up alpha-blending for smooth (anti-aliased) lines in mask window
Screen('BlendFunction', CueL, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% Set the color of the triangle
rectColor = window.black;

% Cue to tell PTB that the polygon is convex (concave polygons require much
% more processing)
isConvex = 1;

left_arrow = [window.centerX-10 window.centerY-30; window.centerX+10 window.centerY-40; window.centerX+10 window.centerY-20];   %create the attention cues

% Draw the rect to the screen
Screen('FillPoly', CueL, rectColor, left_arrow, isConvex);
end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------
%% ~~~ Set task variables for timing and stimuli ~~~
function prefs = getPreferences

prefs.retentionIntervals = [0.50]; %time between target and color wheel (I think)

prefs.rate = 10; % constitutes a 12 Hz rhythmic presentation (12 would be a 10 Hz)
    % refresh cycles before next entrainer (1000msec /120 Hz) = 8.333 msec; 1000msec / 10 Hz  = 100 msec; 100 msec / 8.333 msec = 12 cycles
        %1000/8.33*6 == 20Hz, 1000/8.33*8 == 15Hz, 1000/8.33*10 == 12Hz, 1000/8.33*14 == 8.5Hz, 1000/8.33*30 == 4.0Hz
        %number of refreshes, so formula = 1000/(8.33*desired frequency) desired frequency = 12, 15, 20, etc. 

        
% ---------------
% Target --------
% ---------------
prefs.targ_length = 1; %in refresh cycles; each refresh is 1000msec / 120 Hz = 8.333 msec
prefs.setSizes = [1]; %number of targets presented (should always be 1)
prefs.squareSize = 20; % size of each stimulus object, in pixels
prefs.radius = 0; %how far apart the targets are from each other (if more than 1 target)
% What target will look like
prefs.targ_thresh = 50; %target RGB points darker than grey background (15)
prefs.baseCircleresp = [0,0,10,10]; %Size of target's central circle in pixels
% Used to draw the targets and response wheel
prefs.orientwheel = [sind(1:360).*30; cosd(1:360).*30]; %all possible x & y positions for center line
prefs.orientwheel2 = [cosd(1:360).*2; sind(1:360).*2]; %all possible x & y positions edge lines
prefs.lineWidthPix = 2; %Set the line width for our response stimulus


% ---------------
% Entrainers ----
% ---------------
prefs.entr_grey = 0; %points darker than background (do not want to see entrainers)
prefs.entr_length = 1 ; %refresh cycles of entrainer = refresh cycles of Draw.entrainer
prefs.entr_gap_length = prefs.rate - prefs.entr_length;
prefs.n_entrs = [1]; %number of entrainers on a trial (can be list: entrs = [6:1:8]);
                     %number of entrainers set to 1 because 0 will break the code


% ---------------
% Lags ----------
% ---------------
% number of refreshes (cycles) after last entrainer before target onset (can be list: lags = [4:2:8])
prefs.lags = [prefs.rate/2:prefs.rate/2:prefs.rate*2]; %lag*8.333 = time ms
prefs.n_lags = length(prefs.lags); %number of unique lags
prefs.lagISI = prefs.lags - 1;
prefs.lags_per_block = 105; %this is from when entrainers were used (not so important now)


% ---------------
% Mask ----------
% ---------------
prefs.SOA = 6; %refresh target onset to mask onset (50 ms optimal/8.3333 msec = 6 cycles)
prefs.maskISI = prefs.SOA - 1; %target OFFSET to mask onset
prefs.maskwidth = 20;
prefs.mask_thresh = 50; %points darker than background
prefs.mask_length = 1; %refresh cycles of mask


% ---------------
% Fixation ------
% ---------------
prefs.fixation_length = 84; %700ms
prefs.fixationSize = 4; %size of dot


% ---------------------
% Orientation Wheel ---
% ---------------------
% To estimate orientation based on location of mouse (not actually drawn)
prefs.orientWheelRadius = 35; %size of orientation wheel (not drawn)


% ------------
% Cue --------
% ------------
%in refresh cycles; each refresh is 1000msec / 120 Hz = 8.333 msec
prefs.cue_length = 120; %1000 ms
prefs.cue_size = 20;
prefs.cue_line = 5; %thickness of arrow line
prefs.postcue_length = 24; %200 ms
prefs.p_validity = 1; %proportion of informative cues (all should be informative)

prefs.preblank_length = 24; %200ms


% Other Variables
prefs.trigger_size = [0 0 1 1];
prefs.p_catchtrials = 0; %what proportion of trials will be catch trials
prefs.tilesize = 5; % how big the coloured squares are

prefs.timeLimit = 5; %how long wait for yes/no response (s)

prefs.placeRadius = prefs.orientWheelRadius + 5; %size of target place holder



% -------------------
% Trials & Blocks ---
% -------------------
% Randomize trial order of full factorial design order.
% prefs.fullFactorialDesign = fullfact([length(prefs.setSizes), ...
%     prefs.lags_per_block, ...
%     prefs.n_lags]);
% prefs.order = Shuffle(1:length(prefs.fullFactorialDesign));

% All possible orientations for targets to be at
prefs.degslocs = [15:15:360]; 

% Number of blocks
prefs.nBlocks = 3;
% prefs.nBlocks = 2;

% Trials per block (multiple of # of target orientations)
prefs.trials_per_block = length(prefs.degslocs)*2;

% Total number of experimental trials (384)
prefs.nTrials = prefs.nBlocks*prefs.trials_per_block;
% prefs.nTrials = 12; %for testing

% Random order of orientations on each trial
prefs.selectorientTrial = randi(length(prefs.degslocs),[(prefs.nTrials),1]);



end %end prefs function

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------




































