% OrientationWheelTask_Exo_v2() 
% Runs a color working memory task a la Zhang & Luck (2008). The task 
% requires memory for the color of briefly presented squares. Participants 
% then report the color of a single probed square using a contnuous report 
% task.
%
% 
% Fixation appears for one entrainer length instead of a blank screen.
% 
% -------------------------------------------------------------------------
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% _________________________________________________________________________
% 
% /////////////////////////////////////////////////////////////////////////
%
% ==== TRIGGERS ====
%
% Start eye tracker: 199
% 
% -- Targets Present --
% Trial start (fixation): 1
% Entrainers: 61-68*
% Cue: 
%     left:  2
%     right: 3
%     both:  9
% Target:
%     left:  30
%     right: 20
% Mask: 90
% Response screen target: 5  
%   Response target: 80 
% Response screen Y/N: 25  
%   Response Yes: 180 
%   Response No: 185 
%   Response invalid: 150
%
% -- Targets Not Present --
% Trial start (fixation): 11
% Entrainers: 61-68
% Cue: 
%     left:  12
%     right: 13
%     both:  19
% Target: 40
% Mask: 95
% Response screen target: 15
%   Response target: 70
% Response screen Y/N: 35  
%   Response Yes: 170 
%   Response No: 175 
%   Response invalid: 140
% 
% /////////////////////////////////////////////////////////////////////////


function OrientationWheelTask_Exo_v3()

ccc

try    
    
prepareEnvironment;

% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------    
% Input participant's number, name, date of testing, preferences
part_num = input('Participant Number:','s');
Filename = ['M:\Experiments\OrientWheel_Exo\Orient_Data\' part_num '_Orient_Exo.mat'];

% Do you already have a color value or use default? 
color_yn = input('Do you want to input the target color? [y or n]:','s');
if strncmpi(color_yn,'y',1)
    % If not doing staircasing, input target color (0 black to 128 background grey)
    xtarget_gray = input('Target color [0 to 128]:','s');
    xtarget_gray = str2num(xtarget_gray); %convert input to number
end 

% Is eyetracking going to be used? 
% eyetrack_yn = input('Are you using the eye tracker? [y or n]:','s');
eyetrack_yn = 'n'; %eyetracker not working so always no for now

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% /////////////////////////////////////////////////////////////////////////
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

window = openWindow(); %get info for psychtoolbox
prefs = getPreferences(); %get task variable info

% /////////////////////////////////////////////////////////////////////////
% /////////////////////////////////////////////////////////////////////////
%% Instructions for experimenter to start eye tracking
if strncmpi(eyetrack_yn,'y',1)
    Screen('FillRect',window.onScreen,window.gray);
    DrawFormattedText(window.onScreen,'Calibrated? Start EEG Recording and Press the SPACE BAR','center','center',[]);
    Screen('FillRect',window.onScreen,Vpixx2Vamp(0),prefs.trigger_size);
    Screen('Flip',window.onScreen)

    KbWait;

    Screen('FillRect',window.onScreen,window.gray);
    Screen('FillRect',window.onScreen,Vpixx2Vamp(199),prefs.trigger_size);
    Screen('Flip',window.onScreen);
    % system('C:\Users\user\Downloads\CoreSDK\CoreSDK\samples\Streams\Interaction_Streams_101\bin\Debug\Interaction_Streams_101.exe &');
    system('M:\Experiments\micb_eyetrack\CoreSDK\CoreSDK\samples\Streams\Interaction_Streams_101\bin\Debug\Interaction_Streams_101.exe &');

    WaitSecs(2);
end

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% /////////////////////////////////////////////////////////////////////////
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% counter of trials per block
nblock = prefs.trials_per_block + 20; %including the first 20 practice trials
cnt_blck = 0; %counter of blocks

% total number of trials including practice trials
ntotaltrial = prefs.nTrials + 20; %add 20 practice trials

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
% Center the target oval & cues on the centre of the screen
centeredRectresp = CenterRectOnPointd(prefs.baseCircleresp,window.centerX,window.centerY);
centeredRectrespR = CenterRectOnPointd(prefs.baseCircleresp,window.centerXR,window.centerY); 
centeredRectrespL = CenterRectOnPointd(prefs.baseCircleresp,window.centerXL,window.centerY);

% -------------------------------------------------------------------------
% Set colors of stimuli
% prefs.mask_gray = window.gray - prefs.mask_thresh; %mask color (77.5)
prefs.mask_gray = 0; %default is black
% Use default stimuli colors if not specified above
if strncmpi(color_yn,'y',1)
%    prefs.targ_gray = window.gray - prefs.targ_thresh; %default target color (77.5)
    prefs.targ_gray = xtarget_gray;
else
    prefs.targ_gray = 0; %default is black
end

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
CueN = createCueN(window); %create offscreen window with non-informative cue

% -------------------------------------------------------------------------    
% Put up instructions and wait for keypress.
instruct(window,prefs);

% -------------------------------------------------------------------------    
% Location of orientations on screen
orientWheelLocations = orientwheelLocations(window,prefs);

% -------------------------------------------------------------------------    

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% Set up the random mix of lags and target position and cue for each block

%Now to find a lag you pick the next random index between 1:n_lags from i_lags
%Then you find that index in all_lags, it tells you where to look in lags
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
    
    
%set up the catch trials on every n-lagsth trial
p = 1/prefs.p_catchtrials;
q = 1/prefs.p_catchtrials;
present = [1];
for i_pres = 2:ntotaltrial
    if i_pres == p
        present = [present 0];
        p = p + q;
    else
        present = [present 1];
    end
end
rand_pres = randperm(ntotaltrial);
present = present(rand_pres);
clear rand_pres    


% Pre-allocate some variables
allCoords_targ = cell(1,ntotaltrial);
trialOrient = cell(1,ntotaltrial);   

    
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% START TASK
    for trialIndex = 1:ntotaltrial
        
        lag = prefs.lags(all_lags(trialIndex))/5;
        
        % Determine how many items there are on this trial and the duration (this is left over code from original task)
        nItems = 1; %should always be 1 because there will only be 1 target
        retentionInterval = prefs.retentionIntervals; %(prefs.fullFactorialDesign(prefs.order(trialIndex), 2));
        
        % Pick an item to test.
%         itemToTest(trialIndex) = randsample(1:nItems); %legacy code that no longer works
        itemToTest(trialIndex) = nItems; %nItems should always be 1 because there will only be 1 target

        % Pick the orientation for this trial.
        trialOrient{trialIndex} = prefs.degslocs(prefs.selectorientTrial(trialIndex));
        
        
        if present(trialIndex) == 1
            pause(window);
            %Present Fixation
            Screen('FillRect', window.onScreen, window.gray);
            Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255);
            Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect); %placeholder
            Screen('FillRect',window.onScreen, Vpixx2Vamp(1), prefs.trigger_size);
            t_fixate_onset = Screen('Flip', window.onScreen);
            
        elseif present(trialIndex) == 0
            pause(window);
            %Present Fixation
            Screen('FillRect', window.onScreen, window.gray);
            Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255);
            Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect); %placeholder
            Screen('FillRect',window.onScreen, Vpixx2Vamp(11), prefs.trigger_size);
            t_fixate_onset = Screen('Flip', window.onScreen);
        end
        
        % Interval (refresh screen)
        Screen('FillOval', window.onScreen, window.gray, rects{nItems});
        Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255); %keep fixation till target
        Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect); %placeholder
        Screen('FillRect',window.onScreen,Vpixx2Vamp(0),prefs.trigger_size);
        Screen('Flip', window.onScreen);
        
% =========================================================================            
        %% Cue
        if present(trialIndex) == 1 %two options depending on whether the target is present or absent
            
            if valid(trialIndex) == 0 %non-informative cue
                % present the non-informative cue
                Screen('DrawTexture', window.onScreen, CueN, [], [], [], 0);
                Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255);%draw fixation
                Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect); %placeholder
                Screen('FillRect',window.onScreen,Vpixx2Vamp(9),prefs.trigger_size);
                
            else %informative cue
                if position(trialIndex) == 0 
                    % present the Left cue
                    Screen('DrawTexture', window.onScreen, CueL, [], [], [], 0);
                    Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255);%draw fixation
                    Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect); %placeholder
                    Screen('FillRect',window.onScreen,Vpixx2Vamp(2),prefs.trigger_size);
                end
                if position(trialIndex) == 1
                    % present the Right cue
                    Screen('DrawTexture', window.onScreen, CueR, [], [], [], 0);
                    Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255);%draw fixation
                    Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect); %placeholder
                    Screen('FillRect',window.onScreen,Vpixx2Vamp(3),prefs.trigger_size);
                end
            end
            
        elseif present(trialIndex) == 0 %no target
            
            if valid(trialIndex) == 0 %non-informative cue
                % present the non-informative cue
                Screen('DrawTexture', window.onScreen, CueN, [], [], [], 0);
                Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255);%draw fixation
                Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect); %placeholder
                Screen('FillRect',window.onScreen,Vpixx2Vamp(19),prefs.trigger_size);
                
            else %informative cue
                if position(trialIndex) == 0 
                    % present the Left cue
                    Screen('DrawTexture', window.onScreen, CueL, [], [], [], 0);
                    Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255);%draw fixation
                    Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect); %placeholder
                    Screen('FillRect',window.onScreen,Vpixx2Vamp(12),prefs.trigger_size);
                end
                if position(trialIndex) == 1
                    % present the Right cue
                    Screen('DrawTexture', window.onScreen, CueR, [], [], [], 0);
                    Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255);%draw fixation
                    Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect); %placeholder
                    Screen('FillRect',window.onScreen,Vpixx2Vamp(13),prefs.trigger_size);
                end
            end
        end
        tcue_onset = Screen('Flip',window.onScreen,t_fixate_onset + prefs.fixation_length*refresh - slack);
        prefs.fixation_time(trialIndex) = tcue_onset - t_fixate_onset;
                  
        
        
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
            tentr_onset = Screen(window.onScreen, 'Flip',tblank_onset + prefs.preblank_length*refresh - slack);
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
% =========================================================================  
        if present(trialIndex) == 1 %two options depending on whether the target is present or absent
% =========================================================================  
            %% Draw the target stimulus
            %same for left and right targets
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
                Screen('DrawLines', window.onScreen, allCoords_targ{trialIndex}, prefs.lineWidthPix, prefs.targ_gray, [window.centerXR window.centerY], 2);
                % Draw the central circle to the screen
                Screen('FillOval', window.onScreen, prefs.targ_gray, centeredRectrespR');
                % Draw fixation
                Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255);
                Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect); %placeholder
                % Draw trigger
                Screen('FillRect',window.onScreen,Vpixx2Vamp(20),prefs.trigger_size);
% =========================================================================                 
            else
            %% Left Target
                % Draw the orientation lines, set it to the center of our screen and set good quality antialiasing
                Screen('DrawLines', window.onScreen, allCoords_targ{trialIndex}, prefs.lineWidthPix, prefs.targ_gray, [window.centerXL window.centerY], 2);
                % Draw the central circle to the screen
                Screen('FillOval', window.onScreen, prefs.targ_gray, centeredRectrespL');
                % Draw fixation
                Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255);
                Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect); %placeholder
                % Draw trigger
                Screen('FillRect',window.onScreen,Vpixx2Vamp(30),prefs.trigger_size);
            end       
% =========================================================================            
        else
            %% present the Missing Target
            Screen('FillOval', window.onScreen, window.gray, rects{nItems});
            % Draw fixation
            Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255);
            Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect); %placeholder
            % Draw trigger
            Screen('FillRect',window.onScreen,Vpixx2Vamp(40),prefs.trigger_size);            
        end
% =========================================================================          
        % Post the stimulus
        ttarget_onset = Screen('Flip', window.onScreen,tlag_onset + prefs.lagISI(all_lags(trialIndex))*refresh - slack);
% =========================================================================  
        %% blank Inter stimulus interval
        Screen('FillOval', window.onScreen, window.gray, rects{nItems});
        %draw fixation
        Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255);
        Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect); %placeholder
        Screen('FillRect',window.onScreen,Vpixx2Vamp(0),prefs.trigger_size);
        tISI_onset = Screen('Flip', window.onScreen, ttarget_onset + prefs.targ_length*refresh - slack);
        
% =========================================================================         
        %% Mask
        % Draw the texture 
        if present(trialIndex) == 1
            Screen('DrawTexture', window.onScreen, maskwin, [], [], [], 0);
            %draw fixation
            Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255);
            Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect); %placeholder
            Screen('FillRect',window.onScreen,Vpixx2Vamp(90),prefs.trigger_size);
        elseif present(trialIndex) == 0
            Screen('DrawTexture', window.onScreen, maskwin, [], [], [], 0);
            %draw fixation
            Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255);
            Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect); %placeholder
            Screen('FillRect',window.onScreen,Vpixx2Vamp(95),prefs.trigger_size);
        end
        
        t_maskonset = Screen('Flip', window.onScreen, tISI_onset + prefs.maskISI*refresh - slack);
        
        % Interval
        Screen('FillOval', window.onScreen, window.gray, rects{nItems});
        Screen('DrawDots', window.onScreen, [window.centerX, window.centerY], prefs.fixationSize, 255); %Draw fixation
        Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect); %placeholder
        Screen('FillRect',window.onScreen,Vpixx2Vamp(0),prefs.trigger_size);
        Screen('Flip', window.onScreen, t_maskonset + prefs.mask_length*refresh - slack);
        
        % Retention interval
        WaitSecs(retentionInterval);
        
% ========================================================================= 
        %% Response Target
        
        % Choose a orientation to test, then display the response screen.
        data.presentedOrientRads(trialIndex) = deg2rad(trialOrient{trialIndex});
        data.presentedOrientDegrees(trialIndex) = trialOrient{trialIndex};
        
        if present(trialIndex) == 1
            Screen('FillRect',window.onScreen, Vpixx2Vamp(5), prefs.trigger_size);
            Screen('Flip', window.onScreen);
        else
            Screen('FillRect',window.onScreen, Vpixx2Vamp(15), prefs.trigger_size);
            Screen('Flip', window.onScreen);
        end
        
        %------------------------------------------------------------------
        % Set-up response stimuli.
        %randomize starting location of lines if mouse is not moved 
        startdeg = prefs.degslocs(randi(length(prefs.degslocs),1));
        centeredRectmouse = CenterRectOnPointd([prefs.orientwheel(1,startdeg),prefs.orientwheel(2,startdeg),...
            20,20],window.centerX-6,window.centerY-6); %-6 to not overlap target
        SetMouse(centeredRectmouse(1),centeredRectmouse(2)); %mouse start location
        ShowCursor('Arrow');
        clear centeredRectmouse
        %------------------------------------------------------------------
        % If mouse button is already down, wait for release.
        GetMouse(window.onScreen);
        buttons = 0;
        while any(buttons)
            [x, y, buttons] = GetMouse(window.onScreen);
        end
        %------------------------------------------------------------------
        everMovedFromCenter = false;
        while ~any(buttons) % keep track of mouse location if moved
            
            [x,y,buttons] = GetMouse(window.onScreen);
            [minDistance, minDistanceIndex] = min(sqrt((orientWheelLocations(1, :) - x).^2 + (orientWheelLocations(2, :) - y).^2));
            
            if(minDistance < 3) %make sure subject moved mouse
                everMovedFromCenter = true;
            end
            
            if(everMovedFromCenter) %as the mouse is moved, update coords of response stimulus to match
                % new coordinates of start and end of lines based on mouse position
                yloc = prefs.orientwheel(2,minDistanceIndex); %adjust lines y position 
                xloc = prefs.orientwheel(1,minDistanceIndex); %adjust lines x position 
                yloc2 = prefs.orientwheel2(2,minDistanceIndex); %for edge lines 
                xloc2 = prefs.orientwheel2(1,minDistanceIndex); %for edge lines 
            else
                %starting location of lines if mouse is not moved - randomized 
                yloc = prefs.orientwheel(2,startdeg); %adjust lines y position
                xloc = prefs.orientwheel(1,startdeg); %adjust lines x position 
                yloc2 = prefs.orientwheel2(2,startdeg); %for edge lines 
                xloc2 = prefs.orientwheel2(1,startdeg); %for edge lines 
%                 yloc = -30;
%                 xloc = 0;
%                 yloc2 = 0;
%                 xloc2 = 2;
            end
            
            % Set the new coordinates
            xCoords = [xloc -xloc2 xloc xloc2 xloc 0];
            yCoords = [yloc yloc2 yloc -yloc2 yloc 0];
            allCoords = [xCoords; yCoords];

            % Draw the orient lines, set it to the center of our screen and
            % set good quality antialiasing
            Screen('DrawLines', window.onScreen, allCoords, prefs.lineWidthPix, prefs.mask_gray, [window.centerX window.centerY], 2);
            % Draw the circle part of the stimulus to the screen
            Screen('FillOval', window.onScreen, prefs.mask_gray, centeredRectresp');
            
            Screen('FillRect',window.onScreen,Vpixx2Vamp(0),prefs.trigger_size);
            Screen('Flip', window.onScreen);
            
            clear xCoords yCoords allCoords
            
        end
        %------------------------------------------------------------------
        
        if present(trialIndex) == 1
            Screen('FillRect',window.onScreen,Vpixx2Vamp(80),prefs.trigger_size);
            Screen('Flip', window.onScreen);

        elseif present(trialIndex) == 0
            Screen('FillRect',window.onScreen,Vpixx2Vamp(70),prefs.trigger_size);
            Screen('Flip', window.onScreen);
        end
        
        % calculate response location
        [respDistance, respDistanceIndex] = min(sqrt((orientWheelLocations(1,:) - x).^2 + (orientWheelLocations(2,:) - y).^2));
        data.reportedOrientRads(trialIndex) = deg2rad(respDistanceIndex);
        data.reportedOrientDegrees(trialIndex) = respDistanceIndex;
        data.reportedRespDistance(trialIndex) = respDistance;
        
        % put trial info into data structure
        data.lags(trialIndex) = lag;
        data.present(trialIndex) = present(trialIndex); %target present or absent
        
        
        HideCursor
        
% =========================================================================     
% ========================================================================= 
        %% Response Yes/No
        
%         SetMouse(window.centerX, window.centerY); %mouse at center
%         ShowCursor('Arrow');
%         
%         % Create response screen
%         answer = BinaryQuestion_OrientWheel(window.onScreen,255,window.black,...
%             'Did you see the target?','Yes','No', prefs.trigger_size, present(trialIndex));
%         
%         %------------------------------------------------------------------
%         %Responded Yes
%         if answer == 1
%             data.seen_resp(trialIndex) = 1; %record response
%             if present(trialIndex) == 1
%                 Screen('FillRect',window.onScreen,Vpixx2Vamp(180),prefs.trigger_size);
%                 Screen('Flip', window.onScreen);
%             elseif present(trialIndex) == 0
%                 Screen('FillRect',window.onScreen,Vpixx2Vamp(170),prefs.trigger_size);
%                 Screen('Flip', window.onScreen);
%             end
%         % -------------------------------------------------------------------------------------------   
%         %Responded No
%         elseif answer == 2
%             data.seen_resp(trialIndex) = 0; %record response
%             if present(trialIndex) == 1
%                Screen('FillRect',window.onScreen,Vpixx2Vamp(185),prefs.trigger_size);
%                Screen('Flip', window.onScreen);
%             elseif present(trialIndex) == 0
%                Screen('FillRect',window.onScreen,Vpixx2Vamp(175),prefs.trigger_size);
%                Screen('Flip', window.onScreen);
%             end  
%         % ------------------------------------------------------------------------------------------  
%         %Did not click on text
%         else
%             data.seen_resp(trialIndex) = -2; %record response
%             if present(trialIndex) == 1
%                Screen('FillRect',window.onScreen,Vpixx2Vamp(150),prefs.trigger_size);
%                Screen('Flip', window.onScreen);
%             elseif present(trialIndex) == 0
%                Screen('FillRect',window.onScreen,Vpixx2Vamp(140),prefs.trigger_size);
%                Screen('Flip', window.onScreen);
%             end   
%         end     
%                 
%         HideCursor
%         
% % =========================================================================

        
% =========================================================================        
        if trialIndex == ntotaltrial
            finish(window);
        elseif trialIndex == nblock
            cnt_blck = cnt_blck + 1; %count # of blocks
            rest(window,cnt_blck);
            nblock = nblock + prefs.trials_per_block;
        elseif trialIndex == 20 % first 20 trials are practice trials
            practice(window);
        end
% =========================================================================   
    end %end trials
    
% =========================================================================     
    %% Preliminary analysis of results.
    %errors
    data.errorDegrees = (180/pi) .* (angle(exp(1i*data.reportedOrientRads)./exp(1i*data.presentedOrientRads)));
    data.errorRads = deg2rad(data.errorDegrees);
    data.error_differenceRads = data.reportedOrientRads - data.presentedOrientRads;
    data.error_differenceDegrees = data.reportedOrientDegrees - data.presentedOrientDegrees;
    
    %lags
    prefs.all_lags = all_lags;
    
    %cues
    data.position = position;
    data.valid = valid;
    
    save(Filename, 'data', 'prefs');
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
Screen('Preference', 'SkipSyncTests', 1); %choose 1 for test on an LCD

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
%% Orientation wheel task instructions
function instruct(window,prefs)

%Screen 1
% prefs = getPreferences();
Screen('TextSize', window.onScreen, window.fontsize);
Screen('DrawText', window.onScreen, 'The next few screens will briefly explain the task.',(window.centerX-400),(window.centerY+20), 255);
Screen('DrawText', window.onScreen, 'After the instructions, you will perform several practice trials followed by the experiment.',...
    (window.centerX-400),(window.centerY+50), 255); %line every +30
Screen('DrawText', window.onScreen, 'On each screen click the mouse to continue.',(window.centerX-400),(window.centerY+80), 255); %line every +25
Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
Screen('Flip', window.onScreen);
GetClicks(window.onScreen);


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
GetClicks(window.onScreen);


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
GetClicks(window.onScreen);

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
GetClicks(window.onScreen);


%Screen 2.4
Screen('TextSize', window.onScreen, window.fontsize); 
CueN = createCueN(window); %create offscreen window with non-informative cue
Screen('DrawTexture', window.onScreen, CueN, [], [], [], 0);
Screen('DrawText', window.onScreen, 'Sometimes the arrow will be pointing both left and right.',...
    (window.centerX-400),(window.centerY+50), 255);
Screen('DrawText', window.onScreen, 'This means the target stimulus can appear to either the left OR right.',...
    (window.centerX-400),(window.centerY+80), 255); %line every +30
%place holder
placeRect = createPlaceHolder(window,prefs);
Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect);
%put stimulus and text on screen
Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
Screen('Flip', window.onScreen);
GetClicks(window.onScreen);


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
Screen('DrawLines', window.onScreen, allCoords_targ, prefs.lineWidthPix, prefs.mask_gray, [window.centerXR window.centerY], 2);
centeredRectresp = CenterRectOnPointd([0,0,10,10], window.centerXR, window.centerY); %location of center for oval
Screen('FillOval', window.onScreen, prefs.mask_gray, centeredRectresp');
%left target
Screen('DrawLines', window.onScreen, allCoords_targ, prefs.lineWidthPix, prefs.mask_gray, [window.centerXL window.centerY], 2);
centeredRectresp = CenterRectOnPointd([0,0,10,10], window.centerXL, window.centerY); %location of center for oval
Screen('FillOval', window.onScreen, prefs.mask_gray, centeredRectresp');
%place holder
placeRect = createPlaceHolder(window,prefs);
Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect);
%put stimulus and text on screen
Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
Screen('Flip', window.onScreen);
GetClicks(window.onScreen);



%Screen 4
Screen('TextSize', window.onScreen, window.fontsize);
maskwin = createMasktex(window,prefs); %create offscreen window with mask
Screen('DrawTexture', window.onScreen, maskwin, [], [], [], 0);
Screen('DrawText', window.onScreen, 'A star-like shape will then quickly appear and then disappear.',...
    (window.centerX-400),(window.centerY+50), 255);
Screen('DrawText', window.onScreen, 'Remember that you want to detect the direction of the target, not the star.',...
    (window.centerX-400),(window.centerY+80), 255); %line every +30
%place holder
placeRect = createPlaceHolder(window,prefs);
Screen('FrameOval', window.onScreen, prefs.mask_gray, placeRect);
Screen('FillRect',window.onScreen,Vpixx2Vamp(0),prefs.trigger_size);
Screen('Flip', window.onScreen);
GetClicks(window.onScreen);


%Screen 5
%draw target stimulus
xCoords = [-21.21 -1.41 -21.21 1.41 -21.21 0];
yCoords = [21.21 -1.41 21.21 1.41 21.21 0];
allCoords_targ = [xCoords; yCoords];
Screen('DrawLines', window.onScreen, allCoords_targ, prefs.lineWidthPix, prefs.mask_gray, [window.centerX window.centerY], 2);
centeredRectresp = CenterRectOnPointd([0,0,10,10], window.centerX, window.centerY); 
Screen('FillOval', window.onScreen, prefs.mask_gray, centeredRectresp');
Screen('TextSize', window.onScreen, window.fontsize);
Screen('DrawText', window.onScreen, 'Another target will then appear on screen.',(window.centerX-880),(window.centerY-50), 255);
Screen('DrawText', window.onScreen, 'Using the mouse, move the new target until',(window.centerX-880),(window.centerY), 255);
Screen('DrawText', window.onScreen, 'it is pointing in the exact same direction',(window.centerX-880),(window.centerY+30), 255);
Screen('DrawText', window.onScreen, 'of the original target that appeared in the circle.',(window.centerX-880),(window.centerY+60), 255);
Screen('DrawText', window.onScreen, 'If you are unsure of the direction of the target,',(window.centerX-880),(window.centerY+110), 255);
Screen('DrawText', window.onScreen, 'or if you did not see a target, please provide',(window.centerX-880),(window.centerY+140), 255);
Screen('DrawText', window.onScreen, 'your best guess.',(window.centerX-880),(window.centerY+170), 255);
Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
Screen('Flip', window.onScreen);
GetClicks(window.onScreen);


%Screen 6
Screen('TextSize', window.onScreen, window.fontsize);
Screen('DrawText', window.onScreen, 'You will now complete several practice trials so that you become more familiar with the task.',(window.centerX-400),(window.centerY+20), 255);
Screen('DrawText', window.onScreen, 'Remember to remain focused and fixated on the white central dot that will appear.',(window.centerX-400),(window.centerY+50), 255);
Screen('DrawText', window.onScreen, 'Click the left mouse button when you are ready to begin the practice trials.',(window.centerX-400),(window.centerY+80), 255);
Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
Screen('Flip', window.onScreen);
GetClicks(window.onScreen);
end

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

function practice(window)
prefs = getPreferences();
Screen('TextSize', window.onScreen, window.fontsize);
Screen('DrawText', window.onScreen, 'You have now finished the set of practice trials.',(window.centerX-400),(window.centerY+20), 255);
Screen('DrawText', window.onScreen, 'Please let the experimenter know by using the call box.',(window.centerX-400),(window.centerY+50), 255);
Screen('DrawText', window.onScreen, 'Do NOT begin the experiment.',(window.centerX-400),(window.centerY+80), 255);
Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
Screen('Flip', window.onScreen);
GetClicks(window.onScreen);
end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------

function rest(window,cnt_blck)
prefs = getPreferences();
tot_blcks = num2str(prefs.nBlocks);
rest_str = ['You have completed ' num2str(cnt_blck) ' out of ' tot_blcks ' blocks!'];
Screen('TextSize', window.onScreen, window.fontsize);
Screen('DrawText', window.onScreen, rest_str,...
    (window.centerX-400),(window.centerY+20), 255);
Screen('DrawText', window.onScreen, 'Feel free to take a break at this time. Click the mouse when you are ready to continue.',...
    (window.centerX-400),(window.centerY+50), 255); %line every +30
Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
Screen('Flip', window.onScreen);
GetClicks(window.onScreen);
end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------

function finish(window)
prefs = getPreferences();
Screen('TextSize', window.onScreen, window.fontsize);
Screen('DrawText', window.onScreen, 'You are now finished the experiment. Thank you for your time.',...
    (window.centerX-400),(window.centerY+20), 255);
Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
Screen('Flip', window.onScreen);
GetClicks(window.onScreen);
end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------

% function drawFixation(window, fixationX, fixationY, fixationSize)
% prefs = getPreferences();
% Screen('DrawDots', window.onScreen, [fixationX, fixationY], fixationSize, 255);
% Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
% end

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
%% Draw the color wheel
function drawOrientWheel(window)
prefs = getPreferences();
orientWheelLocations = [sind(1:360).*prefs.orientWheelRadius + window.centerX;...
    cosd(1:360).*prefs.orientWheelRadius + window.centerY];
% colorWheelSizes = 60;

% Now draws "rotated" color wheel
% Screen('DrawDots', window.onScreen, orientWheelLocations, colorWheelSizes, trial_color_wheel', [], 1);
Screen('FillRect',window.onScreen, Vpixx2Vamp(0), prefs.trigger_size);
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

prefs = getPreferences();

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

Screen('FillRect',CueR, window.black, prefs.trigger_size);

end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------
%% *****Set-up off screen window with LEFT cue*****
function CueL = createCueL(window)

prefs = getPreferences();

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

Screen('FillRect',CueL, window.black, prefs.trigger_size);
end

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------
%% *****Set-up off screen window with NON-INFORMATIVE cue*****
function CueN = createCueN(window)

prefs = getPreferences();

% Open a new window off screen
CueN = Screen('OpenOffscreenwindow', window.onScreen, window.gray); 

% Set up alpha-blending for smooth (anti-aliased) lines in mask window
Screen('BlendFunction', CueN, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% Set the color of the triangle
rectColor = window.black;

% Cue to tell PTB that the polygon is convex (concave polygons require much
% more processing)
isConvex = 1;

right_arrow = [window.centerX+10 window.centerY-30; window.centerX-10 window.centerY-40; window.centerX-10 window.centerY-20];
left_arrow = [window.centerX-10 window.centerY-30; window.centerX+10 window.centerY-40; window.centerX+10 window.centerY-20];   %create the attention cues

% Draw the rect to the screen
Screen('FillPoly', CueN, rectColor, right_arrow, isConvex);
Screen('FillPoly', CueN, rectColor, left_arrow, isConvex);

Screen('FillRect',CueN, window.black, prefs.trigger_size);
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
prefs.p_validity = 0.67; %proportion of informative cues (33 non-informative cues)

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
prefs.nBlocks = 8;
% prefs.nBlocks = 2;

% Trials per block (multiple of # of target orientations)
prefs.trials_per_block = length(prefs.degslocs)*2;

% Total number of experimental trials (384)
prefs.nTrials = prefs.nBlocks*prefs.trials_per_block;
% prefs.nTrials = 12; %for testing

% Random order of orientations on each trial
prefs.selectorientTrial = randi(length(prefs.degslocs),[(prefs.nTrials+20),1]);



end %end prefs function

% -------------------------------------------------------------------------
% #########################################################################
% -------------------------------------------------------------------------












































