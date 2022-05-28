%% ANALYSIS WRAPPER

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



%clear and close everything
ccc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        >>>>> Description of the saved dataset settings <<<<<
% 
% -> byCues_v1: FFT; winsize is 512, no ERSP baseline. Epoched to cues;
%    postocularthresh = [-500 500]. Epoch limit [-1.2 2.2]. ERP baseline 
%    [-200 0]. Filter on [0.1 50]. TF freq range [1 40]. Timesout = 300. 
%    Padratio = 4.
% 
%            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% -> byCues_v2: FFT; winsize is 512, no ERSP baseline. Epoched to cues;
%    postocularthresh = [-500 500]. Epoch limit [-1.2 2.2]. ERP baseline 
%    [-200 0]. Filter on [0.1 50]. TF freq range [1 40]. Timesout = 300. 
%    Padratio = 4. Rejected trials w/horizontal eye movements, thresh = 20,
%    timewin = [1240 1950]
% 
%            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% -> byCues_v3: FFT; winsize is 512, no ERSP baseline. Epoched to cues;
%    postocularthresh = [-500 500]. Epoch limit [-1.2 2.2]. ERP baseline 
%    [-200 0]. Filter on [0.1 50]. TF freq range [1 40]. Timesout = 300. 
%    Padratio = 4. 
%    Rejected trials w/horizontal eye movements, thresh = 25, timewin = 
%    [-200 1950], winSize = 100, winStep = 10
%    No visual inspection to reject trials.
% 
% 
%            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
%            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% -> byTargets_v1: FFT; winsize is 512, no ERSP baseline. Epoched to targets;
%    postocularthresh = [-500 500]. Epoch limit [-2.4 1.4]. ERP baseline 
%    [-1600 -1400]. Filter on [0.1 50]. Freq range [1 40]. Timesout = 300. 
%    Padratio = 4.
% 
%            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% -> byTargets_v2: FFT; winsize is 512, no ERSP baseline. Epoched to targets;
%    postocularthresh = [-500 500]. Epoch limit [-2.4 1.4]. ERP baseline 
%    [-1600 -1400]. Filter on [0.1 50]. Freq range [1 40]. Timesout = 300. 
%    Padratio = 4. Rejected trials w/horizontal eye movements, thresh = 20,
%    timewin = [0 600]
% 
%            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% -> byTargets_v3: FFT; winsize is 512, no ERSP baseline. Epoched to targets;
%    postocularthresh = [-500 500]. Epoch limit [-2.4 1.4]. ERP baseline 
%    [-1600 -1400]. Filter on [0.1 50]. Freq range [1 40]. Timesout = 400. 
%    Padratio = 4. 
%    Rejected trials w/horizontal eye movements, thresh = 25, timewin = 
%    [-1600 600], winSize = 100, winStep = 10
%    No visual inspection to reject trials.
% 
%            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% -> byTargets_v4: FFT; winsize is 512, no ERSP baseline. Epoched to targets;
%    postocularthresh = [-500 500]. Epoch limit [-2.4 1.4]. ERP baseline 
%    [-1600 -1400]. Filter on [0.1 50]. Freq range [1 40]. Timesout = 400. 
%    Padratio = 4. 
%    Rejected trials w/horizontal eye movements, thresh = 25, timewin = 
%    [-1600 600], winSize = 100, winStep = 10
%    No visual inspection to reject trials.
% 
%            ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
%    EXAMPLE FROM ORIGINAL ORIENTATION TASK
% -> filt_byTargets_v4: winsize is 256, no ERSP baseline, epoched to targets; 
%    Epoch limit [-1.5 1.5]. ERP baseline [-200 0]. Filter on [0.1 50]. 
%    Cycles [2 0.6] (2 cycles at lowest freq & 12 at highest). Freq range [2 40]. 
%    Freq increase in steps of 1.027 Hz. Timesout = 300. Padratio = 4. 
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load data
exp.name = 'OrientWheel_Exo';
exp.conds = ''; %conds is usually used for comparing the same type of trials
                %under different conditions (e.g., stimulation vs sham)
                
% exp.pathname = 'M:\Data\OrientWheel_Exo\'; %path of data on server
% exp.pathname = 'C:\Users\ssshe\Documents\MathLab\Data\OrientWheel_Exo\'; %path of data (personal computer)
exp.pathname = 'D:\MathLab\Data\OrientWheel_Exo\'; %data location (personal data storage)

% exp.setname = {'byTargets_v3'}; % name each epoched set
exp.setname = {'byCues_v3'}; % name each epoched set
% note: the meaning of set here is to identify the type of analysis done.
%       set is usually used to identify different trial types (e.g., standards
%       vs targets) within the same experimental condition.

%% List of participants' ids 
% exclude 20 because trials unknown - all
%
% *byCues_v1 & *byTargets_v1:
% exclude 12 because data failed to fit model
% exclude 6 because extreme outlier in BEH data
% 
% exp.participants = {'004';'005';'006';'007';'008';'009';'010';'011';...
%     '012';'013';'014';'015';'017';'018';'019';'021';'022';'023';'024';...
%     '025';'026';'027';'028';'029';'030';'031';'032';'033';'034';'035';...
%     '036';'037';'038';'039';'040'};

% *byTargets_v2 & *byCues_v2:
% 28 has 38% of trials removed due to horizontal eye movements
% 30 has 72% of trials removed due to horizontal eye movements
% 33 has 54% of trials removed due to horizontal eye movements
% 8 has 33% of trials removed due to horizontal eye movements
% 11 has 35% of trials removed due to horizontal eye movements
% 
% exclude 12 because data failed to fit model
% exclude 6 because extreme outlier in BEH data
% 
% exp.participants = {'004';'005';'007';'009';'010';'013';...
%     '014';'015';'016';'017';'018';'019';'021';'022';'023';'024';'025';'026';...
%     '027';'029';'031';'032';'034';'035';'036';'037';'038';'039';'040'};

% *byTargets_v3:
% 6 has 36% of trials removed due to horizontal eye movements
% 8 has 35% of trials removed due to horizontal eye movements
% 11 has >30% of trials removed due to horizontal eye movements
% 23 has 33% of trials removed due to horizontal eye movements
% 28 has >30% of trials removed due to horizontal eye movements
% 30 has >30% of trials removed due to horizontal eye movements
% 33 has >30% of trials removed due to horizontal eye movements

% exclude 12 because data failed to fit model

exp.participants = {'004';'005';'007';'009';'010';'013';'014';'015';...
    '017';'018';'019';'021';'022';'024';'025';'026';'027';'029';'031';...
    '032';'034';'035';'036';'037';'038';'039';'040'};


% *byCues_v3:
% 6 has 35% of trials removed due to horizontal eye movements
% 8 has 35.6% of trials removed due to horizontal eye movements
% 11 has 53% of trials removed due to horizontal eye movements
% 23 has 34% of trials removed due to horizontal eye movements
% 28 has 55.6% of trials removed due to horizontal eye movements
% 30 has 73% of trials removed due to horizontal eye movements
% 33 has 64% of trials removed due to horizontal eye movements

% exclude 12 because data failed to fit model

exp.participants = {'004';'005';'007';'009';'010';'013';'014';'015';...
    '017';'018';'019';'021';'022';'024';'025';'026';'027';'029';'031';...
    '032';'034';'035';'036';'037';'038';'039';'040'};

% 
% *byTargets_v4:
% 6 has 49% of trials removed due to horizontal eye movements
% 8 has 43% of trials removed due to horizontal eye movements
% 11 has >30% of trials removed due to horizontal eye movements
% 23 has 44% of trials removed due to horizontal eye movements
% % 24 has 30.6% of trials removed due to horizontal eye movements
% 26 has 45% of trials removed due to horizontal eye movements
% 28 has >30% of trials removed due to horizontal eye movements
% 30 has >30% of trials removed due to horizontal eye movements
% 32 has 36% of trials removed due to horizontal eye movements
% 33 has >30% of trials removed due to horizontal eye movements
% 35 has 34% of trials removed due to horizontal eye movements
% 39 has 48% of trials removed due to horizontal eye movements
% 40 has 33% of trials removed due to horizontal eye movements
% 
% exclude 12 because data failed to fit model
% 
% exp.participants = {'004';'005';'007';'009';'010';'012';...
%     '013';'014';'015';'017';'018';'019';'021';'022';'024';'025';...
%     '027';'029';'031';'034';'036';'037';'038'};



%% Blink Correction
% the Blink Correction wants dissimilar events (different erps) seperated by 
% commas and similar events (similar erps) seperated with spaces. See 'help gratton_emcp'
% exp.selection_cards = {'11 21','13 23'};
%%%indicates where you want to center your data (where time zero is)
%%%%must be list == length(exp.setname) 
% exp.selection_cards = {'30','20'}; %target-aligned 
exp.selection_cards = {'2 3 9'}; %cue-align

%% Artifact rejection. 
% Choose the threshold to reject trials. More lenient threshold followed by an (optional) stricter threshold 
exp.preocularthresh = [-1000 1000]; %First happens before the ocular correction.
% exp.postocularthresh = [ ]; %Second happens after. Leave blank [] to skip
exp.postocularthresh = [-500 500]; %Second happens after. Leave blank [] to skip

%% Events and event labels
%events are what happen within each trial (e.g., fixation, target, response, etc...) 
%%%for each condition (lag 1-4 in this case), numbers correspond to
%%%triggers that will be kept for each condition. All other triggers will
%%%be removed
% You can collect multiple triggers into one event with square brackets [].

% exp.events = {[30],[20]}; %can be list or matrix (sets x events)    
% exp.event_names = {'LT','RT'}; %must be list or matrix (sets x events)
% exp.suffix = {'byTarg'};

% exp.events = {[2 3 9]}; %can be list or matrix (sets x events)    
exp.event_names = {'Cues'}; %must be list or matrix (sets x events)
exp.suffix = {'byCue'};

%% Electrode location
%Where are your electrodes? (.ced file)
exp.electrode_locs = [pwd '\EOG-electrode-locs-32_orientwheel.ced'];
% exp.electrode_locs = 'M:\Analysis\OrientWheel_Exo\EOG-electrode-locs-32_orientwheel.ced';
% electrode information
exp.electrode = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32];
exp.elec_names = {'Oz';'Pz';'Cz';'FCz';'Fz';'O1';'O2';'PO3';'PO4';'P7';'P8';'P5';'P6';'P3';'P4';'CP5';...
    'CP6';'CP1';'CP2';'C3';'C4';'FC5';'FC6';'FC1';'FC2';'F7';'F8';'F3';'F4';'Fp1';'Fp2'};

%% Re-referencing the data
exp.refelec = 1; %which electrode do you want to re-reference to?
exp.brainelecs = [2:32]; %list of every electrode collecting brain data (exclude mastoid reference, EOGs, HR, EMG, etc.

%% Filter the data?
exp.filter = 'on'; %filter all files
exp.hicutoff = 50; %higher edge of the frequency pass band (Hz)
exp.locutoff = 0.1; %lower edge of the frequency pass band (Hz)

%% FFT/Wavelet settings
% How long is your window going to be? (Longer window == BETTER frequency 
% resolution & WORSE time resolution)
exp.winsize = 512; %use numbers that are 2^x, e.g. 2^10 == 1024ms

% Baseline will be subtracted from the power variable. It is relative to 
% your window size. Can use just NaN for no baseline
%e.g., [-200 0] will use [-200-exp.winsize/2 0-exp.winsize/2]; 
exp.erspbaseline = NaN;
% exp.erspbaseline = [-500 -200];

% Instead of choosing a windowsize, you can choose a number of cycles per 
% frequency for standard wavelet analysis: usually [3] for better temporal
% precision or [6] for better frequency precision.
% If [wavecycles factor], wavelet cycles increase with frequency beginning 
% at wavecyles. See "help popnewtimef"
exp.cycles = [0]; %leave it at 0 to use FFT
% exp.cycles = [2 0.7]; %number of cycles 

% Choose number of output times
exp.timesout = 300; %200 is usually used

% Set sampling factor for frequencies. 
% when exp.cycles==0, frequency spacing is (low_freq/padratio). For wavelet,
% multiplies the # of output freqs by dividing their spacing (2 is default).
% higher values give smooth looking plot but at computational cost (16 is
% very high)
exp.padratio = 4;

% What frequencies to consider?
exp.freqrange = [1 40]; 
% exp.freqrange = [exp.cycles(1) 40]; %when doing wavelet

%% Epoching the data
exp.epoch = 'off'; %on to epoch data; off to load previous data

%%%indicates where you want to center your data (where time zero is)
exp.epochs = {}; %must be list == length(exp.setname)
exp.epochs_name = {};

%in seconds; epoched trigger is 0 e.g. [-1 2]
exp.epochslims = [-1.2 2.2]; %cue-aligned
% exp.epochslims = [-2.4 1.4]; %target-aligned

%remove the baseline for each epoched set, in ms. e.g. [-200 0] 
%baseline is from after start of trial & before cue onset
exp.epochbaseline = [-200 0]; %cue-aligned
% exp.epochbaseline = [-1600 -1400]; %target-aligned


%% Time-Frequency settings
%Do you want to run time-frequency analyses? (on/off)
exp.tf = 'on';

%Do you want to use all the electrodes or just a few? Leave blank [] for 
% all (will use same as exp.brainelecs)
exp.tfelecs = [];

%Do you want to save the single-trial data? (on/off) (Memory intensive!!!)
exp.singletrials = 'on';
%Saving the single trial data is memory intensive. Just use the electrodes
% you need. 
exp.singletrialselecs = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32];
exp.singtrlelec_name = {'Oz';'Pz';'Cz';'FCz';'Fz';'O1';'O2';'PO3';'PO4';'P7';'P8';'P5';'P6';'P3';'P4';'CP5';...
    'CP6';'CP1';'CP2';'C3';'C4';'FC5';'FC6';'FC1';'FC2';'F7';'F8';'F3';'F4';'Fp1';'Fp2'};


%//////////////////////////////////////////////////////////////////////////
%% Save your pipeline settings
% The settings will be saved as a new folder. It lets you save multiple datasets with different preprocessing parameters.
exp.settings = char(exp.setname); %name settings
% `````````````````````````````````````````````````````````````````````````
% Saving will help you remember what settings were used in each dataset
save([exp.settings '_Settings'],'exp') %save these settings as a .mat file. 
%//////////////////////////////////////////////////////////////////////////


%% Run preprocessing code
tic %start timer
Preprocessing_OrientWheel_Exo(exp)
toc %end timer


