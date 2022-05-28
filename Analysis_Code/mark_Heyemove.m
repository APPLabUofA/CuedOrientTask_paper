function [EEG,rejtrial] = mark_Heyemove(EEG,rejtrial,exp,i_rej)
% Procedure is adapted from methods describe in this paper: 
% http://journals.sagepub.com/doi/10.1177/0956797619830384
% Code from paper was downloaded from here: https://osf.io/ws3j9/
% 
% Created/adapted by SSS (12/2020)
% #########################################################################

fprintf('checking for eye movements... \n')

i_set = 1; %code changed so now always 1


%% EEG information
srate = EEG.srate; %sampling rate
times = EEG.times; %time in ms

%% Time window
% Only care about eye movements during relevant time period
%this needs to be updated if epochs &/or event alignment are changed
if strcmpi(exp.setname,'byTargets_v2')  
    timewin = [0 600]; %target-aligned: target-onset 2 response window onset
elseif strcmpi(exp.setname,'byTargets_v3')   
    timewin = [-1600 600]; %target-aligned: baseline 2 response window onset
elseif strcmpi(exp.setname,'byTargets_v4')   
    timewin = [-1600 600]; %target-aligned: baseline 2 response window onset
elseif strcmpi(exp.setname,'byCues_v2') 
    timewin = [1240 1950]; %cue-aligned: earliest target onset 2 response window onset of latest cue
elseif strcmpi(exp.setname,'byCues_v3') 
    timewin = [-200 1950]; %cue-aligned: baseline 2 response window onset of latest cue
else
    timewin = [min(times) max(times)]; %whole epoch
end

timelim = (times>=timewin(1) & times<=timewin(2));

%% Settings for moving window and rejection
winStep = 10; %size of steps time window moves across trial, ms (e.g., 50)
winSize = 100; %size of sliding time window, ms (e.g., 100)
thresh = 25; %threshold for what is considered an eye movement (e.g., 40)

%% Extract HEOG data
% Get electrode labels
[electrodes{1:size(EEG.chanlocs,2)}] = deal(EEG.chanlocs.labels);
electrodes = electrodes';

% Extract HEOG data
heog_data = squeeze(EEG.data(ismember(electrodes,'HEOG'),timelim,:));


%% Function that checks for horizontal eye movements 
markH = zeros(1,EEG.trials); %set-up for keeping track of rejected trials 
parfor i_trial = 1:EEG.trials
    rawTS = heog_data(:,i_trial); %extract single trial data
    eMoveH = step_t(rawTS,srate,winStep,winSize,thresh);   %actual function     
    markH(i_trial) = markDetect(eMoveH);  %mark trial as 0 or 1
end
clear i_trial


%% Actual rejection of trials using eeglab function
rejtrial(i_set,i_rej).ids = find(markH==1);
EEG = pop_rejepoch(EEG,markH,0);
        
% save rejected trials
EEG.rejtrial = rejtrial;


%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end %end mark_Heyemove function



% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
% `````````````````````````````````````````````````````````````````````````  
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 











