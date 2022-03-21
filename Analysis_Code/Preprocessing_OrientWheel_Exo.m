function Preprocessing_OrientWheel_Exo(exp)

% Not set-up to process data from the staircasing task

try
%     if parpool('poolsize') == 0
%         parpool OPEN 3;
%         parpool(3)
%     end
    
    nparts = length(exp.participants);
    nsets = length(exp.setname);
    
    % Replicating event names when exp.events is a matrix
%     if any(size(exp.event_names) ~= size(exp.events))
%         repfactor = int8(size(exp.events)./size(exp.event_names));
%         exp.event_names = repmat(exp.event_names, repfactor);
%     end

    % Replicating event triggers when exp.events is a matrix
    if isempty(exp.epochs) == 1
    %     exp.epochs = exp.events;
        exp.epochs = cellstr(num2str(cell2mat(reshape(exp.events,1,size(exp.events,1)*size(exp.events,2)) )'))';
    else
        exp.events = exp.epochs;
    end

    % Is epoch names not specified, use event names
    if isempty(exp.epochs_name) == 1
        exp.epochs_name = exp.event_names;
    else
        exp.event_names = exp.epochs_name;
    end
    
    
    for i_set = 1:nsets
        
        sprintf(exp.setname{i_set})
        
        % if folder doesn't exist yet, create one
        if ~exist([exp.pathname  '\' exp.setname{i_set} '\Segments\'])
            mkdir([exp.pathname  '\' exp.setname{i_set} '\Segments\']);
        end
        
        %initialize EEGLAB
        [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
        %subject numbers to analyze
        
        nevents = length(exp.events(i_set,:));
        
        %% Load data and channel locations
        for i_part = 1:nparts
%         for i_part = 18:nparts   
            
            part_name = exp.participants{i_part};
            
            if strcmp('on',exp.epoch) == 1
                
                sprintf(['Participant ' num2str(exp.participants{i_part})])
                
                %% Load a data file
                EEG = pop_loadbv(exp.pathname, [part_name '_Orient_Exo_Exp.vhdr']); %orientation wheel data only
                
                %% Load channel information
                EEG = pop_chanedit(EEG, 'load',{exp.electrode_locs 'filetype' 'autodetect'});
               
                %% Arithmetically re-reference to linked mastoid (M1 + M2)/2
                % only for brain electrodes (exclude mastoids & EOGs)
                for ii = exp.brainelecs(1):length(exp.brainelecs)
                    EEG.data(ii,:) = (EEG.data(ii,:)-((EEG.data(exp.refelec,:))*.5));
                end
                clear ii
                
                %% Filter the data
                if strcmpi(exp.filter,'on')
                   EEG = pop_eegfiltnew(EEG, exp.locutoff, exp.hicutoff); % filter function
                end
                
                %% Change markers so they can be used by the gratton_emcp script
                allevents = length(EEG.event);
                for i_event = 2:allevents %skip the first
                    EEG.event(i_event).type = num2str(str2num(EEG.event(i_event).type(2:end)));
                end

                %% The triggers are early
                [EEG] = VpixxEarlyTriggerFix(EEG);
                
                %% Extract epochs of data time locked to event
                %Extract data time locked to targets and remove all other events
                EEG = pop_epoch(EEG, exp.epochs, exp.epochslims, 'newname', [part_name '_epochs'], 'epochinfo', 'yes');
                %subtract baseline
                EEG = pop_rmbase(EEG, exp.epochbaseline);
  
                %% Remove practice trials from data
                rej_practice = zeros(1,length(EEG.epoch));
                rej_practice(1,1:20) = 1; %mark the first 20 trials for removal
                EEG = pop_rejepoch(EEG, rej_practice, 0);
                clear rej_practice
                
                %% Remove trials due to recording errors
                   % (see paper doc for details)
               if strcmpi(part_name,'009') 
%                     %also need BEH trial removal 
                    rej_oops = zeros(1,length(EEG.epoch));
                    rej_oops(1,344:end) = 1; %mark mismatched trials
%                     rej_oops(1,345:347) = 1; %mark mismatched trials
                    EEG = pop_rejepoch(EEG, rej_oops, 0);
                    clear rej_oops
               end
%                 if (strcmpi(exp.setname,'byTargets_v1') ||...
%                         strcmpi(exp.setname,'byTargets_v2') ||...
%                         strcmpi(exp.setname,'byTargets_v3')) &&...
%                         strcmpi(part_name,'009')
%                     %also need BEH trial removal 
%                     rej_oops = zeros(1,length(EEG.epoch));
%                     rej_oops(1,344:346) = 1; %mark mismatched trials
%                     EEG = pop_rejepoch(EEG, rej_oops, 0);
%                     clear rej_oops  
%                 elseif strcmpi(exp.setname,'byCues_v1') && strcmpi(part_name,'009') 
%                     %also need BEH trial removal 
%                     rej_oops = zeros(1,length(EEG.epoch));
%                     rej_oops(1,344:347) = 1; %mark mismatched trials
% %                     rej_oops(1,345:347) = 1; %mark mismatched trials
%                     EEG = pop_rejepoch(EEG, rej_oops, 0);
%                     clear rej_oops
%                 elseif strcmpi(exp.setname,'byCues_v2') && strcmpi(part_name,'009') 
%                     %also need BEH trial removal 
%                     rej_oops = zeros(1,length(EEG.epoch));
%                     rej_oops(1,344:347) = 1; %mark mismatched trials
% %                     rej_oops(1,345:347) = 1; %mark mismatched trials
%                     EEG = pop_rejepoch(EEG, rej_oops, 0);
%                     clear rej_oops    
%                 end
                
                %% Get behavior data and add to EEG structure
                [error_deg] = getBEHdata_OrientWheel_Exo(part_name,exp); %load BEH data
                
                % add error deg to epoch structure
                if (length(EEG.epoch) == length(error_deg)) ||... %make sure right BEH file
                   strcmpi(part_name,'009')
                    EEG.error_deg = error_deg;
                end
                clear error_deg

% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                %% Artifact Rejection, EMCP Correction, then 2nd Rejection
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  

                % :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                % Artifact rejection 1, trials with range >exp.preocularthresh uV
                if isempty(exp.preocularthresh) == 0
                    rejtrial = struct([]);
                    [EEG Indexes] = pop_eegthresh(EEG,1,[1:size(EEG.data,1)],exp.preocularthresh(1),exp.preocularthresh(2),EEG.xmin,EEG.xmax,0,1);
                    % ````````````````````````````````````````````````````` 
                    rejtrial(i_set,1).ids = find(EEG.reject.rejthresh==1);
                end
                % :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                
                % :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                % EMCP occular correction
                temp_ocular = EEG.data(end-1:end,:,:); %to save the EYE data for after
                EEG = gratton_emcp(EEG, exp.selection_cards, {'VEOG'},{'HEOG'}); %this assumes the eye channels are called this
                EEG.emcp.table %this prints out the regression coefficients
                EEG.data(end-1:end,:,:) = temp_ocular; %replace the eye data
                % :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                
                % :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                % Baseline again since EMCP changed it
                EEG = pop_rmbase(EEG,exp.epochbaseline);
                % :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                
                % :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                % Artifact rejection 2, trials with range >exp.postocularthresh uV
                if isempty(exp.postocularthresh) == 0
                    [EEG Indexes] = pop_eegthresh(EEG,1,[1:size(EEG.data,1)-2],exp.postocularthresh(1),exp.postocularthresh(2),EEG.xmin,EEG.xmax,0,1);
                    rejtrial(i_set,2).ids = find(EEG.reject.rejthresh==1);
                end
                % :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
                %% Keep track of order of rejection procedures
                i_rej = 3; %this will always be 3 even if the threshold rejects are not used
                
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
% ````````````````````````````````````````````````````````````````````````````````````````````  
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
                %% Reject trials with horizontal eye movements
                [EEG,rejtrial] = mark_Heyemove(EEG,rejtrial,exp,i_rej);
                
                i_rej = i_rej + 1; %keep track of order of rejection procedures
                
                % :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                %% Additional rejection of trials (reviewed visually) 
                [EEG,rejtrial] = reject_trials_inspect_Exo(EEG,part_name,exp,rejtrial,i_rej);
                
                i_rej = i_rej + 1; %keep track of order of rejection procedures
                
                % :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                %% Replace the stored data with this new set
                tempEEG = EEG;             
                
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
% ````````````````````````````````````````````````````````````````````````````````````````````  
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                 
                %% Select individual events
                for i_event = 1:nevents
                    EEG = pop_selectevent( tempEEG, 'type', exp.events{i_set,i_event}, 'deleteevents','on','deleteepochs','on','invertepochs','off');
                    EEG = pop_editset(EEG, 'setname', [part_name '_' exp.event_names{i_set,i_event} '_' exp.setname{i_set}] );
                    EEG = pop_saveset( EEG, 'filename',[part_name '_' exp.event_names{i_set,i_event} '_' exp.setname{i_set} '.set'],'filepath',[exp.pathname '\' exp.setname{i_set} '\Segments\']);
                end
% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::                 
            end %create epoch loop end
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
            %% Time-Frequency Data
            if strcmp('on',exp.tf) == 1 || strcmp('on',exp.singletrials) == 1
                
                if ~exist([exp.pathname '\' exp.setname{i_set} '\TimeFrequency\'])
                    mkdir([exp.pathname '\' exp.setname{i_set} '\TimeFrequency\']);
                end

                
                for i_event = 1:nevents
                    
                    if strcmp('on',exp.epoch) == 0 %loading previous epochs if not created this session
                        filename = [part_name '_' exp.event_names{i_set,i_event} '_' exp.setname{i_set}];
                        EEG = pop_loadset('filename',[filename '.set'],'filepath',[exp.pathname  '\' exp.setname{i_set} '\Segments\']);
                    end
                    
                    if isempty(exp.tfelecs) %if TF electrodes not specified, same as exp.brainelecs
                       exp.tfelecs = exp.brainelecs; 
                    end
                    
                    for i_tf = 1:length(exp.tfelecs)
                        i_chan = exp.tfelecs(i_tf);
                        EEG = eeg_checkset(EEG);
                        [ersp(i_chan,:,:),itc(i_chan,:,:),powbase,times,freqs,dum1,dum2,all_ersp(i_chan).trials] =...
                            pop_newtimef(EEG, 1, i_chan, exp.epochslims*1000, exp.cycles, ...
                            'topovec', i_chan, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo,...
                            'baseline', exp.erspbaseline, 'freqs', exp.freqrange, 'freqscale', 'linear', ...
                            'padratio', exp.padratio, 'plotphase','off','plotitc','off','plotersp','off',...
                            'winsize',exp.winsize,'timesout',exp.timesout);
                    end
                    clear i_chan i_tf

                    if strcmp('on',exp.tf) == 1 %if TF was done already, do not save
                        save([exp.pathname '\' exp.setname{i_set} '\TimeFrequency\' part_name '_' exp.event_names{i_set,i_event} '_' exp.setname{i_set} '.mat'],'ersp','itc','times','freqs','powbase','exp')
                    end
                        
                    
                     % Save single trial data
                    if strcmp('on',exp.singletrials) == 1
                        
                        % Create folder for single trial data
                        if ~exist([exp.pathname '\' exp.setname{i_set} '\SingleTrials\' part_name '\'],'dir')
                            mkdir([exp.pathname '\' exp.setname{i_set} '\SingleTrials\' part_name '\']);
                        end
                        
                        % File path name
                        Filepath_Trials = [exp.pathname '\' exp.setname{i_set} '\SingleTrials\' part_name '\'];
                        
                        % Save single trial data from the selected electrodes
                        for zzz = 1:length(exp.singletrialselecs)
                            i_chan = exp.singletrialselecs(zzz);
                            elec_all_ersp = all_ersp(i_chan).trials;
                            save([Filepath_Trials exp.singtrlelec_name{zzz} '_SingleTrials_' exp.event_names{i_set,i_event} '_' exp.setname{i_set} '.mat'],...
                                'elec_all_ersp','times','freqs','powbase','exp')
                        end
                        clear i_chan elec_all_ersp zzz
                    end
                    clear Filepath_Trials

                end
                eeglab redraw
                
            end
            clear part_name i_rej
% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::             
        end %i_part loop end
% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::         
    end %i_set loop end
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 

    
% !!!!!!!!!!!!!!!!    
catch ME % !!!!!!!
    save('dump') 
    throw(ME)
end % !!!!!!!!!!!!
% !!!!!!!!!!!!!!!! 


% ####################
end % end function ###
% ####################
