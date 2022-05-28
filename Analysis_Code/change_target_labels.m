function [EEG] = change_target_labels(EEG)


allevents = length(EEG.event);
for i_event = 1:allevents 
    
    if strcmpi(EEG.event(i_event).type,'30') %trials where target left
        
       ii = 0;
       
       %find closest cue trigger prior to target trigger
        while ~strcmpi(EEG.event(i_event+ii).type,'61') &&...
              ~strcmpi(EEG.event(i_event+ii).type,'70') &&...
              ~strcmpi(EEG.event(i_event+ii).type,'81')
          
              ii = ii + 1;
              
              if strcmpi(EEG.event(i_event+ii).type,'61') ||...
                 strcmpi(EEG.event(i_event+ii).type,'70') ||...
                 strcmpi(EEG.event(i_event+ii).type,'81')
                 %combine labels to separate correct & incorrect turn trials
                 EEG.event(i_event).type = strcat('10',EEG.event(i_event+ii).type);
              end
        end
    end
end


















