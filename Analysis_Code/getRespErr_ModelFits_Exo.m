function [resp_errdeg, model_out] = getRespErr_ModelFits_Exo(exp,ALLEEG)
% Must load data using LoadProcData_OrientWheel.m for function to work.
% ALLEEG structure should contain all the subjects' data.
% Returns degree response error on each trial with trials rejected during
% the EEG pre-processing pipeline removed. 
% Returns parameters from fitting the response errors to the mixed model.



% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------
%% Remove rejected trials
% Get BEH data for trials excluding trials that were rejected in the EEG
% preprocessing of the epochs
resp_errdeg = cell(length(exp.participants),1); %pre-allocate
for i_part = 1:length(exp.participants)
    
%     load(['M:\Data\OrientWheel_Exo\BEH\' num2str(exp.participants{i_part}) '_Orient_Exo.mat'])
    load(['C:\Users\ssshe\Documents\MathLab\Data\OrientWheel_Exo\BEH\' num2str(exp.participants{i_part}) '_Orient_Exo.mat'])

    
    [n,m] = size(ALLEEG(i_part).rejtrial);
    % Get list of rejected trials
    pip = 1;
    for ni = 1:n %for when there are more than 1 column
        for mi = 1:m
            if ~isempty(ALLEEG(i_part).rejtrial(ni,mi).ids)
                rejlist{pip} = ALLEEG(i_part).rejtrial(ni,mi).ids;
                pip = 1 + pip;
            end
        end
        clear mi
    end
    if pip > 1 %if trials were rejected
        err_deg_tmp = ALLEEG(i_part).error_deg; %start with all the errors
        % each set of rejected trials needs to be removed in order
        % sequentially
        for mi = 1:length(rejlist)
            tmplist = [rejlist{mi}];
            err_deg_tmp(tmplist) = []; %removes the trials
            clear tmplist
        end
        clear mi
    elseif pip == 1 %if no trials were rejected, rejlist variable not created
        err_deg_tmp = ALLEEG(i_part).error_deg;
    end
    % create variable with selected BEH 
    resp_errdeg{i_part} = err_deg_tmp;
    
    clear rejlist n m err_deg_tmp pip ni
end
clear i_part




% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------
%% Fit all errors to mixed model
% model1 = StandardMixtureModel(); %standard 2 parameter model
model2 = WithBias(StandardMixtureModel); % model with mu
model_out = cell(1,length(exp.participants)); %pre-allocate
for ii = 1:length(exp.participants)
%         model_out{ii} = MemFit(errordeg{ii},model2); %plots fits  
        model_out{ii} = MLE(errordeg{ii},model2); %fits without plotting - w/bias
%         model_out{ii} = model_out{ii}(1:3);
end
clear ii error_deg model model1 model2


% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------
%% Fit errors by target side to mixed model
% model1 = StandardMixtureModel(); %standard 2 parameter model
model2 = WithBias(StandardMixtureModel); % model with mu
model_out_L = cell(1,length(exp.participants)); %pre-allocate
model_out_R = cell(1,length(exp.participants)); %pre-allocate
for ii = 1:length(exp.participants)
%         model_out_L{ii} = MemFit(errordeg{ii}(position{ii}==0),model2); %plots fits  
        model_out_L{ii} = MLE(errordeg{ii}(position{ii}==0),model2); %fits without plotting - w/bias
        
%         model_out_R{ii} = MemFit(errordeg{ii}(position{ii}==1),model2); %plots fits  
        model_out_R{ii} = MLE(errordeg{ii}(position{ii}==1),model2); %fits without plotting - w/bias

end
clear ii error_deg model model1 model2


% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------
%% Fit errors by cue to mixed model
% model1 = StandardMixtureModel(); %standard 2 parameter model
model2 = WithBias(StandardMixtureModel); % model with mu
model_out_invalid = cell(1,length(exp.participants)); %pre-allocate
model_out_valid = cell(1,length(exp.participants)); %pre-allocate
for ii = 1:length(exp.participants)
%         model_out_invalid{ii} = MemFit(errordeg{ii}(valid{ii}==0),model2); %plots fits  
        model_out_invalid{ii} = MLE(errordeg{ii}(valid{ii}==0),model2); %fits without plotting - w/bias
        
%         model_out_valid{ii} = MemFit(errordeg{ii}(valid{ii}==1),model2); %plots fits  
        model_out_valid{ii} = MLE(errordeg{ii}(valid{ii}==1),model2); %fits without plotting - w/bias

end
clear ii error_deg model model1 model2



% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------
%% Fit errors by target side & cue to mixed model
% model1 = StandardMixtureModel(); %standard 2 parameter model
model2 = WithBias(StandardMixtureModel); % model with mu
model_out_invalidL = cell(1,length(exp.participants)); %pre-allocate
model_out_validL = cell(1,length(exp.participants)); %pre-allocate
model_out_invalidR = cell(1,length(exp.participants)); %pre-allocate
model_out_validR = cell(1,length(exp.participants)); %pre-allocate
for ii = 1:length(exp.participants)
%         model_out_invalidL{ii} = MemFit(errordeg{ii}(valid{ii}==0 & position{ii}==0),model2); %plots fits  
        model_out_invalidL{ii} = MLE(errordeg{ii}(valid{ii}==0 & position{ii}==0),model2); %fits without plotting - w/bias
%         model_out_invalidR{ii} = MemFit(errordeg{ii}(valid{ii}==0 & position{ii}==1),model2); %plots fits  
        model_out_invalidR{ii} = MLE(errordeg{ii}(valid{ii}==0 & position{ii}==1),model2); %fits without plotting - w/bias
        
%         model_out_validL{ii} = MemFit(errordeg{ii}(valid{ii}==1 & position{ii}==0),model2); %plots fits  
        model_out_validL{ii} = MLE(errordeg{ii}(valid{ii}==1 & position{ii}==0),model2); %fits without plotting - w/bias
%         model_out_validR{ii} = MemFit(errordeg{ii}(valid{ii}==1 & position{ii}==1),model2); %plots fits  
        model_out_validR{ii} = MLE(errordeg{ii}(valid{ii}==1 & position{ii}==1),model2); %fits without plotting - w/bias

end
clear ii error_deg model model1 model2



% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////











