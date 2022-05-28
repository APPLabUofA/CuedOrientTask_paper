% /////////////////////////////////////////////////////////////////////////
%% Load processing settings
load('byTargets_v3_Settings.mat');
% load('byCues_v3_Settings.mat');

% /////////////////////////////////////////////////////////////////////////
% Load specific EEG dataset to make EEGLab happy
% EEG = pop_loadset('040_RT_byTargets_v3.set');
% eeglab redraw

% /////////////////////////////////////////////////////////////////////////
%% Load behavioral data
load([pwd '/BEH_' exp.settings '.mat']);

% /////////////////////////////////////////////////////////////////////////
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% Load & convert processed fooof data to .mat file
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% /////////////////////////////////////////////////////////////////////////

% Folder to save the power spectra out to mat files
% fooof_dir = [pwd '\fooof_trl\fix\' exp.settings];
% fooof_dir = [pwd '\fooof_trl\mid\' exp.settings];
% fooof_dir = [pwd '\fooof_trl\late\' exp.settings];
fooof_dir = [pwd '\fooof_trl\post\' exp.settings];

% Remember current folder
prime_dir = pwd;

% Folder the saved fooof are located
results_dir = [fooof_dir '\fooof_results\'];


% for i_part = length(exp.participants)
for i_part = 1:length(exp.participants)    
    
    % Extract by electrode
    for ii = 1:length(exp.electrode)
%         i_elect = exp.electrode(ii); %for doing only a selection of electrodes
        
        % Create new folder for each participant
        tmp_loc = [results_dir '\' exp.participants{i_part} '\' exp.elec_names{ii} '\'];
        
        %.....................................................................
        
        % Create new folder for condition
        tmpfile_loc = [tmp_loc '\L_v\'];
        
        cd(tmpfile_loc) % Go to folder with saved results
        list = cellstr(ls('*.json')); % Get list of result file names

        % Loop through file list
        for ifile = 1:length(list)
            fname = list{ifile};
            dat = importdata(fname);

            fooof_results = []; %pre-allocate
            for ind = 1:length(dat)
                fooof_results = [fooof_results, jsondecode(dat{ind})];
            end
            clear dat ind

            save(replace(fname,'.json','.mat'),'fooof_results') % save results
            clear fname fooof_results
        end
        clear ifile tmpfile_loc list
        
        %.....................................................................
        
        % Create new folder for condition
        tmpfile_loc = [tmp_loc '\L_iv\'];
        
        cd(tmpfile_loc) % Go to folder with saved results
        list = cellstr(ls('*.json')); % Get list of result file names

        % Loop through file list
        for ifile = 1:length(list)
            fname = list{ifile};
            dat = importdata(fname);

            fooof_results = []; %pre-allocate
            for ind = 1:length(dat)
                fooof_results = [fooof_results, jsondecode(dat{ind})];
            end
            clear dat ind

            save(replace(fname,'.json','.mat'),'fooof_results') % save results
            clear fname fooof_results
        end
        clear ifile tmpfile_loc list
        
        %.....................................................................
        
        % Create new folder for condition
        tmpfile_loc = [tmp_loc '\R_v\'];
        
        cd(tmpfile_loc) % Go to folder with saved results
        list = cellstr(ls('*.json')); % Get list of result file names

        % Loop through file list
        for ifile = 1:length(list)
            fname = list{ifile};
            dat = importdata(fname);

            fooof_results = []; %pre-allocate
            for ind = 1:length(dat)
                fooof_results = [fooof_results, jsondecode(dat{ind})];
            end
            clear dat ind

            save(replace(fname,'.json','.mat'),'fooof_results') % save results
            clear fname fooof_results
        end
        clear ifile tmpfile_loc list
        
        %.....................................................................
        
        % Create new folder for condition
        tmpfile_loc = [tmp_loc '\R_iv\'];
        
        cd(tmpfile_loc) % Go to folder with saved results
        list = cellstr(ls('*.json')); % Get list of result file names

        % Loop through file list
        for ifile = 1:length(list)
            fname = list{ifile};
            dat = importdata(fname);

            fooof_results = []; %pre-allocate
            for ind = 1:length(dat)
                fooof_results = [fooof_results, jsondecode(dat{ind})];
            end
            clear dat ind

            save(replace(fname,'.json','.mat'),'fooof_results') % save results
            clear fname fooof_results
        end
        clear ifile tmpfile_loc list
        
        %.....................................................................
        
    end
    clear ii i_elect tmp_loc
end
clear i_part

cd(prime_dir) %return to main directory

clear fooof_dir results_dir



% -------------------------------------------------------------------------
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% Load fooof band .mat files
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% -------------------------------------------------------------------------
% Each band peak contains threee values: 
%              [CF (center frequency), PW (power), BW (bandwidth)]

% Folder to save the power spectra out to mat files
% fooof_dir = [pwd '\fooof_trl\fix\' exp.settings];
% fooof_dir = [pwd '\fooof_trl\mid\' exp.settings];
% fooof_dir = [pwd '\fooof_trl\late\' exp.settings];
fooof_dir = [pwd '\fooof_trl\post\' exp.settings];


% Remember current folder
prime_dir = pwd;

% Folder the saved fooof are located
results_dir = [fooof_dir '\fooof_results\'];

% Condition names
cond_name = {'R_v'; 'R_iv'; 'L_v'; 'L_iv'};

% Saving data to cells
beta1_CF = cell(length(cond_name),length(exp.participants)); % pre-allocate
alpha_CF = cell(length(cond_name),length(exp.participants)); % pre-allocate
theta_CF = cell(length(cond_name),length(exp.participants)); % pre-allocate
delta_CF = cell(length(cond_name),length(exp.participants)); % pre-allocate
beta1_PW = cell(length(cond_name),length(exp.participants)); % pre-allocate
alpha_PW = cell(length(cond_name),length(exp.participants)); % pre-allocate
theta_PW = cell(length(cond_name),length(exp.participants)); % pre-allocate
delta_PW = cell(length(cond_name),length(exp.participants)); % pre-allocate
for i_part = 1:length(exp.participants)
    
    for ielect = 1:length(exp.elec_names)
        
        % Results folder for each electrode
        tmp_loc = [results_dir exp.participants{i_part} '\' exp.elec_names{ielect} '\'];
        
        for icond = 1:length(cond_name)
        
            % Folder of files in condition
            tmpfile_loc = [tmp_loc cond_name{icond} '\'];

            list = cellstr(ls([tmpfile_loc 'bands_*'])); % Get list of result file names
            list = list(~contains(list,'json')); % Only the .mat files

            % Loop through file list
            for ifile = 1:length(list)
                load([tmpfile_loc list{ifile}]) %load data

                % Create CF (center frequency) variables w/cond names 
                eval(strcat('beta1_CF_', cond_name{icond}, '{ielect,i_part}(:,ifile)=squeeze(beta1(:,1));'))
                eval(strcat('alpha_CF_', cond_name{icond}, '{ielect,i_part}(:,ifile)=squeeze(alphas(:,1));'))
                eval(strcat('theta_CF_', cond_name{icond}, '{ielect,i_part}(:,ifile)=squeeze(thetas(:,1));'))
                eval(strcat('delta_CF_', cond_name{icond}, '{ielect,i_part}(:,ifile)=squeeze(delta(:,1));'))

                beta1_CF{icond,i_part}(ielect,ifile,:) = squeeze(beta1(:,1));
                alpha_CF{icond,i_part}(ielect,ifile,:) = squeeze(alphas(:,1));
                theta_CF{icond,i_part}(ielect,ifile,:) = squeeze(thetas(:,1));
                delta_CF{icond,i_part}(ielect,ifile,:) = squeeze(delta(:,1));

                % Create PW (power) variables w/cond names 
                eval(strcat('beta1_PW_', cond_name{icond}, '{ielect,i_part}(:,ifile)=squeeze(beta1(:,2));'))
                eval(strcat('alpha_PW_', cond_name{icond}, '{ielect,i_part}(:,ifile)=squeeze(alphas(:,2));'))
                eval(strcat('theta_PW_', cond_name{icond}, '{ielect,i_part}(:,ifile)=squeeze(thetas(:,2));'))
                eval(strcat('delta_PW_', cond_name{icond}, '{ielect,i_part}(:,ifile)=squeeze(delta(:,2));'))

                beta1_PW{icond,i_part}(ielect,ifile,:) = squeeze(beta1(:,2));
                alpha_PW{icond,i_part}(ielect,ifile,:) = squeeze(alphas(:,2));
                theta_PW{icond,i_part}(ielect,ifile,:) = squeeze(thetas(:,2));
                delta_PW{icond,i_part}(ielect,ifile,:) = squeeze(delta(:,2));

                clear alphas thetas beta1 delta
                
            end
            clear ifile list tmpfile_loc
        end
        clear icond tmp_loc
    end
    clear ielect
end
clear i_part cond_name



% Save extracted data - each save below corresponds to a time window
% save(['bands_trl_fix_' exp.settings],'alpha_*','beta1_*','delta_*','theta_*')
% save(['bands_trl_mid_' exp.settings],'alpha_*','beta1_*','delta_*','theta_*')
% save(['bands_trl_late_' exp.settings],'alpha_*','beta1_*','delta_*','theta_*')
save(['bands_trl_post_' exp.settings],'alpha_*','beta1_*','delta_*','theta_*')


clear alpha_* beta1_* delta_* theta_* fooof_dir





% -------------------------------------------------------------------------
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%% Load fooof results .mat files
% '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% -------------------------------------------------------------------------
% Aperiodic parameters: [Offset, Exponent]
% Peak parameters:
%         CF (center frequency) - same as mean parameter of the Gaussian
%         PW (power) - height of the model fit above the aperiodic component
%         BW (bandwidth) - 2 * the standard deviation of the Gaussian


% Folder to save the power spectra out to mat files
% fooof_dir = [pwd '\fooof_trl\fix\' exp.settings];
% fooof_dir = [pwd '\fooof_trl\mid\' exp.settings];
% fooof_dir = [pwd '\fooof_trl\late\' exp.settings];
fooof_dir = [pwd '\fooof_trl\post\' exp.settings];

% Folder the saved fooof are located
results_dir = [fooof_dir '\fooof_results\'];

% Condition names
cond_name = {'R_v'; 'R_iv'; 'L_v'; 'L_iv'};

% Saving data to cells
aperiodic = cell(length(cond_name),length(exp.participants)); %pre-allocate
r_square  = cell(length(cond_name),length(exp.participants)); %pre-allocate
fit_error = cell(length(cond_name),length(exp.participants)); %pre-allocate
peaks  = cell(length(cond_name),length(exp.participants),length(exp.elec_names)); %pre-allocate
for i_part = 1:length(exp.participants)
    
    for ielect = 1:length(exp.elec_names)
        
        % Results folder for each electrode
        tmp_loc = [results_dir exp.participants{i_part} '\' exp.elec_names{ielect} '\'];
        
        for icond = 1:length(cond_name)
        
            % Folder of files in condition
            tmpfile_loc = [tmp_loc cond_name{icond} '\'];

            list = cellstr(ls([tmpfile_loc 'f_results_*'])); % Get list of result file names
            list = list(~contains(list,'json')); % Only the .mat files

            % Loop through file list
            peak_tmp = NaN(length(list),8,2); %pre-allocate
            for ifile = 1:length(list)
                load([tmpfile_loc list{ifile}]) %load data
                
                aperiodic{icond,i_part}(ielect,ifile,:) = fooof_results.aperiodic_params_;
                r_square{icond,i_part}(ielect,ifile,:) = fooof_results.r_squared_;
                fit_error{icond,i_part}(ielect,ifile,:) = fooof_results.error_;
                
                % Adjusted for variable outputs (can be 1-8 peaks)
                if isempty(fooof_results.peak_params_)
                    %skip if no peaks found
                else
                    peak_tmp(ifile,1:size(fooof_results.peak_params_,1),:) = fooof_results.peak_params_(:,1:2);
                end


                eval(strcat('aperiodic_', cond_name{icond}, '{ielect,i_part}(1:2,ifile)=fooof_results.aperiodic_params_;'))
                eval(strcat('r_square_', cond_name{icond}, '{ielect,i_part}(:,ifile)=fooof_results.r_squared_;'))
                eval(strcat('fit_error_', cond_name{icond}, '{ielect,i_part}(:,ifile)=fooof_results.error_;'))

                clear fooof_results
                
            end
            peaks{icond,i_part,ielect} = peak_tmp;
            
            clear ifile list tmpfile_loc peak_tmp
        end
        clear icond tmp_loc
    end
    clear ielect
end
clear i_part cond_name



% Save extracted data - each save below corresponds to a time window
% save(['fits_trl_fix_' exp.settings],'peaks','aperiodic*','r_square*','fit_error*')
% save(['fits_trl_mid_' exp.settings],'peaks','aperiodic*','r_square*','fit_error*')
% save(['fits_trl_late_' exp.settings],'peaks','aperiodic*','r_square*','fit_error*')
save(['fits_trl_post_' exp.settings],'peaks','aperiodic*','r_square*','fit_error*')


clear peaks aperiodic* r_square* fit_error* fooof_dir results_dir




