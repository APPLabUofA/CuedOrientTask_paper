
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%% Load previously processed target-aligned epoch data

% load settings
load('byTargets_v3_Settings.mat');

% load specific EEG dataset to make EEGLab happy
% EEG = pop_loadset('040_RT_byTargets_v3.set');
% eeglab redraw

% load behavioral data
load([pwd '/BEH_' exp.settings '.mat']);
load([pwd '/BEHerr_M_' exp.settings '.mat']);

% exp is a Matlab function that might get used below
exp2 = exp;
clear exp

% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

% Load saved phase data created in Anal_PhaseOpp_byTarget.m
load(['phase_cond_' exp2.settings '.mat'])

% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

% Load previously analyzed data (usually done just to make figures)
% load('witpcz_out_1e4perm.mat')

% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::




% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
%% Calculate wITPC_z maps per subject
% wITPC is response error modulating the length of phase angles

% specify values  
n_permutes = 1e3; %number of permutations used to estimate null distribution (takes really long time with increasing number)
% voxel_pval = 0.05; %uncorrect pixel-level threshold
mcc_voxel_pval = 0.05/4; % mcc = multiple comparisons correction

% ------------------------------------------------------------------------- 
% Test only relevant times
% timewin = [-1000 0]; 
timewin = [0 600]; %target onset to after response screen onset
timelim = find(times>=timewin(1) & times<=timewin(2));
% ------------------------------------------------------------------------- 

% pre-allocate variables
witpc_Liv = NaN(length(exp2.participants),length(exp2.singletrialselecs),length(freqs),length(timelim)); %pre-allocate
witpc_Riv = NaN(length(exp2.participants),length(exp2.singletrialselecs),length(freqs),length(timelim)); %pre-allocate
witpc_Lv  = NaN(length(exp2.participants),length(exp2.singletrialselecs),length(freqs),length(timelim)); %pre-allocate
witpc_Rv  = NaN(length(exp2.participants),length(exp2.singletrialselecs),length(freqs),length(timelim)); %pre-allocate

zmap_out_Liv = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
zmap_out_Riv = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
zmap_out_Lv  = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
zmap_out_Rv  = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate

zmap_witpc_thresh_Liv = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
zmap_witpc_thresh_Riv = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
zmap_witpc_thresh_Lv  = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
zmap_witpc_thresh_Rv  = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate

witpc_threshold_out_Liv = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
witpc_threshold_out_Riv = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
witpc_threshold_out_Lv  = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
witpc_threshold_out_Rv  = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate

threshold_out_Liv = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
threshold_out_Riv = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
threshold_out_Lv  = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
threshold_out_Rv  = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate

p_witpc_z_Liv = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
p_witpc_z_Riv = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
p_witpc_z_Lv  = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
p_witpc_z_Rv  = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate

p_witpcz_thresh_Liv = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
p_witpcz_thresh_Riv = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
p_witpcz_thresh_Lv  = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
p_witpcz_thresh_Rv  = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate

p_witpc_thresh_Liv = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
p_witpc_thresh_Riv = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
p_witpc_thresh_Lv  = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate
p_witpc_thresh_Rv  = cell(length(exp2.participants),length(exp2.singletrialselecs)); %pre-allocate

% Perform analysis on each participant
for i_part = 1:length(exp2.participants)
   
    % abs response errors
%     tmp_err_Liv = abs(err_invalidL{i_part}); %abs value of response errors
%     tmp_err_Riv = abs(err_invalidR{i_part}); %abs value of response errors
%     tmp_err_Lv  = abs(err_validL{i_part}); %abs value of response errors
%     tmp_err_Rv  = abs(err_validR{i_part}); %abs value of response errors

    % response errors
    tmp_err_Liv = err_invalidL{i_part}; 
    tmp_err_Riv = err_invalidR{i_part}; 
    tmp_err_Lv  = err_validL{i_part}; 
    tmp_err_Rv  = err_validR{i_part}; 

    % probability of guess trial
%     tmp_err_Liv = M_invalidL{i_part}(:,2)'; 
%     tmp_err_Riv = M_invalidR{i_part}(:,2)'; 
%     tmp_err_Lv  = M_validL{i_part}(:,2)'; 
%     tmp_err_Rv  = M_validR{i_part}(:,2)'; 
    
    n_resp_Liv = length(tmp_err_Liv); %number of trials
    n_resp_Riv = length(tmp_err_Riv); %number of trials
    n_resp_Lv  = length(tmp_err_Lv); %number of trials
    n_resp_Rv  = length(tmp_err_Rv); %number of trials
    
    for ii = 1:length(exp2.singletrialselecs)
        i_elect = exp2.singletrialselecs(ii); %for doing only a selection of electrodes
        
        %create a freq x time matrix of the response error values
        tmp_resperr_Liv = permute(repmat(tmp_err_Liv,[length(timelim) 1 length(freqs)]),[3 1 2]);
        tmp_resperr_Riv = permute(repmat(tmp_err_Riv,[length(timelim) 1 length(freqs)]),[3 1 2]);
        tmp_resperr_Lv  = permute(repmat(tmp_err_Lv,[length(timelim) 1 length(freqs)]),[3 1 2]);
        tmp_resperr_Rv  = permute(repmat(tmp_err_Rv,[length(timelim) 1 length(freqs)]),[3 1 2]);
        
        %compute wITPC
        witpc_Liv(i_part,ii,:,:) = squeeze(abs(nanmean(tmp_resperr_Liv.*squeeze(exp(1i*phase_Liv{i_part,i_elect}(:,timelim,:))),3)));
        witpc_Riv(i_part,ii,:,:) = squeeze(abs(nanmean(tmp_resperr_Riv.*squeeze(exp(1i*phase_Riv{i_part,i_elect}(:,timelim,:))),3)));
        witpc_Lv(i_part,ii,:,:)  = squeeze(abs(nanmean(tmp_resperr_Lv.*squeeze(exp(1i*phase_Lv{i_part,i_elect}(:,timelim,:))),3)));
        witpc_Rv(i_part,ii,:,:)  = squeeze(abs(nanmean(tmp_resperr_Rv.*squeeze(exp(1i*phase_Rv{i_part,i_elect}(:,timelim,:))),3)));
        
        clear tmp_resperr_*
        
        
        %% Create H0 distribution
        perm_witpc_Liv  = NaN(n_permutes,length(freqs),length(timelim)); % initialize null hypothesis matrices
        perm_witpc_Riv  = NaN(n_permutes,length(freqs),length(timelim)); % initialize null hypothesis matrices
        perm_witpc_Lv   = NaN(n_permutes,length(freqs),length(timelim)); % initialize null hypothesis matrices
        perm_witpc_Rv   = NaN(n_permutes,length(freqs),length(timelim)); % initialize null hypothesis matrices
        max_pixel_vals_Liv = NaN(n_permutes,1); % initialize null hypothesis matrices
        max_pixel_vals_Riv = NaN(n_permutes,1); % initialize null hypothesis matrices
        max_pixel_vals_Lv  = NaN(n_permutes,1); % initialize null hypothesis matrices
        max_pixel_vals_Rv  = NaN(n_permutes,1); % initialize null hypothesis matrices
        % generate pixel-specific null hypothesis parameter distributions
        for permi = 1:n_permutes
            
            fake_err_map_Liv = tmp_err_Liv(randperm(n_resp_Liv)); %randomize order of response errors
            fake_err_map_Riv = tmp_err_Riv(randperm(n_resp_Riv)); %randomize order of response errors
            fake_err_map_Lv = tmp_err_Lv(randperm(n_resp_Lv)); %randomize order of response errors
            fake_err_map_Rv = tmp_err_Rv(randperm(n_resp_Rv)); %randomize order of response errors

            %create a freq x time matrix of the response error values
            tmp_resperr_Liv = permute(repmat(fake_err_map_Liv,[length(timelim) 1 length(freqs)]),[3 1 2]);
            tmp_resperr_Riv = permute(repmat(fake_err_map_Riv,[length(timelim) 1 length(freqs)]),[3 1 2]);
            tmp_resperr_Lv = permute(repmat(fake_err_map_Lv,[length(timelim) 1 length(freqs)]),[3 1 2]);
            tmp_resperr_Rv = permute(repmat(fake_err_map_Rv,[length(timelim) 1 length(freqs)]),[3 1 2]);

            % save all permuted values
            perm_witpc_Liv(permi,:,:) = squeeze(abs(nanmean(tmp_resperr_Liv.*squeeze(exp(1i*phase_Liv{i_part,i_elect}(:,timelim,:))),3)));
            perm_witpc_Riv(permi,:,:) = squeeze(abs(nanmean(tmp_resperr_Riv.*squeeze(exp(1i*phase_Riv{i_part,i_elect}(:,timelim,:))),3)));
            perm_witpc_Lv(permi,:,:) = squeeze(abs(nanmean(tmp_resperr_Lv.*squeeze(exp(1i*phase_Lv{i_part,i_elect}(:,timelim,:))),3)));
            perm_witpc_Rv(permi,:,:) = squeeze(abs(nanmean(tmp_resperr_Rv.*squeeze(exp(1i*phase_Rv{i_part,i_elect}(:,timelim,:))),3)));

            max_pixel_vals_Liv(permi) = max(max(perm_witpc_Liv(permi,:,:)));
            max_pixel_vals_Riv(permi) = max(max(perm_witpc_Riv(permi,:,:)));
            max_pixel_vals_Lv(permi) = max(max(perm_witpc_Lv(permi,:,:)));
            max_pixel_vals_Rv(permi) = max(max(perm_witpc_Rv(permi,:,:)));
            
            clear fake_err_map_* tmp_resperr_*
        end
        clear permi
        
        %% Real z-map
        % now compute Z-map (standardized units from each value away from
        % the distribution of null-hypothesis values)
        permmean_Liv = squeeze(mean(perm_witpc_Liv,1));
        permmean_Riv = squeeze(mean(perm_witpc_Riv,1));
        permmean_Lv = squeeze(mean(perm_witpc_Lv,1));
        permmean_Rv = squeeze(mean(perm_witpc_Rv,1));
        
        permstd_Liv  = squeeze(std(perm_witpc_Liv,[],1));
        permstd_Riv  = squeeze(std(perm_witpc_Riv,[],1));
        permstd_Lv  = squeeze(std(perm_witpc_Lv,[],1));
        permstd_Rv  = squeeze(std(perm_witpc_Rv,[],1));
        
        witpc_z_Liv = (squeeze(squeeze(witpc_Liv(i_part,ii,:,:)))-permmean_Liv)./permstd_Liv;
        witpc_z_Riv = (squeeze(squeeze(witpc_Riv(i_part,ii,:,:)))-permmean_Riv)./permstd_Riv;
        witpc_z_Lv = (squeeze(squeeze(witpc_Lv(i_part,ii,:,:)))-permmean_Lv)./permstd_Lv;
        witpc_z_Rv = (squeeze(squeeze(witpc_Rv(i_part,ii,:,:)))-permmean_Rv)./permstd_Rv;
        
        zmap_out_Liv{i_part,ii} = witpc_z_Liv; %save values of plot
        zmap_out_Riv{i_part,ii} = witpc_z_Riv; %save values of plot
        zmap_out_Lv{i_part,ii} = witpc_z_Lv; %save values of plot
        zmap_out_Rv{i_part,ii} = witpc_z_Rv; %save values of plot
        
        % Z-value to p %(p_z method)
        p_witpc_z_Liv{i_part,ii} = 2*(1-normcdf(abs(witpc_z_Liv)));  %times 2 for 2-tailed test
        p_witpc_z_Riv{i_part,ii} = 2*(1-normcdf(abs(witpc_z_Riv)));  %times 2 for 2-tailed test
        p_witpc_z_Lv{i_part,ii} = 2*(1-normcdf(abs(witpc_z_Lv)));  %times 2 for 2-tailed test
        p_witpc_z_Rv{i_part,ii} = 2*(1-normcdf(abs(witpc_z_Rv)));  %times 2 for 2-tailed test
        
        % just wITPC
        witpc_threshold_out_Liv{i_part,ii} = prctile(max_pixel_vals_Liv(:),100-mcc_voxel_pval*100);
        witpc_threshold_out_Riv{i_part,ii} = prctile(max_pixel_vals_Riv(:),100-mcc_voxel_pval*100);
        witpc_threshold_out_Lv{i_part,ii} = prctile(max_pixel_vals_Lv(:),100-mcc_voxel_pval*100);
        witpc_threshold_out_Rv{i_part,ii} = prctile(max_pixel_vals_Rv(:),100-mcc_voxel_pval*100);
        
        zmapthresh_Liv = witpc_z_Liv;
        zmapthresh_Riv = witpc_z_Riv;
        zmapthresh_Lv = witpc_z_Lv;
        zmapthresh_Rv = witpc_z_Rv;
        
        zmapthresh_Liv(squeeze(squeeze(witpc_Liv(i_part,ii,:,:)))<witpc_threshold_out_Liv{i_part,ii})=0;
        zmapthresh_Riv(squeeze(squeeze(witpc_Riv(i_part,ii,:,:)))<witpc_threshold_out_Riv{i_part,ii})=0;
        zmapthresh_Lv(squeeze(squeeze(witpc_Lv(i_part,ii,:,:)))<witpc_threshold_out_Lv{i_part,ii})=0;
        zmapthresh_Rv(squeeze(squeeze(witpc_Rv(i_part,ii,:,:)))<witpc_threshold_out_Rv{i_part,ii})=0;
        
        zmap_witpc_thresh_Liv{i_part,ii} = zmapthresh_Liv; %save values of plot
        zmap_witpc_thresh_Riv{i_part,ii} = zmapthresh_Riv; %save values of plot
        zmap_witpc_thresh_Lv{i_part,ii} = zmapthresh_Lv; %save values of plot
        zmap_witpc_thresh_Rv{i_part,ii} = zmapthresh_Rv; %save values of plot
        
        %p-values after correction (p_z method)
        p_witpc_thresh_Liv{i_part,ii} = 2*(1-normcdf(abs(zmapthresh_Liv)));  %times 2 for 2-tailed test
        p_witpc_thresh_Riv{i_part,ii} = 2*(1-normcdf(abs(zmapthresh_Riv)));  %times 2 for 2-tailed test
        p_witpc_thresh_Lv{i_part,ii} = 2*(1-normcdf(abs(zmapthresh_Lv)));  %times 2 for 2-tailed test
        p_witpc_thresh_Rv{i_part,ii} = 2*(1-normcdf(abs(zmapthresh_Rv)));  %times 2 for 2-tailed test
        
        
        clear max_pixel_vals_* upper_threshold_* zmapthresh_*
        
        
        %% Create pixel-level corrected threshold
        max_pixel_vals_Liv = NaN(n_permutes,2); % initialize null hypothesis matrices
        max_pixel_vals_Riv = NaN(n_permutes,2); % initialize null hypothesis matrices
        max_pixel_vals_Lv = NaN(n_permutes,2); % initialize null hypothesis matrices
        max_pixel_vals_Rv = NaN(n_permutes,2); % initialize null hypothesis matrices
        for permi = 1:n_permutes
            
            % null z-map
            fakez_Liv = (squeeze(perm_witpc_Liv(permi,:,:))-permmean_Liv) ./ permstd_Liv;
            fakez_Riv = (squeeze(perm_witpc_Riv(permi,:,:))-permmean_Riv) ./ permstd_Riv;
            fakez_Lv = (squeeze(perm_witpc_Lv(permi,:,:))-permmean_Lv) ./ permstd_Lv;
            fakez_Rv = (squeeze(perm_witpc_Rv(permi,:,:))-permmean_Rv) ./ permstd_Rv;
            
            % save maximum pixel values
            max_pixel_vals_Liv(permi,:) = [ min(fakez_Liv(:)) max(fakez_Liv(:)) ];
            max_pixel_vals_Riv(permi,:) = [ min(fakez_Riv(:)) max(fakez_Riv(:)) ];
            max_pixel_vals_Lv(permi,:) = [ min(fakez_Lv(:)) max(fakez_Lv(:)) ];
            max_pixel_vals_Rv(permi,:) = [ min(fakez_Rv(:)) max(fakez_Rv(:)) ];

            clear fakez_*
        end
        clear permi permmean_* permstd_*
        
        % apply pixel-level corrected threshold
        lower_threshold_Liv = prctile(max_pixel_vals_Liv(:,1), mcc_voxel_pval*100/2);
        lower_threshold_Riv = prctile(max_pixel_vals_Riv(:,1), mcc_voxel_pval*100/2);
        lower_threshold_Lv = prctile(max_pixel_vals_Lv(:,1), mcc_voxel_pval*100/2);
        lower_threshold_Rv = prctile(max_pixel_vals_Rv(:,1), mcc_voxel_pval*100/2);
        
        upper_threshold_Liv = prctile(max_pixel_vals_Liv(:,2),100-mcc_voxel_pval*100/2);
        upper_threshold_Riv = prctile(max_pixel_vals_Riv(:,2),100-mcc_voxel_pval*100/2);
        upper_threshold_Lv = prctile(max_pixel_vals_Lv(:,2),100-mcc_voxel_pval*100/2);
        upper_threshold_Rv = prctile(max_pixel_vals_Rv(:,2),100-mcc_voxel_pval*100/2);
        
        threshold_out_Liv{i_part,ii} = [lower_threshold_Liv upper_threshold_Liv];
        threshold_out_Riv{i_part,ii} = [lower_threshold_Riv upper_threshold_Riv];
        threshold_out_Lv{i_part,ii} = [lower_threshold_Lv upper_threshold_Lv];
        threshold_out_Rv{i_part,ii} = [lower_threshold_Rv upper_threshold_Rv];
        
        zmapthresh_Liv = witpc_z_Liv;
        zmapthresh_Riv = witpc_z_Riv;
        zmapthresh_Lv = witpc_z_Lv;
        zmapthresh_Rv = witpc_z_Rv;
        
        zmapthresh_Liv(witpc_z_Liv>lower_threshold_Liv & witpc_z_Liv<upper_threshold_Liv)=0;
        zmapthresh_Riv(witpc_z_Riv>lower_threshold_Riv & witpc_z_Riv<upper_threshold_Riv)=0;
        zmapthresh_Lv(witpc_z_Lv>lower_threshold_Lv & witpc_z_Lv<upper_threshold_Lv)=0;
        zmapthresh_Rv(witpc_z_Rv>lower_threshold_Rv & witpc_z_Rv<upper_threshold_Rv)=0;
        
        %p-values after correction (p_z method)
        p_witpcz_thresh_Liv{i_part,ii} = 2*(1-normcdf(abs(zmapthresh_Liv)));  %times 2 for 2-tailed test
        p_witpcz_thresh_Riv{i_part,ii} = 2*(1-normcdf(abs(zmapthresh_Riv)));  %times 2 for 2-tailed test
        p_witpcz_thresh_Lv{i_part,ii} = 2*(1-normcdf(abs(zmapthresh_Lv)));  %times 2 for 2-tailed test
        p_witpcz_thresh_Rv{i_part,ii} = 2*(1-normcdf(abs(zmapthresh_Rv)));  %times 2 for 2-tailed test
        
        
        clear perm_witpc_* max_pixel_vals_* lower_threshold_* upper_threshold_*...
            zmapthresh_* i_elect witpc_z_*
    end
    clear ii n_resp_* tmp_err_*
end
clear i_part voxel_pval mcc_cluster_pval mcc_voxel_pval n_permutes


% /////////////////////////////////////////////////////////////////////////
%% Save Data
% this is large file so it will take some time to save
% save(['wITPCz_abs_err_post_' exp2.settings '.mat'],'p_*','threshold_*','witpc_*',...
%     'zmap_*','freqs','times','timewin','timelim','-v7.3')
% save(['wITPCz_M_post_' exp2.settings '.mat'],'p_*','threshold_*','witpc_*',...
%     'zmap_*','freqs','times','timewin','timelim','-v7.3')
save(['wITPCz_err_post_' exp2.settings '.mat'],'p_*','threshold_*','witpc_*',...
    'zmap_*','freqs','times','timewin','timelim','-v7.3')

% /////////////////////////////////////////////////////////////////////////




% /////////////////////////////////////////////////////////////////////////
%% Combine pvals across subjects
% /////////////////////////////////////////////////////////////////////////
%Stouffer
% result = squeeze(1-normcdf(sum(norminv(1-pmatrix),dim)./sqrt(size(pmatrix,dim))));
% .........................................................................

clear pcorr_cat_* pcorr_cat_avg_* pcorr_cat

% .........................................................................

% tmp_pcorr_Liv = p_witpcz_thresh_Liv; %pvals to combine
% tmp_pcorr_Riv = p_witpcz_thresh_Riv; %pvals to combine
% tmp_pcorr_Lv  = p_witpcz_thresh_Lv; %pvals to combine
% tmp_pcorr_Rv  = p_witpcz_thresh_Rv; %pvals to combine

% tmp_pcorr_Liv = p_witpc_thresh_Liv; %pvals to combine
% tmp_pcorr_Riv = p_witpc_thresh_Riv; %pvals to combine
% tmp_pcorr_Lv  = p_witpc_thresh_Lv; %pvals to combine
% tmp_pcorr_Rv  = p_witpc_thresh_Rv; %pvals to combine

tmp_pcorr_Liv = p_witpc_z_Liv; %pvals to combine (used in paper)
tmp_pcorr_Riv = p_witpc_z_Riv; %pvals to combine (used in paper)
tmp_pcorr_Lv  = p_witpc_z_Lv; %pvals to combine (used in paper)
tmp_pcorr_Rv  = p_witpc_z_Rv; %pvals to combine (used in paper)

pcorr_cat_Liv = NaN([length(exp2.singletrialselecs),length(freqs),length(timelim)]); %pre-allocate
pcorr_cat_Riv = NaN([length(exp2.singletrialselecs),length(freqs),length(timelim)]); %pre-allocate
pcorr_cat_Lv  = NaN([length(exp2.singletrialselecs),length(freqs),length(timelim)]); %pre-allocate
pcorr_cat_Rv  = NaN([length(exp2.singletrialselecs),length(freqs),length(timelim)]); %pre-allocate

% zmap_cat_Liv = NaN([length(exp2.singletrialselecs),length(freqs),length(times)]); %pre-allocate
for ii = 1:length(exp2.singletrialselecs)    
    i_elect = exp2.singletrialselecs(ii); %for doing only a selection of electrodes
    
    combo_p_z_Liv = NaN([length(exp2.participants),length(freqs),length(timelim)]); %pre-allocate
    combo_p_z_Riv = NaN([length(exp2.participants),length(freqs),length(timelim)]); %pre-allocate
    combo_p_z_Lv  = NaN([length(exp2.participants),length(freqs),length(timelim)]); %pre-allocate
    combo_p_z_Rv  = NaN([length(exp2.participants),length(freqs),length(timelim)]); %pre-allocate
    
    combo_z_Liv = NaN([length(exp2.participants),length(freqs),length(timelim)]); %pre-allocate
    combo_z_Riv = NaN([length(exp2.participants),length(freqs),length(timelim)]); %pre-allocate
    combo_z_Lv  = NaN([length(exp2.participants),length(freqs),length(timelim)]); %pre-allocate
    combo_z_Rv  = NaN([length(exp2.participants),length(freqs),length(timelim)]); %pre-allocate
    
    %extract each participant's data
    for i_part = 1:length(exp2.participants)
        
        combo_p_z_Liv(i_part,:,:) = tmp_pcorr_Liv{i_part,ii};
        combo_p_z_Riv(i_part,:,:) = tmp_pcorr_Riv{i_part,ii};
        combo_p_z_Lv(i_part,:,:)  = tmp_pcorr_Lv{i_part,ii};
        combo_p_z_Rv(i_part,:,:)  = tmp_pcorr_Rv{i_part,ii};
        
        combo_z_Liv(i_part,:,:) = zmap_out_Liv{i_part,ii};
        combo_z_Riv(i_part,:,:) = zmap_out_Riv{i_part,ii};
        combo_z_Lv(i_part,:,:)  = zmap_out_Lv{i_part,ii};
        combo_z_Rv(i_part,:,:)  = zmap_out_Rv{i_part,ii};
        
    end
    
    % function that combines p-values across participants for each electrode separately
    pcorr_cat_Liv(i_elect,:,:) = combine_pvalues(combo_p_z_Liv,1,1);
    pcorr_cat_Riv(i_elect,:,:) = combine_pvalues(combo_p_z_Riv,1,1);
    pcorr_cat_Lv(i_elect,:,:)  = combine_pvalues(combo_p_z_Lv,1,1);
    pcorr_cat_Rv(i_elect,:,:)  = combine_pvalues(combo_p_z_Rv,1,1);
    
    %trying combine pvalue procedure with z-values
%     zmap_cat_Liv(i_elect,:,:) = squeeze(sum(norminv(1-(1-normcdf(combo_z_Liv))),1)./sqrt(size(combo_z_Liv,1)));
    
    clear combo_p_z_* i_part combo_z_*   
end
clear ii i_elect tmp_pcorr_*

% Put into cell for plotting
pcorr_cat{1} = pcorr_cat_Liv;
pcorr_cat{2} = pcorr_cat_Riv;
pcorr_cat{3} = pcorr_cat_Lv;
pcorr_cat{4} = pcorr_cat_Rv;


% -------------------------------------------------------------------------
%% Calculate grand average across electrodes

tmp_out_Liv = pcorr_cat_Liv(2:32,:,:);
tmp_out_Riv = pcorr_cat_Riv(2:32,:,:);
tmp_out_Lv  = pcorr_cat_Lv(2:32,:,:);
tmp_out_Rv  = pcorr_cat_Rv(2:32,:,:);

% pcorr_cat_avg = combine_pvalues(tmp_out,1,1); %grand average

pcorr_cat_avg_Liv = squeeze(mean(tmp_out_Liv,1)); %grand average (used in paper)
pcorr_cat_avg_Riv = squeeze(mean(tmp_out_Riv,1)); %grand average (used in paper)
pcorr_cat_avg_Lv  = squeeze(mean(tmp_out_Lv,1)); %grand average (used in paper)
pcorr_cat_avg_Rv  = squeeze(mean(tmp_out_Rv,1)); %grand average (used in paper)

clear tmp_out_*



% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------
%% Plots
% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////

% Location to save results
saveFig = [pwd '\Figures\witpcz_err2_post_' exp2.settings '\']; % set save directory of data set


% if folder doesn't exist yet, create one
if ~exist(saveFig)
    mkdir(saveFig);
end

% /////////////////////////////////////////////////////////////////////////

alpha = 0.05;
pcorr_tmp = pcorr_cat;
cmap = makeColorMap([0.1098 0.5216 0.1922],[.68 .83 .17],[1 1 1],80);
for ii = 1:length(exp2.singletrialselecs)
    i_elect = exp2.singletrialselecs(ii); %for doing only a selection of electrodes
    
    figure; colormap(cmap); %open a new figure
    CLim = [0 0.05]; %cuz plotting p-vals

    subplot(2,2,1)
    tmp_plot = squeeze(pcorr_tmp{1}(i_elect,:,:)); %extract data
    [~,~,~,tmp_plot] = fdr_bh(tmp_plot,alpha,'pdep','yes'); %multi-comp correction
    imagesc(times(timelim),freqs,tmp_plot,CLim)
    set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
    line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
    xlim([0 600]); xticks(0:200:600)
    ylim([1 40]); yticks(5:5:40)
    ylabel('Freqency (Hz)'); xlabel('Time (ms)');
    t = colorbar('peer',gca,'Ticks',[0:.01:max(CLim)]);
    set(get(t,'ylabel'),'String', 'p-value');
    title(['Left, No Info: ' exp2.singtrlelec_name{ii}],'FontSize',10.5);
    clear tmp_plot
    
    subplot(2,2,2)
    tmp_plot = squeeze(pcorr_tmp{2}(i_elect,:,:)); %extract data
    [~,~,~,tmp_plot] = fdr_bh(tmp_plot,alpha,'pdep','yes'); %multi-comp correction
    imagesc(times(timelim),freqs,tmp_plot,CLim)
    set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
    line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
    xlim([0 600]); xticks(0:200:600)
    ylim([1 40]); yticks(5:5:40)
    ylabel('Freqency (Hz)'); xlabel('Time (ms)');
    t = colorbar('peer',gca,'Ticks',[0:.01:max(CLim)]);
    set(get(t,'ylabel'),'String', 'p-value');
    title(['Right, No Info: ' exp2.singtrlelec_name{ii}],'FontSize',10.5);
    clear tmp_plot
    
    subplot(2,2,3)
    tmp_plot = squeeze(pcorr_tmp{3}(i_elect,:,:)); %extract data
    [~,~,~,tmp_plot] = fdr_bh(tmp_plot,alpha,'pdep','yes'); %multi-comp correction
    imagesc(times(timelim),freqs,tmp_plot,CLim)
    set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
    line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
    xlim([0 600]); xticks(0:200:600)
    ylim([1 40]); yticks(5:5:40)
    ylabel('Freqency (Hz)'); xlabel('Time (ms)');
    t = colorbar('peer',gca,'Ticks',[0:.01:max(CLim)]);
    set(get(t,'ylabel'),'String', 'p-value');
    title(['Left, Info: ' exp2.singtrlelec_name{ii}],'FontSize',10.5);
    clear tmp_plot
    
    subplot(2,2,4)
    tmp_plot = squeeze(pcorr_tmp{4}(i_elect,:,:)); %extract data
    [~,~,~,tmp_plot] = fdr_bh(tmp_plot,alpha,'pdep','yes'); %multi-comp correction
    imagesc(times(timelim),freqs,tmp_plot,CLim)
    set(gca,'Ydir','Normal')
    line([0 0],[min(freqs) max(freqs)],'Color','k','LineStyle','--','LineWidth',1.5) %vertical line
    line([567 567],[min(freqs) max(freqs)],'color','m','LineStyle','--','LineWidth',1.5)  %vertical line for response screen onset
    xlim([0 600]); xticks(0:200:600)
    ylim([1 40]); yticks(5:5:40)
    ylabel('Freqency (Hz)'); xlabel('Time (ms)');
    t = colorbar('peer',gca,'Ticks',[0:.01:max(CLim)]);
    set(get(t,'ylabel'),'String', 'p-value');
    title(['Right, Info: ' exp2.singtrlelec_name{ii}],'FontSize',10.5);
    clear tmp_plot
    
    
    savefig([saveFig 'p_witpcz_05_' exp2.singtrlelec_name{ii}])

   clear i_elect i_cond
end
clear ii ncond conds CLim cmap POS_tmp alpha
























