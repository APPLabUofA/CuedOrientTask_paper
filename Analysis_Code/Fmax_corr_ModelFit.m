%Calculate Fmax correction for a F-observed and the permutation
%distribution
%
%EXAMPLE USAGE
% >> [h, p, Fmax_crit, est_alpha] = Fmax_corr(F_obs, F_dist, 0.05)
%
%REQUIRED INPUTS
% F_obs          - An electrode x time points array of observed F-values
% F-dist         - A permutation x electrode x time point array of the
%                  permutation F-distribution
% alpha          - Alpha level to use for hypothesis test
%
%OUTPUT
% h              - electrode x time point array indicating which locations
%                  are part of a statistically significant cluster
% p              - electrode x time point array of p-values
% Fmax_crit      - Critical value for statistical significance
% est_alpha      - estimated achieved alpha level of the test; may not be
%                  accurate if F_dist was generated by an approximate
%                  permutation method
%
%
%VERSION DATE: 4 April 2019
%AUTHOR: Eric Fields
%
%NOTE: This function is provided "as is" and any express or implied warranties 
%are disclaimed. 

%Copyright (c) 2017-2019, Eric Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.
% 
% Adapted by SSS
% Date: April 2021

function [h, p, Fmax_crit, est_alpha] = Fmax_corr_ModelFit(F_obs, F_dist, alpha)
    
    global VERBLEVEL;
    
    %Some useful numbers
    [~, n_electrodes] = size(F_dist);
    
    %Calculate Fmax distribution and Fmax critical value
    Fmax_dist = max(F_dist, [], 2);
    Fmax_dist = sort(Fmax_dist);
    Fmax_crit = Fmax_dist(ceil((1-alpha) * length(Fmax_dist)));

    %Null hypothesis test
    h = F_obs > Fmax_crit;
    est_alpha = mean(Fmax_dist > Fmax_crit);
    if VERBLEVEL
        fprintf('Estimated alpha level is %f\n', est_alpha);
    end

    %Calculate p-value
    p = NaN(n_electrodes);
    for e = 1:n_electrodes
        p(e) = mean(Fmax_dist >= F_obs(e));
    end

    assert(isequal(h, p<=alpha));

end