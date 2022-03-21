function [error_deg] = getBEHdata_OrientWheel_Exo(subj_id,exp)
% Used in the pre-processing code


% Load each subject's data
load([exp.pathname 'BEH\' subj_id '_Orient_Exo.mat'])
% load(['C:\Users\ssshe\Documents\MathLab\Data\OrientWheel_Exo\BEH\' subj_id '_Orient_Exo.mat']) %my laptop

% Remove practice trials (first 20 trials)
error_deg = data.errorDegrees(1,21:end);
valid = data.valid(1,21:end);
position = data.position(1,21:end);

% Remove trials where no target appeared
% error_deg = error_deg(targets == 1);
% error_rad = error_rad(targets == 1);
% report_deg = report_deg(targets == 1);
% report_rad = report_rad(targets == 1);


% Need to remove trials from subject due to recording error
if strcmpi(subj_id,'009')
    error_deg(344:end) = [];
end





















