% List of open inputs
% Smooth: Images to Smooth - cfg_files
% Smooth: Images to Smooth - cfg_files
nrun = X; % enter the number of runs here
jobfile = {'/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/module_smooth_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(2, nrun);
for crun = 1:nrun
    inputs{1, crun} = MATLAB_CODE_TO_FILL_INPUT; % Smooth: Images to Smooth - cfg_files
    inputs{2, crun} = MATLAB_CODE_TO_FILL_INPUT; % Smooth: Images to Smooth - cfg_files
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});
