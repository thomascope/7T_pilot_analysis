% List of open inputs
% Estimate TIV and global tissue volumes: XML files - cfg_files
% Estimate TIV and global tissue volumes: Output file - cfg_entry
nrun = X; % enter the number of runs here
jobfile = {'/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/module_cat12_tiv_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(2, nrun);
for crun = 1:nrun
    inputs{1, crun} = MATLAB_CODE_TO_FILL_INPUT; % Estimate TIV and global tissue volumes: XML files - cfg_files
    inputs{2, crun} = MATLAB_CODE_TO_FILL_INPUT; % Estimate TIV and global tissue volumes: Output file - cfg_entry
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});
