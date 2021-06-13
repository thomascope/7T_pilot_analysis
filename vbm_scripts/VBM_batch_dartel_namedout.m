% List of open inputs
% Run Dartel (create Templates): Images - cfg_files
% Run Dartel (create Templates): Images - cfg_files
% Run Dartel (create Templates): Template basename - cfg_entry
nrun = X; % enter the number of runs here
jobfile = {'/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/vbm_scripts/VBM_batch_dartel_namedout_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(3, nrun);
for crun = 1:nrun
    inputs{1, crun} = MATLAB_CODE_TO_FILL_INPUT; % Run Dartel (create Templates): Images - cfg_files
    inputs{2, crun} = MATLAB_CODE_TO_FILL_INPUT; % Run Dartel (create Templates): Images - cfg_files
    inputs{3, crun} = MATLAB_CODE_TO_FILL_INPUT; % Run Dartel (create Templates): Template basename - cfg_entry
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});
