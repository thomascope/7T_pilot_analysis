% List of open inputs
% CAT12: Segmentation: Volumes - cfg_files
nrun = X; % enter the number of runs here
jobfile = {'/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/module_cat12_segment_expert_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(1, nrun);
for crun = 1:nrun
    inputs{1, crun} = MATLAB_CODE_TO_FILL_INPUT; % CAT12: Segmentation: Volumes - cfg_files
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});
