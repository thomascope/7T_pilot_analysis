% List of open inputs
% Factorial design specification: Scans - cfg_files
nrun = X; % enter the number of runs here
jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_factorial_TIV_age_ERF_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(1, nrun);
for crun = 1:nrun
    inputs{1, crun} = MATLAB_CODE_TO_FILL_INPUT; % Factorial design specification: Scans - cfg_files
end
spm('defaults', 'PET');
spm_jobman('run', jobs, inputs{:});
