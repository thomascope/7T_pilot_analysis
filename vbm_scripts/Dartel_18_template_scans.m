% List of open inputs
% Segment: Volumes - cfg_files
nrun = X; % enter the number of runs here
jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/Dartel_18_template_scans_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(1, nrun);
for crun = 1:nrun
    inputs{1, crun} = MATLAB_CODE_TO_FILL_INPUT; % Segment: Volumes - cfg_files
end
spm('defaults', 'PET');
spm_jobman('run', jobs, inputs{:});
