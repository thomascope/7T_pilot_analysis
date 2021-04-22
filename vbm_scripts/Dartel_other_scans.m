% List of open inputs
% Segment: Volumes - cfg_files
% Run Dartel (existing Templates): Template - cfg_files
% Run Dartel (existing Templates): Template - cfg_files
% Run Dartel (existing Templates): Template - cfg_files
% Run Dartel (existing Templates): Template - cfg_files
% Run Dartel (existing Templates): Template - cfg_files
% Run Dartel (existing Templates): Template - cfg_files
nrun = X; % enter the number of runs here
jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/Dartel_other_scans_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(7, nrun);
for crun = 1:nrun
    inputs{1, crun} = MATLAB_CODE_TO_FILL_INPUT; % Segment: Volumes - cfg_files
    inputs{2, crun} = MATLAB_CODE_TO_FILL_INPUT; % Run Dartel (existing Templates): Template - cfg_files
    inputs{3, crun} = MATLAB_CODE_TO_FILL_INPUT; % Run Dartel (existing Templates): Template - cfg_files
    inputs{4, crun} = MATLAB_CODE_TO_FILL_INPUT; % Run Dartel (existing Templates): Template - cfg_files
    inputs{5, crun} = MATLAB_CODE_TO_FILL_INPUT; % Run Dartel (existing Templates): Template - cfg_files
    inputs{6, crun} = MATLAB_CODE_TO_FILL_INPUT; % Run Dartel (existing Templates): Template - cfg_files
    inputs{7, crun} = MATLAB_CODE_TO_FILL_INPUT; % Run Dartel (existing Templates): Template - cfg_files
end
spm('defaults', 'PET');
spm_jobman('run', jobs, inputs{:});
