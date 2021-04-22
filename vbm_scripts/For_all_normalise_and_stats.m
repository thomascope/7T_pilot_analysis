% List of open inputs
% Normalise to MNI Space: Dartel Template - cfg_files
% Normalise to MNI Space: Flow fields - cfg_files
% Normalise to MNI Space: Images - cfg_repeat
% Factorial design specification: Group 1 scans - cfg_files
% Factorial design specification: Group 2 scans - cfg_files
nrun = X; % enter the number of runs here
jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/For_all_normalise_and_stats_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(5, nrun);
for crun = 1:nrun
    inputs{1, crun} = MATLAB_CODE_TO_FILL_INPUT; % Normalise to MNI Space: Dartel Template - cfg_files
    inputs{2, crun} = MATLAB_CODE_TO_FILL_INPUT; % Normalise to MNI Space: Flow fields - cfg_files
    inputs{3, crun} = MATLAB_CODE_TO_FILL_INPUT; % Normalise to MNI Space: Images - cfg_repeat
    inputs{4, crun} = MATLAB_CODE_TO_FILL_INPUT; % Factorial design specification: Group 1 scans - cfg_files
    inputs{5, crun} = MATLAB_CODE_TO_FILL_INPUT; % Factorial design specification: Group 2 scans - cfg_files
end
spm('defaults', 'PET');
spm_jobman('run', jobs, inputs{:});
