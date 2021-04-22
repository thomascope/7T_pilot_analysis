% List of open inputs
% Factorial design specification: Directory - cfg_files
% Factorial design specification: Scans - cfg_files
% Factorial design specification: Scans - cfg_files
% Factorial design specification: Scans - cfg_files
% Factorial design specification: Scans - cfg_files
% Factorial design specification: Vector - cfg_entry
% Factorial design specification: Vector - cfg_entry
% Factorial design specification: Explicit Mask - cfg_files
nrun = X; % enter the number of runs here
jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/scratchpad_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(8, nrun);
for crun = 1:nrun
    inputs{1, crun} = MATLAB_CODE_TO_FILL_INPUT; % Factorial design specification: Directory - cfg_files
    inputs{2, crun} = MATLAB_CODE_TO_FILL_INPUT; % Factorial design specification: Scans - cfg_files
    inputs{3, crun} = MATLAB_CODE_TO_FILL_INPUT; % Factorial design specification: Scans - cfg_files
    inputs{4, crun} = MATLAB_CODE_TO_FILL_INPUT; % Factorial design specification: Scans - cfg_files
    inputs{5, crun} = MATLAB_CODE_TO_FILL_INPUT; % Factorial design specification: Scans - cfg_files
    inputs{6, crun} = MATLAB_CODE_TO_FILL_INPUT; % Factorial design specification: Vector - cfg_entry
    inputs{7, crun} = MATLAB_CODE_TO_FILL_INPUT; % Factorial design specification: Vector - cfg_entry
    inputs{8, crun} = MATLAB_CODE_TO_FILL_INPUT; % Factorial design specification: Explicit Mask - cfg_files
end
spm('defaults', 'PET');
spm_jobman('run', jobs, inputs{:});
