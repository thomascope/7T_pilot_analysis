% List of open inputs
nrun = 1; % enter the number of runs here
jobfile = {'/imaging/mlr/users/tc02/7T_full_paradigm_pilot_second_am/Preprocessed_Images/AM_univariate_basisderivatives_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
    for sess = 1:4
    inputs{(8*(sess-1))+1, crun} = cat(2, tempDesign{sess}{1:9})';
    inputs{(8*(sess-1))+2, crun} = cat(2, tempDesign{sess}{10:18})';
    inputs{(8*(sess-1))+3, crun} = cat(2, tempDesign{sess}{19:27})';
    inputs{(8*(sess-1))+4, crun} = cat(2, tempDesign{sess}{28:36})';
    inputs{(8*(sess-1))+5, crun} = cat(2, tempDesign{sess}{37:45})';
    inputs{(8*(sess-1))+6, crun} = cat(2, tempDesign{sess}{46:54})';
    inputs{(8*(sess-1))+7, crun} = cat(2, tempDesign{sess}{[55:63, 73]})';
    inputs{(8*(sess-1))+8, crun} = cat(2, tempDesign{sess}{[64:72, 74]})';
    end
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});
