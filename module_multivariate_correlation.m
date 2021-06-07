function covariatesecondlevelworkedcorrectly = module_multivariate_correlation(covariate,subject1scanofinterest,subjects,age_lookup,outpath)
%Input a covariate, and the path to the first level scan of the first
%subject, assuming all the others have the same path simply with the
%subject ID replaced. Corrects for age.

%Check covariate is a column
if size(covariate,1)==1
    covariate = covariate';
end

nrun = 1; % enter the number of runs here
jobfile = {'/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/module_covariate_secondlevel_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(4, nrun);
this_scan = {};
remove_these = [];
this_age = [];
for crun = 1:size(subjects,2)
    if isnan(covariate(crun))
        disp(['Warning, covariate is not a number for ' subjects{crun} ' so removing them - please check'])
        continue
    else
        this_age(end+1) = age_lookup.Age(strcmp(age_lookup.Study_ID,subjects{crun}));
        this_scan(end+1) = cellstr(strrep(subject1scanofinterest,subjects{1},subjects{crun}));
    end
end

inputs{1, 1} = cellstr(outpath); % Factorial design specification: Directory - cfg_files
inputs{2, 1} = this_scan'; % Factorial design specification: Scans - cfg_files
inputs{3, 1} = covariate(~isnan(covariate)); % Factorial design specification: Vector - cfg_entry
inputs{4, 1} = this_age'; % Factorial design specification: Age - cfg_entry

spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});

covariatesecondlevelworkedcorrectly = 1;

