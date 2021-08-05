function covariatesecondlevelworkedcorrectly = module_multivariate_all_correlations(covariate,GLMDir,subjects,group,age_lookup,outpath,downsamp_ratio)
%Input a covariate, and the path to the first level scan of the first
%subject, assuming all the others have the same path simply with the
%subject ID replaced. Corrects for age.

%Check covariate is a column
if size(covariate,1)==1
    covariate = covariate';
end

assert(size(covariate,1)==length(subjects),'Number of covariate measures and subject IDs does not match')

do_smoothed_maps = 0;  % if you want to do smoothed maps change this
if ~exist('downsamp_ratio','var')
    downsamp_ratio = 1;
end

versionCurrent = 'spearman';

% Gather images for current subject
if downsamp_ratio == 1
    images = cellstr(spm_select('FPList', [GLMDir '/TDTcrossnobis/' versionCurrent '/'], '^whireseffect-map_.*.nii'));
else
    images = cellstr(spm_select('FPList', [GLMDir '/TDTcrossnobis_downsamp_' num2str(downsamp_ratio) '/' versionCurrent '/'], '^whireseffect-map_.*.nii'));
end

nrun = size(images,1); % enter the number of runs here - if want to do smoothed as well
%jobfile = {'/group/language/data/thomascope/vespa/SPM12version/Standalone preprocessing pipeline/tc_source/batch_forwardmodel_job_noheadpoints.m'};

jobfile = {'/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/module_covariate_secondlevel_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(4, nrun);

for this_condition = 1:(nrun)
    this_scan = {};
    remove_these = [];
    this_age = [];
    
    condition_name = strsplit(images{this_condition},'whireseffect-map_');
    condition_name = condition_name{2}(1:end-4);
    
    inputs{1, this_condition} = cellstr([outpath filesep condition_name '_hires']);
    
    for crun = 1:size(subjects,2)
        if isnan(covariate(crun))
            disp(['Warning, covariate is not a number for ' subjects{crun} ' so removing them - please check'])
            continue
        else
            this_age(end+1) = age_lookup.Age(strcmp(age_lookup.Study_ID,subjects{crun}));
            this_scan(end+1) = cellstr(strrep(images{this_condition},subjects{1},subjects{crun}));
        end
    end
    
    inputs{2, this_condition} = this_scan'; % Factorial design specification: Scans - cfg_files
    inputs{3, this_condition} = covariate(~isnan(covariate)); % Factorial design specification: Vector - cfg_entry
    inputs{4, this_condition} = this_age'; % Factorial design specification: Age - cfg_entry
end

covariatesecondlevelworkedcorrectly = zeros(1,nrun);

parfor crun = 1:nrun
    if exist(fullfile(char(inputs{1, crun}),'SPM.mat'),'file')
        disp([fullfile(char(inputs{1, crun}),'SPM.mat') ' already exists, delete it if you want to re-make, otherwise moving on.'])
    elseif ~all(cellfun(@exist,inputs{2, crun}))
        disp(['Missing input files for ' char(inputs{1, crun}) ' moving on'])
    else
        spm('defaults', 'fMRI');
        spm_jobman('initcfg')
        try
            spm_jobman('run', jobs{crun}, inputs{:,crun});
            covariatesecondlevelworkedcorrectly(crun) = 1;
        catch
            covariatesecondlevelworkedcorrectly(crun) = 0;
        end
    end
end