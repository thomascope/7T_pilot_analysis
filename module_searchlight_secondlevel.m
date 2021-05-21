function secondlevelworkedcorrectly = module_searchlight_secondlevel(GLMDir,subjects,group,age_lookup,outpath,downsamp_ratio)
% Normalise effect-maps to MNI template

do_smoothed_maps = 0;  % if you want to do smoothed maps change this
if ~exist('downsamp_ratio','var')
    downsamp_ratio = 1;
end

versionCurrent = 'spearman';

% Gather images for current subject
if downsamp_ratio == 1
    images = cellstr(spm_select('FPList', [GLMDir '/TDTcrossnobis/' versionCurrent '/'], '^weffect-map_.*.nii'));
else
    images = cellstr(spm_select('FPList', [GLMDir '/TDTcrossnobis_downsamp_' num2str(downsamp_ratio) '/' versionCurrent '/'], '^weffect-map_.*.nii'));
end

nrun = 2*size(images,1); % enter the number of runs here - if want to do smoothed as well
%jobfile = {'/group/language/data/thomascope/vespa/SPM12version/Standalone preprocessing pipeline/tc_source/batch_forwardmodel_job_noheadpoints.m'};

this_scan = {};
this_t_scan = {};

jobfile = {'/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/module_secondlevel_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(4, nrun);

for this_condition = 1:(nrun/2)
    group1_mrilist = {}; %NB: Patient MRIs, so here group 2 (sorry)
    group1_ages = [];
    group2_mrilist = {};
    group2_ages = [];
    
    condition_name = strsplit(images{this_condition},'weffect-map_');
    condition_name = condition_name{2}(1:end-4);
    
    inputs{1, this_condition} = cellstr([outpath filesep condition_name]);
    
    for crun = 1:size(subjects,2)
        this_age = age_lookup.Age(strcmp(age_lookup.Study_ID,subjects{crun}));
        this_scan(crun) = cellstr(strrep(images{this_condition},subjects{1},subjects{crun}));
        
        if group(crun) == 1 % Controls
            group2_mrilist(end+1) = this_scan(crun);
            group2_ages(end+1) = this_age;
        elseif group(crun) == 2 % Patients
            group1_mrilist(end+1) = this_scan(crun);
            group1_ages(end+1) = this_age;
        end
    end
    inputs{2, this_condition} = group1_mrilist';
    inputs{3, this_condition} = group2_mrilist';
    inputs{4, this_condition} = [group1_ages';group2_ages'];
end

%% if you want to do smoothed maps uncomment this

if ~do_smoothed_maps
    nrun = nrun/2;
else
    for this_condition = (1+(nrun/2)):nrun
        group1_mrilist = {}; %NB: Patient MRIs, so here group 2 (sorry)
        group1_ages = [];
        group2_mrilist = {};
        group2_ages = [];
        
        images{this_condition-(nrun/2)} = strrep(images{this_condition-(nrun/2)},'weffect-map_','sweffect-map_');
        
        condition_name = strsplit(images{this_condition-(nrun/2)},'sweffect-map_');
        condition_name = condition_name{2}(1:end-4);
        
        inputs{1, this_condition} = cellstr([outpath filesep 'sm_' condition_name]);
        
        for crun = 1:size(subjects,2)
            this_age = age_lookup.Age(strcmp(age_lookup.Study_ID,subjects{crun}));
            this_scan(crun) = cellstr(strrep(images{this_condition-(nrun/2)},subjects{1},subjects{crun}));
            
            if group(crun) == 1 % Controls
                group2_mrilist(end+1) = this_scan(crun);
                group2_ages(end+1) = this_age;
            elseif group(crun) == 2 % Patients
                group1_mrilist(end+1) = this_scan(crun);
                group1_ages(end+1) = this_age;
            end
        end
        inputs{2, this_condition} = group1_mrilist';
        inputs{3, this_condition} = group2_mrilist';
        inputs{4, this_condition} = [group1_ages';group2_ages'];
    end
end

secondlevelworkedcorrectly = zeros(1,nrun);

parfor crun = 1:nrun
    if exist(fullfile(char(inputs{1, crun}),'SPM.mat'),'file')
        disp([fullfile(char(inputs{1, crun}),'SPM.mat') ' already exists, delete it if you want to re-make, otherwise moving on.'])
    elseif ~all(cellfun(@exist,inputs{2, crun})) || ~all(cellfun(@exist,inputs{3, crun}))
        disp(['Missing input files for ' char(inputs{1, crun}) ' moving on'])
    else
        spm('defaults', 'fMRI');
        spm_jobman('initcfg')
        try
            spm_jobman('run', jobs{crun}, inputs{:,crun});
            secondlevelworkedcorrectly(crun) = 1;
        catch
            secondlevelworkedcorrectly(crun) = 0;
        end
    end
end