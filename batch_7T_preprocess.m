% Batch script for preprocessing of pilot 7T data
% Written by TEC Feb 2018 and updated 2021
%
%% Setup environment
clear all
rmpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/'))
%addpath /imaging/local/software/spm_cbu_svn/releases/spm12_fil_r6906
addpath /group/language/data/thomascope/spm12_fil_r6906/
%spm fmri
scriptdir = '/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/';

%% Define parameters
setup_file = 'PINFA_subjects_parameters';
eval(setup_file)
tr=2.5;

%% Options to skip steps
applytopup = 1;
opennewanalysispool = 0;

%% Open a worker pool
if opennewanalysispool == 1
    if size(subjects,2) > 64
        workersrequested = 64;
        fprintf([ '\n\nUnable to ask for a worker per run; asking for 64 instead\n\n' ]);
    else
        workersrequested = size(subjects,2);
    end
    
    %Open a parallel pool
    if numel(gcp('nocreate')) == 0
        Poolinfo = cbupool(workersrequested,'--mem-per-cpu=12G --time=167:00:00 --exclude=node-i[01-15]');
        parpool(Poolinfo,Poolinfo.NumWorkers);
    end
end

%% Skullstrip structural
nrun = size(subjects,2); % enter the number of runs here
%jobfile = {'/group/language/data/thomascope/vespa/SPM12version/Standalone preprocessing pipeline/tc_source/batch_forwardmodel_job_noheadpoints.m'};
jobfile = {[scriptdir 'module_skullstrip_INV2_job.m']};
inputs = cell(2, nrun);

for crun = 1:nrun
    inputs{1, crun} = cellstr([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '/' blocksin{crun}{find(strcmp(blocksout{crun},'structural'))} ',1']);
    inputs{2, crun} = cellstr([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'INV2'))} '/' blocksin{crun}{find(strcmp(blocksout{crun},'INV2'))}]);
    inputs{3, crun} = cellstr([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '/' blocksin{crun}{find(strcmp(blocksout{crun},'structural'))}]);
    inputs{4, crun} = 'structural';
    inputs{5, crun} = cellstr([preprocessedpathstem subjects{crun} '/']);
    if ~exist(inputs{5, crun}{1})
        mkdir(inputs{5, crun}{1});
    end
    inputs{6, crun} = 'structural_csf';
    inputs{7, crun} = cellstr([preprocessedpathstem subjects{crun} '/']);
    inputs{8, crun} = 'c1structural';
    inputs{9, crun} = cellstr([preprocessedpathstem subjects{crun} '/']);
    inputs{10, crun} = 'c2structural';
    inputs{11, crun} = cellstr([preprocessedpathstem subjects{crun} '/']);
end

skullstripworkedcorrectly = zeros(1,nrun);
jobs = repmat(jobfile, 1, 1);

parfor crun = 1:nrun
    spm('defaults', 'fMRI');
    spm_jobman('initcfg')
    try
        spm_jobman('run', jobs, inputs{:,crun});
        skullstripworkedcorrectly(crun) = 1;
    catch
        skullstripworkedcorrectly(crun) = 0;
    end
end

if ~all(skullstripworkedcorrectly)
    error('failed at skullstrip');
end

%% Now realign the EPIs

realignworkedcorrectly = zeros(1,nrun);
parfor crun = 1:nrun
    base_image_path = [rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'Pos_topup'))} '/' blocksin{crun}{find(strcmp(blocksout{crun},'Pos_topup'))}];
    reversed_image_path = [rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'Neg_topup'))} '/' blocksin{crun}{find(strcmp(blocksout{crun},'Neg_topup'))}];
    theseepis = find(strncmp(blocksout{crun},'Run',3));
    filestorealign = cell(1,length(theseepis));
    for i = 1:length(theseepis)
        inpath = [rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{theseepis(i)} '/'];
        filestorealign{i} = spm_select('ExtFPList',inpath,['^' blocksin{crun}{theseepis(i)}],1:minvols(crun));
    end
    filestorealign{i+1} = base_image_path
    filestorealign{i+2} = reversed_image_path
    
    flags = struct;
    flags.fhwm = 3;
    flags.interp = 5; % Improve quality by using 5th degree B-spline interpolation
    try
        spm_realign(filestorealign,flags)
        realignworkedcorrectly(crun) = 1;
    catch
        realignworkedcorrectly(crun) = 0;
    end
    for i = 1:length(theseepis) %Now move movement parameters
        inpath = [rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{theseepis(i)} '/'];
        outpath = [preprocessedpathstem subjects{crun} '/'];
        copyfile([inpath 'rp_' blocksin{crun}{theseepis(i)}(1:end-4) '.txt'],[outpath 'rp_' blocksin{crun}{theseepis(i)}(1:end-4) '.txt'])
    end
end

if ~all(realignworkedcorrectly)
    error('failed at realign');
end

%% Now apply topup to distortion correct the EPI

if applytopup == 1
    topupworkedcorrectly = zeros(1,nrun);
    parfor crun = 1:nrun
        base_image_path = [rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'Pos_topup'))} '/' blocksin{crun}{find(strcmp(blocksout{crun},'Pos_topup'))}];
        reversed_image_path = [rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'Neg_topup'))} '/' blocksin{crun}{find(strcmp(blocksout{crun},'Neg_topup'))}];
        outpath = [preprocessedpathstem subjects{crun} '/'];
        theseepis = find(strncmp(blocksout{crun},'Run',3));
        filestocorrect = cell(1,length(theseepis));
        for i = 1:length(theseepis)
            filestocorrect{i} = [rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{theseepis(i)} '/' blocksin{crun}{theseepis(i)}];
        end
        
        try
            module_topup_job(base_image_path, reversed_image_path, outpath, minvols(crun), filestocorrect)
            topupworkedcorrectly(crun) = 1;
        catch
            topupworkedcorrectly(crun) = 0;
        end
        
    end
end

if ~all(topupworkedcorrectly)
    error('failed at topup');
end

% %% Now realign the EPIs % Now moved to before topup
%
% realignworkedcorrectly = zeros(1,nrun);
% parfor crun = 1:nrun
%     theseepis = find(strncmp(blocksout{crun},'Run',3))
%     filestorealign = cell(1,length(theseepis));
%     outpath = [preprocessedpathstem subjects{crun} '/'];
%     for i = 1:length(theseepis)
%         filestorealign{i} = spm_select('ExtFPList',outpath,['^topup_' blocksin{crun}{theseepis(i)}],1:minvols(crun));
%     end
%     flags = struct;
%     flags.fhwm = 3;
%     try
%         spm_realign(filestorealign,flags)
%         realignworkedcorrectly(crun) = 1;
%     catch
%         realignworkedcorrectly(crun) = 0;
%     end
% end
%
% if ~all(realignworkedcorrectly)
%     error('failed at realign');
% end

%% Now reslice all the images

resliceworkedcorrectly = zeros(1,nrun);
parfor crun = 1:nrun
    theseepis = find(strncmp(blocksout{crun},'Run',3))
    filestorealign = cell(1,length(theseepis));
    outpath = [preprocessedpathstem subjects{crun} '/'];
    for i = 1:length(theseepis)
        filestorealign{i} = spm_select('ExtFPList',outpath,['^topup_' blocksin{crun}{theseepis(i)}],1:minvols(crun));
    end
    flags = struct
    flags.which = 2;
    try
        spm_reslice(filestorealign,flags)
        resliceworkedcorrectly(crun) = 1;
    catch
        resliceworkedcorrectly(crun) = 0;
    end
end

if ~all(resliceworkedcorrectly)
    error('failed at reslice');
end

%% Now smooth the realigned, undistorted, resliced native space functionals at 3 and 8
nrun = size(subjects,2); % enter the number of runs here
%jobfile = {'/group/language/data/thomascope/vespa/SPM12version/Standalone preprocessing pipeline/tc_source/batch_forwardmodel_job_noheadpoints.m'};
jobfile = {[scriptdir 'module_smooth_job.m']};
inputs = cell(2, nrun);

for crun = 1:nrun
    outpath = [preprocessedpathstem subjects{crun} '/'];
    
    theseepis = find(strncmp(blocksout{crun},'Run',3));
    filestosmooth = cell(1,length(theseepis));
    filestosmooth_list = [];
    for i = 1:length(theseepis)
        filestosmooth{i} = spm_select('ExtFPList',outpath,['^rtopup_' blocksin{crun}{theseepis(i)}],1:minvols(crun));
        filestosmooth_list = [filestosmooth_list; filestosmooth{i}];
    end
    inputs{1, crun} = cellstr(filestosmooth_list); % Needs to be twice, once for each smoothing kernel
    inputs{2, crun} = cellstr(filestosmooth_list);
end

smoothworkedcorrectly = zeros(1,nrun);
jobs = repmat(jobfile, 1, 1);

parfor crun = 1:nrun
    spm('defaults', 'fMRI');
    spm_jobman('initcfg')
    try
        spm_jobman('run', jobs, inputs{:,crun});
        smoothworkedcorrectly(crun) = 1;
    catch
        smoothworkedcorrectly(crun) = 0;
    end
end

if ~all(smoothworkedcorrectly)
    error('failed at native space smooth');
end

%% Now co-register estimate, using structural as reference, mean as source and epi as others, then reslice only the mean

coregisterworkedcorrectly = zeros(1,nrun);

parfor crun = 1:nrun
    job = struct
    job.eoptions.cost_fun = 'nmi'
    job.eoptions.tol = [repmat(0.02,1,3), repmat(0.01,1,6), repmat(0.001,1,3)];
    job.eoptions.sep = [4 2];
    job.eoptions.fwhm = [7 7];
    
    outpath = [preprocessedpathstem subjects{crun} '/'];
    job.ref = {[outpath 'structural.nii,1']};
    theseepis = find(strncmp(blocksout{crun},'Run',3));
    job.source = {[outpath 'meantopup_' blocksin{crun}{theseepis(1)} ',1']};
    
    filestocoregister = cell(1,length(theseepis));
    filestocoregister_list = [];
    for i = 1:length(theseepis)
        filestocoregister{i} = spm_select('ExtFPList',outpath,['^topup_' blocksin{crun}{theseepis(i)}],1:minvols(crun));
        filestocoregister_list = [filestocoregister_list; filestocoregister{i}]
    end
    filestocoregister = cellstr(filestocoregister_list);
    
    job.other = filestocoregister
    
    try
        spm_run_coreg(job)
        
        % Now co-register reslice the mean EPI
        P = char(job.ref{:},job.source{:});
        spm_reslice(P)
        
        coregisterworkedcorrectly(crun) = 1;
    catch
        coregisterworkedcorrectly(crun) = 0;
    end
end

if ~all(coregisterworkedcorrectly)
    error('failed at coregister');
end


%% Now do cat12 normalisation of the structural to create deformation fields (works better than SPM segment deformation fields, which sometimes produce too-small brains)

nrun = size(subjects,2); % enter the number of runs here
jobfile = {'/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/module_cat12_normalise_job.m'};
inputs = cell(1, nrun);
for crun = 1:nrun
    outpath = [preprocessedpathstem subjects{crun} '/'];
    inputs{1, crun} = cellstr([outpath 'structural_csf.nii']);
end

cat12workedcorrectly = zeros(1,nrun);
jobs = repmat(jobfile, 1, 1);

parfor crun = 1:nrun
    spm('defaults', 'fMRI');
    spm_jobman('initcfg')
    try
        spm_jobman('run', jobs, inputs{:,crun});
        cat12workedcorrectly(crun) = 1;
    catch
        cat12workedcorrectly(crun) = 0;
    end
end

if ~all(cat12workedcorrectly)
    error('failed at cat12');
end


%% Now run a VBM on the structurals
age_lookup = readtable('Pinfa_ages.csv');
visual_check = 0;
nrun = size(subjects,2); % enter the number of runs here
%jobfile = {'/group/language/data/thomascope/vespa/SPM12version/Standalone preprocessing pipeline/tc_source/batch_forwardmodel_job_noheadpoints.m'};
group1_mrilist = {}; %NB: Patient MRIs, so here group 2 (sorry)
group1_ages = [];
group2_mrilist = {};
group2_ages = [];
this_scan = {};
this_segmented = {};
segmented = 1;

for crun = 1:nrun
    this_age = age_lookup.Age(strcmp(age_lookup.Study_ID,subjects{crun}));
    this_scan(crun) = cellstr([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '/' blocksin{crun}{find(strcmp(blocksout{crun},'structural'))}]);
    %this_scan(crun) = cellstr([preprocessedpathstem subjects{crun} '/structural.nii']);
    if segmented
        this_segmented(crun) = cellstr([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '/c1' blocksin{crun}{find(strcmp(blocksout{crun},'structural'))}]);
    end
    if group(crun) == 1 % Controls
        group2_mrilist(end+1) = this_scan(crun);
        group2_ages(end+1) = this_age;
    elseif group(crun) == 2 % Patients
        group1_mrilist(end+1) = this_scan(crun);
        group1_ages(end+1) = this_age;
    end
end
if visual_check
    if segmented
        spm_check_registration(this_segmented{:})
    else
        spm_check_registration(this_scan{:}) % Optional visual check of your input images (don't need to be aligned or anything, just to see they're all structurals and exist)
    end
end

module_vbm_job(group1_mrilist, group1_ages, group2_mrilist, group2_ages, preprocessedpathstem, segmented)

%% Now normalise write for visualisation and smooth at 3 and 8
nrun = size(subjects,2); % enter the number of runs here
%jobfile = {'/group/language/data/thomascope/vespa/SPM12version/Standalone preprocessing pipeline/tc_source/batch_forwardmodel_job_noheadpoints.m'};
jobfile = {[scriptdir 'module_normalise_smooth_job.m']};
inputs = cell(2, nrun);

for crun = 1:nrun
    outpath = [preprocessedpathstem subjects{crun} '/'];
    
    % % First is for SPM segment, second for CAT12
    %inputs{1, crun} = cellstr([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '/y_' blocksin{crun}{find(strcmp(blocksout{crun},'structural'))}]);
    inputs{1, crun} = cellstr([outpath 'mri/y_structural_csf.nii']);
    
    
    theseepis = find(strncmp(blocksout{crun},'Run',3));
    filestonormalise = cell(1,length(theseepis));
    filestonormalise_list = [];
    for i = 1:length(theseepis)
        filestonormalise{i} = spm_select('ExtFPList',outpath,['^topup_' blocksin{crun}{theseepis(i)}],1:minvols(crun));
        filestonormalise_list = [filestonormalise_list; filestonormalise{i}];
    end
    inputs{2, crun} = cellstr(filestonormalise_list);
    % % First is for SPM segment, second for CAT12
    %inputs{3, crun} = cellstr([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '/y_' blocksin{crun}{find(strcmp(blocksout{crun},'structural'))}]);
    inputs{3, crun} = cellstr([outpath 'mri/y_structural_csf.nii']);
    inputs{4, crun} = cellstr([outpath 'structural.nii,1']);
    % % First is for SPM segment, second for CAT12
    %inputs{5, crun} = cellstr([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '/y_' blocksin{crun}{find(strcmp(blocksout{crun},'structural'))}]);
    inputs{5, crun} = cellstr([outpath 'mri/y_structural_csf.nii']);
    inputs{6, crun} = cellstr([outpath 'structural.nii,1']);
end

normalisesmoothworkedcorrectly = zeros(1,nrun);
jobs = repmat(jobfile, 1, 1);

parfor crun = 1:nrun
    spm('defaults', 'fMRI');
    spm_jobman('initcfg')
    try
        spm_jobman('run', jobs, inputs{:,crun});
        normalisesmoothworkedcorrectly(crun) = 1;
    catch
        normalisesmoothworkedcorrectly(crun) = 0;
    end
end

if ~all(normalisesmoothworkedcorrectly)
    error('failed at normalise and smooth');
end

%% Now do a univariate SPM analysis (currently only implemented for 3 or 4 runs)
for this_smooth = [3,8];
    nrun = size(subjects,2); % enter the number of runs here
    jobfile = {};
    jobfile{3} = {[scriptdir 'module_univariate_3runs_noneutral_job.m']};
    jobfile{4} = {[scriptdir 'module_univariate_4runs_noneutral_job.m']};
    inputs = cell(0, nrun);
    
    for crun = 1:nrun
        theseepis = find(strncmp(blocksout{crun},'Run',3));
        outpath = [preprocessedpathstem subjects{crun} '/'];
        filestoanalyse = cell(1,length(theseepis));
        
        tempDesign = module_get_event_times_AFC4(subjects{crun},dates{crun},length(theseepis),minvols(crun));
        
        inputs{1, crun} = cellstr([outpath 'stats4_' num2str(this_smooth)]);
        for sess = 1:length(theseepis)
            filestoanalyse{sess} = spm_select('ExtFPList',outpath,['^s' num2str(this_smooth) 'wtopup_' blocksin{crun}{theseepis(sess)}],1:minvols(crun));
            inputs{(8*(sess-1))+2, crun} = cellstr(filestoanalyse{sess});
            inputs{(8*(sess-1))+3, crun} = cat(2, tempDesign{sess}{1:16})';
            inputs{(8*(sess-1))+4, crun} = cat(2, tempDesign{sess}{17:32})';
            inputs{(8*(sess-1))+5, crun} = cat(2, tempDesign{sess}{33:48})';
            inputs{(8*(sess-1))+6, crun} = cat(2, tempDesign{sess}{49:64})';
            %         inputs{(8*(sess-1))+7, crun} = cat(2, tempDesign{sess}{65:80})';
            %         inputs{(8*(sess-1))+8, crun} = cat(2, tempDesign{sess}{81:96})';
            inputs{(8*(sess-1))+7, crun} = cat(2, tempDesign{sess}{[97:112, 129]})';
            inputs{(8*(sess-1))+8, crun} = cat(2, tempDesign{sess}{[113:128, 130]})';
            inputs{(8*(sess-1))+9, crun} = cellstr([outpath 'rp_topup_' blocksin{crun}{theseepis(sess)}(1:end-4) '.txt']);
        end
        jobs{crun} = jobfile{length(theseepis)};
        
    end
    
    SPMworkedcorrectly = zeros(1,nrun);
    parfor crun = 1:nrun
        spm('defaults', 'fMRI');
        spm_jobman('initcfg')
        try
            spm_jobman('run', jobs{crun}, inputs{:,crun});
            SPMworkedcorrectly(crun) = 1;
        catch
            SPMworkedcorrectly(crun) = 0;
        end
    end
    
    if ~all(SPMworkedcorrectly)
        error(['failed at SPM ' num2str(this_smooth) 'mm']);
    end
end

%% Now create a univariate second level SPM with Age as a covariate - one per condition of interest
for this_smooth = [3,8];
    exclude_bad = 0; % I have now excluded only the relevant EPI sequences above.
    bad_scans = {
        'P7P16' % Left frontal hole
        'P7C18' % Bilateral frontal holes
        };
    
    age_lookup = readtable('Pinfa_ages.csv');
    all_conditions = {
        'con_0005.nii','Match > Mismatch';
        'con_0010.nii','Mismatch > Match';
        'con_0015.nii','Normal > Written';
        'con_0020.nii','Written > Normal';
        'con_0025.nii','Normal > Silence';
        'con_0030.nii','Clear > Unclear';
        'con_0035.nii','Unclear > Clear';
        'con_0035.nii','Clarity Congruency Interaction'};
    expected_sessions = 4;
    
    visual_check = 0;
    nrun = size(all_conditions,1); % enter the number of runs here
    %jobfile = {'/group/language/data/thomascope/vespa/SPM12version/Standalone preprocessing pipeline/tc_source/batch_forwardmodel_job_noheadpoints.m'};
    
    this_scan = {};
    this_t_scan = {};
    firstlevel_folder = ['stats4_' num2str(this_smooth)];
    
    jobfile = {'/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/module_secondlevel_job.m'};
    jobs = repmat(jobfile, 1, nrun);
    inputs = cell(4, nrun);
    
    for this_condition = 1:nrun
        group1_mrilist = {}; %NB: Patient MRIs, so here group 2 (sorry)
        group1_ages = [];
        group2_mrilist = {};
        group2_ages = [];
        
        if exclude_bad
            inputs{1, this_condition} = cellstr([preprocessedpathstem firstlevel_folder '_nobad' filesep all_conditions{this_condition,2}]);
        else
            inputs{1, this_condition} = cellstr([preprocessedpathstem firstlevel_folder filesep all_conditions{this_condition,2}]);
        end
        for crun = 1:size(subjects,2)
            if exclude_bad
                if any(strcmp(bad_scans,subjects{crun}))
                    continue
                end
            end
            this_age = age_lookup.Age(strcmp(age_lookup.Study_ID,subjects{crun}));
            this_spm_temp = load([preprocessedpathstem subjects{crun} filesep firstlevel_folder filesep 'SPM.mat']);
            if length(this_spm_temp.SPM.Sess)~=expected_sessions
                disp([subjects{crun} ' has ' num2str(length(this_spm_temp.SPM.Sess)) ' sessions when ' num2str(expected_sessions) ' expected. Check this is what you want'])
                con_num = str2num(all_conditions{this_condition,1}(5:8));
                new_con_num = (con_num/(expected_sessions+1))*(length(this_spm_temp.SPM.Sess)+1);
                disp(['Replacing contrast ' num2str(con_num,'%04.f') ' with ' num2str(new_con_num,'%04.f')])
                new_contrast_name = strrep(all_conditions{this_condition,1},num2str(con_num,'%04.f'),num2str(new_con_num,'%04.f'));
                this_scan(crun) = cellstr([preprocessedpathstem subjects{crun} filesep firstlevel_folder filesep new_contrast_name]);
                this_t_scan(crun) = cellstr([preprocessedpathstem subjects{crun} filesep firstlevel_folder filesep strrep(new_contrast_name,'con','spmT')]);
            else
            this_scan(crun) = cellstr([preprocessedpathstem subjects{crun} filesep firstlevel_folder filesep all_conditions{this_condition,1}]);
            this_t_scan(crun) = cellstr([preprocessedpathstem subjects{crun} filesep firstlevel_folder filesep strrep(all_conditions{this_condition,1},'con','spmT')]);
            end
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
        if visual_check
            spm_check_registration(this_t_scan{~cellfun(@isempty,this_t_scan)}) % Optional visual check of your input images (don't need to be aligned or anything, just to see they're all structurals and exist)
            input('Press any key to proceed to second level with these scans')
        end
    end
    
    secondlevelworkedcorrectly = zeros(1,nrun);
    parfor crun = 1:nrun
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



% %% Now create a more complex SPM for future multivariate analysis (currently only implemented for 3 or 4 runs) - Native space (s3r)
% nrun = size(subjects,2); % enter the number of runs here
% jobfile = {};
% jobfile{3} = {[scriptdir 'module_univariate_3runs_complex_job.m']};
% jobfile{4} = {[scriptdir 'module_univariate_4runs_complex_job.m']};
% inputs = cell(0, nrun);
% 
% for crun = 1:nrun
%     theseepis = find(strncmp(blocksout{crun},'Run',3));
%     outpath = [preprocessedpathstem subjects{crun} '/'];
%     filestoanalyse = cell(1,length(theseepis));
%     
%     tempDesign = module_get_complex_event_times_AFC4(subjects{crun},dates{crun},length(theseepis),minvols(crun));
%     
%     inputs{1, crun} = cellstr([outpath 'stats4_multi_3']);
%     for sess = 1:length(theseepis)
%         filestoanalyse{sess} = spm_select('ExtFPList',outpath,['^s3rtopup_' blocksin{crun}{theseepis(sess)}],1:minvols(crun)); %Native space image is s3r, standard space is s3w
%         inputs{(100*(sess-1))+2, crun} = cellstr(filestoanalyse{sess});
%         for cond_num = 1:80
%             inputs{(100*(sess-1))+2+cond_num, crun} = cat(2, tempDesign{sess}{cond_num})';
%         end
%         for cond_num = 81:96 %Response trials
%             inputs{(100*(sess-1))+2+cond_num, crun} = cat(2, tempDesign{sess}{cond_num+32})';
%         end
%         for cond_num = 97 %Button press
%             inputs{(100*(sess-1))+2+cond_num, crun} = cat(2, tempDesign{sess}{81})';
%         end
%         for cond_num = 98 %Absent sound (written only)
%             inputs{(100*(sess-1))+2+cond_num, crun} = cat(2, tempDesign{sess}{129})';
%         end
%         inputs{(100*(sess-1))+101, crun} = cellstr([outpath 'rp_topup_' blocksin{crun}{theseepis(sess)}(1:end-4) '.txt']);
%     end
%     jobs{crun} = jobfile{length(theseepis)};
%     if any(cellfun(@isempty,inputs(:,crun))) % In case of trunkated run where an event did not occur, put it at the very end of the run so it isn't modelled but SPM doesn't crash
%         inputs{find(cellfun(@isempty,inputs(:,crun))),crun} = tr*(length(filestoanalyse{sess})-1);
%     end
% end
% 
% SPMworkedcorrectly = zeros(1,nrun);
% parfor crun = 1:nrun
%     spm('defaults', 'fMRI');
%     spm_jobman('initcfg')
%     if exist([inputs{1,crun}{1} '/SPM.mat']) && ~SPMworkedcorrectly(crun)
%         delete([inputs{1,crun}{1} '/SPM.mat'])
%     elseif exist([inputs{1,crun}{1} '/SPM.mat']) && SPMworkedcorrectly(crun)
%         continue
%     end
%     
%     try
%         spm_jobman('run', jobs{crun}, inputs{:,crun});
%         SPMworkedcorrectly(crun) = 1;
%     catch
%         SPMworkedcorrectly(crun) = 0;
%     end
% end
% 
% %Now repeat with 8mm smoothing
% 
% for crun = 1:nrun
%     theseepis = find(strncmp(blocksout{crun},'Run',3));
%     outpath = [preprocessedpathstem subjects{crun} '/'];
%     filestoanalyse = cell(1,length(theseepis));
%     
%     tempDesign = module_get_complex_event_times_AFC4(subjects{crun},dates{crun},length(theseepis),minvols(crun));
%     
%     inputs{1, crun} = cellstr([outpath 'stats4_multi_8']);
%     for sess = 1:length(theseepis)
%         filestoanalyse{sess} = spm_select('ExtFPList',outpath,['^s8rtopup_' blocksin{crun}{theseepis(sess)}],1:minvols(crun));
%         inputs{(100*(sess-1))+2, crun} = cellstr(filestoanalyse{sess});
%         for cond_num = 1:80
%             inputs{(100*(sess-1))+2+cond_num, crun} = cat(2, tempDesign{sess}{cond_num})';
%         end
%         for cond_num = 81:96 %Response trials
%             inputs{(100*(sess-1))+2+cond_num, crun} = cat(2, tempDesign{sess}{cond_num+32})';
%         end
%         for cond_num = 97 %Button press
%             inputs{(100*(sess-1))+2+cond_num, crun} = cat(2, tempDesign{sess}{81})';
%         end
%         for cond_num = 98 %Absent sound (written only)
%             inputs{(100*(sess-1))+2+cond_num, crun} = cat(2, tempDesign{sess}{129})';
%         end
%         inputs{(100*(sess-1))+101, crun} = cellstr([outpath 'rp_topup_' blocksin{crun}{theseepis(sess)}(1:end-4) '.txt']);
%     end
%     jobs{crun} = jobfile{length(theseepis)};
%     if any(cellfun(@isempty,inputs(:,crun))) % In case of trunkated run where an event did not occur, put it at the very end of the run so it isn't modelled but SPM doesn't crash
%         inputs{find(cellfun(@isempty,inputs(:,crun))),crun} = tr*(length(filestoanalyse{sess})-1);
%     end
% end
% 
% SPMworkedcorrectly = zeros(1,nrun);
% parfor crun = 1:nrun
%     spm('defaults', 'fMRI');
%     spm_jobman('initcfg')
%     try
%         spm_jobman('run', jobs{crun}, inputs{:,crun});
%         SPMworkedcorrectly(crun) = 1;
%     catch
%         SPMworkedcorrectly(crun) = 0;
%     end
% end

% %% Now create a more complex SPM with variable levels of AR whitening, with word omissions specified
%
% jobfile = {};
% jobfile{3} = {[scriptdir 'module_univariate_3runs_complex_AR_job.m']};
% jobfile{4} = {[scriptdir 'module_univariate_4runs_complex_AR_job.m']};
%
% all_aros = [1 3 6 12]; %Autoregressive model order
% nrun = size(subjects,2)*length(all_aros); % enter the number of runs here
% inputs = cell(0, size(subjects,2),length(all_aros));
% for this_aro = 1:length(all_aros);
% for crun = 1:size(subjects,2)
%     aro = all_aros(this_aro);
%
%     theseepis = find(strncmp(blocksout{crun},'Run',3));
%     outpath = [preprocessedpathstem subjects{crun} '/'];
%     filestoanalyse = cell(1,length(theseepis));
%
%     tempDesign = module_get_complex_event_times(subjects{crun},dates{crun},length(theseepis),minvols(crun));
%
%     inputs{1, crun, this_aro} = cellstr([outpath 'stats4_multi_AR' num2str(aro)]);
%     for sess = 1:length(theseepis)
%         filestoanalyse{sess} = spm_select('ExtFPList',outpath,['^s3topup_' blocksin{crun}{theseepis(sess)}],1:minvols(crun));
%         inputs{(100*(sess-1))+2, crun, this_aro} = cellstr(filestoanalyse{sess});
%         for cond_num = 1:80
%             inputs{(100*(sess-1))+2+cond_num, crun, this_aro} = cat(2, tempDesign{sess}{cond_num})';
%         end
%         for cond_num = 81:96 %Response trials
%             inputs{(100*(sess-1))+2+cond_num, crun, this_aro} = cat(2, tempDesign{sess}{cond_num+32})';
%         end
%         for cond_num = 97 %Button press
%             inputs{(100*(sess-1))+2+cond_num, crun, this_aro} = cat(2, tempDesign{sess}{81})';
%         end
%         for cond_num = 98 %Absent sound (written only)
%             inputs{(100*(sess-1))+2+cond_num, crun, this_aro} = cat(2, tempDesign{sess}{129})';
%         end
%         inputs{(100*(sess-1))+101, crun, this_aro} = cellstr([outpath 'rp_topup_' blocksin{crun}{theseepis(sess)}(1:end-4) '.txt']);
%     end
%     %inputs{(100*(sess-1))+102, crun, this_aro} = 'AR(1)';
%     inputs{(100*(sess-1))+102, crun, this_aro} = aro;
%     jobs{crun} = jobfile{length(theseepis)};
%
% end
% end
%
%
% try
%     matlabpool 'close'
% catch
%     delete(gcp)
% end
%
%
% workersrequested = 24;
% workerpool = cbupool(workersrequested);
% workerpool.ResourceTemplate=['-l nodes=^N^,mem=768GB,walltime=168:00:00'];
% try
%     matlabpool(workerpool)
% catch
%     parpool(workerpool,workerpool.NumWorkers)
% end
%
% all_combs = combvec(1:size(subjects,2),1:length(all_aros))';
% SPMworkedcorrectly = zeros(1,size(all_combs,1));
% for thisone = 1:size(all_combs,1)
%     crun = all_combs(thisone,1);
%     this_aro = all_combs(thisone,2);
%     spm('defaults', 'fMRI');
%     spm_jobman('initcfg')
%     try
%         spm_jobman('run', jobs{crun}, inputs{:,crun, this_aro});
%         SPMworkedcorrectly(thisone) = 1;
%     catch
%         SPMworkedcorrectly(thisone) = 0;
%     end
% end
 
%% Now create a ReML SPM without modelling the written word separately for future multivariate analysis (currently only implemented for 3 or 4 runs) - Native space (s3r)
nrun = size(subjects,2); % enter the number of runs here
jobfile = {};
jobfile{3} = {[scriptdir 'module_univariate_3runs_noabsent_job.m']};
jobfile{4} = {[scriptdir 'module_univariate_4runs_noabsent_job.m']};
inputs = cell(0, nrun);

for crun = 1:nrun
    theseepis = find(strncmp(blocksout{crun},'Run',3));
    outpath = [preprocessedpathstem subjects{crun} '/'];
    filestoanalyse = cell(1,length(theseepis));
    
    tempDesign = module_get_complex_event_times_nowritten_AFC4(subjects{crun},dates{crun},length(theseepis),minvols(crun));
    
    inputs{1, crun} = cellstr([outpath 'stats4_multi_3_nowritten2']);
    for sess = 1:length(theseepis)
        filestoanalyse{sess} = spm_select('ExtFPList',outpath,['^s3rtopup_' blocksin{crun}{theseepis(sess)}],1:minvols(crun));
        inputs{(99*(sess-1))+2, crun} = cellstr(filestoanalyse{sess});
        for cond_num = 1:80
            inputs{(99*(sess-1))+2+cond_num, crun} = cat(2, tempDesign{sess}{cond_num})';
        end
        for cond_num = 81:96 %Response trials
            inputs{(99*(sess-1))+2+cond_num, crun} = cat(2, tempDesign{sess}{cond_num+32})';
        end
        for cond_num = 97 %Button press
            inputs{(99*(sess-1))+2+cond_num, crun} = cat(2, tempDesign{sess}{81})';
        end
        inputs{(99*(sess-1))+100, crun} = cellstr([outpath 'rp_topup_' blocksin{crun}{theseepis(sess)}(1:end-4) '.txt']);
        if any(cellfun(@isempty,inputs(1:(99*(sess-1))+100,crun))) % In case of trunkated run where an event did not occur, put it at the very end of the run so it isn't modelled but SPM doesn't crash
            inputs{find(cellfun(@isempty,inputs(1:(99*(sess-1))+100,crun))),crun} = tr*(length(filestoanalyse{sess})-1);
        end
    end
    jobs{crun} = jobfile{length(theseepis)};
end

SPMworkedcorrectly = zeros(1,nrun);
parfor crun = 1:nrun
    spm('defaults', 'fMRI');
    spm_jobman('initcfg')
    try
        spm_jobman('run', jobs{crun}, inputs{:,crun});
        SPMworkedcorrectly(crun) = 1;
    catch
        SPMworkedcorrectly(crun) = 0;
    end
end

nrun = size(subjects,2); % enter the number of runs here
jobfile = {};
jobfile{3} = {[scriptdir 'module_univariate_3runs_noabsent_job.m']};
jobfile{4} = {[scriptdir 'module_univariate_4runs_noabsent_job.m']};
inputs = cell(0, nrun);

for crun = 1:nrun
    theseepis = find(strncmp(blocksout{crun},'Run',3));
    outpath = [preprocessedpathstem subjects{crun} '/'];
    filestoanalyse = cell(1,length(theseepis));
    
    tempDesign = module_get_complex_event_times_nowritten_AFC4(subjects{crun},dates{crun},length(theseepis),minvols(crun));
    
    inputs{1, crun} = cellstr([outpath 'stats4_multi_8_nowritten2']);
    for sess = 1:length(theseepis)
        filestoanalyse{sess} = spm_select('ExtFPList',outpath,['^s8rtopup_' blocksin{crun}{theseepis(sess)}],1:minvols(crun));
        inputs{(99*(sess-1))+2, crun} = cellstr(filestoanalyse{sess});
        for cond_num = 1:80
            inputs{(99*(sess-1))+2+cond_num, crun} = cat(2, tempDesign{sess}{cond_num})';
        end
        for cond_num = 81:96 %Response trials
            inputs{(99*(sess-1))+2+cond_num, crun} = cat(2, tempDesign{sess}{cond_num+32})';
        end
        for cond_num = 97 %Button press
            inputs{(99*(sess-1))+2+cond_num, crun} = cat(2, tempDesign{sess}{81})';
        end
        inputs{(99*(sess-1))+100, crun} = cellstr([outpath 'rp_topup_' blocksin{crun}{theseepis(sess)}(1:end-4) '.txt']);
        if any(cellfun(@isempty,inputs(1:(99*(sess-1))+100,crun))) % In case of trunkated run where an event did not occur, put it at the very end of the run so it isn't modelled but SPM doesn't crash
            inputs{find(cellfun(@isempty,inputs(1:(99*(sess-1))+100,crun))),crun} = tr*(length(filestoanalyse{sess})-1);
        end
    end
    jobs{crun} = jobfile{length(theseepis)};
end

SPMworkedcorrectly = zeros(1,nrun);
parfor crun = 1:nrun
    spm('defaults', 'fMRI');
    spm_jobman('initcfg')
    try
        spm_jobman('run', jobs{crun}, inputs{:,crun});
        SPMworkedcorrectly(crun) = 1;
    catch
        SPMworkedcorrectly(crun) = 0;
    end
end

%% Now run the cross validated Mahalanobis distance and RSM on each subject 
nrun = size(subjects,2); % enter the number of runs here

parfor crun = 1:nrun
    addpath(genpath('./RSA_scripts'))
    GLMDir = [preprocessedpathstem subjects{crun} '/stats4_multi_3_nowritten2'];
    TDTCrossnobisAnalysis_1Subj(GLMDir)
end



%% Now create univariate masks for later MVPA

t_thresh = 3.11; % p<0.001 uncorrected
smoothing_kernels = [3, 8];

for smoo = smoothing_kernels
    for crun = 1:nrun
        spmpath = [preprocessedpathstem subjects{crun} '/stats4_multi_' num2str(smoo) '/'];
        outpath = [preprocessedpathstem subjects{crun} '/'];
        thisSPM = load([spmpath 'SPM.mat']);
        writtenindex = structfind(thisSPM.SPM.xCon,'name','Normal<Written - All Sessions');
        if numel(writtenindex) ~= 1
            error('Something went wrong with finding the written mask condition')
        end
        soundindex = structfind(thisSPM.SPM.xCon,'name','Normal>Written - All Sessions');
        if numel(soundindex) ~= 1
            error('Something went wrong with finding the sound mask condition')
        end
        spm_imcalc([spmpath 'spmT_' sprintf('%04d',soundindex) '.nii'],[outpath 'mask_' num2str(smoo) '_sound_001.nii'],'i1>3.11')
        spm_imcalc([spmpath 'spmT_' sprintf('%04d',writtenindex) '.nii'],[outpath 'mask_' num2str(smoo) '_written_001.nii'],'i1>3.11')
    end
end

%% Now create anatomical masks for later MVPA
nrun = size(subjects,2); % enter the number of runs here

% First create masks
search_labels = {
    'Left STG'
    'Left PT'
    'Left PrG'
    'Left FO'
    'Left TrIFG'
    };

for crun = 1:nrun
    outpath = [preprocessedpathstem subjects{crun} '/'];
    currentdir = pwd;
    cd(outpath)
    xA=spm_atlas('load','Neuromorphometrics');
    for i = 1:size(xA.labels,2)
        all_labels{i} = xA.labels(i).name;
    end
    
    S = cell(1,length(search_labels));
    for i = 1:length(search_labels)
        S{i} = find(strncmp(all_labels,search_labels{i},size(search_labels{i},2)));
    end
    
    for i = 1:size(S,2)
        fname=strcat(strrep(search_labels{i}, ' ', '_'),'.nii');
        VM=spm_atlas('mask',xA,xA.labels(S{i}).name);
        VM.fname=fname;
        spm_write_vol(VM,spm_read_vols(VM));
    end
    
    fname='atlas_all.nii';
    VM=spm_atlas('mask',xA,all_labels);
    VM.fname=fname;
    spm_write_vol(VM,spm_read_vols(VM));
    
    cd(currentdir)
    
end

maskcoregisterworkedcorrectly = zeros(1,nrun);

parfor crun = 1:nrun
    job = struct
    job.eoptions.cost_fun = 'nmi'
    job.eoptions.tol = [repmat(0.02,1,3), repmat(0.01,1,6), repmat(0.001,1,3)];
    job.eoptions.sep = [4 2];
    job.eoptions.fwhm = [7 7];
    
    outpath = [preprocessedpathstem subjects{crun} '/'];
    job.ref = {[outpath 'wstructural.nii']};
    theseepis = find(strncmp(blocksout{crun},'Run',3));
    job.source = {[outpath 'atlas_all.nii']};
    
    filestocoregister = cell(1,length(theseepis));
    filestocoregister_list = [];
    for i = 1:length(search_labels)
        filestocoregister{i} = strcat(outpath,strrep(search_labels{i}, ' ', '_'),'.nii');
        filestocoregister_list = strvcat(filestocoregister_list, filestocoregister{i});
    end
    filestocoregister = cellstr(filestocoregister_list);
    
    job.other = filestocoregister
    
    try
        spm_run_coreg(job)
        
        P = char(job.ref{:},job.source{:},job.other{:});
        %inflate the ROIs a bit to account for smaller brains than template
        for thisone=3:size(P,1)
            dilate_image_spm(P(thisone,:),5)
            %             spm_imcalc(P(thisone,:), P(thisone,:), 'i1*10');
            %             spm_smooth(P(thisone,:),P(thisone,:),10);
            %             spm_imcalc(P(thisone,:),P(thisone,:),'i1>1');
        end
        flags=struct;
        flags.interp = 0;
        spm_reslice(P,flags)
        
        
        maskcoregisterworkedcorrectly(crun) = 1;
    catch
        maskcoregisterworkedcorrectly(crun) = 0;
    end
end

%% Now create anatomical masks for later MVPA with fat Neuromorphometrics
nrun = size(subjects,2); % enter the number of runs here

% First create masks
search_labels = {
    'Left Superior Temporal Gyrus'
    'Left Angular Gyrus'
    'Left Precentral Gyrus'
    'Left Frontal Operculum'
    'Left Inferior Frontal Angular Gyrus'
    'Right Superior Temporal Gyrus'
    'Right Angular Gyrus'
    'Right Precentral Gyrus'
    'Right Frontal Operculum'
    'Right Inferior Frontal Angular Gyrus'
    'Left Cerebellar Lobule Cerebellar Vermal Lobules VI-VII'
    'Right Cerebellar Lobule Cerebellar Vermal Lobules VI-VII'
    };


cat_install_atlases

for crun = 1:nrun
    outpath = [preprocessedpathstem subjects{crun} '/'];
    currentdir = pwd;
    cd(outpath)
    xA=spm_atlas('load','dartel_neuromorphometrics');
    for i = 1:size(xA.labels,2)
        all_labels{i} = xA.labels(i).name;
    end
    
    S = cell(1,length(search_labels));
    for i = 1:length(search_labels)
        S{i} = find(strcmp(all_labels,search_labels{i}));
    end
    
    for i = 1:size(S,2)
        fname=strcat(strrep(search_labels{i}, ' ', '_'),'.nii');
        VM=spm_atlas('mask',xA,xA.labels(S{i}).name);
        VM.fname=fname;
        spm_write_vol(VM,spm_read_vols(VM));
    end
    
    fname='atlas_all.nii';
    VM=spm_atlas('mask',xA,all_labels);
    VM.fname=fname;
    spm_write_vol(VM,spm_read_vols(VM));
    
    cd(currentdir)
    
end

maskcoregisterworkedcorrectly = zeros(1,nrun);

parfor crun = 1:nrun
    job = struct
    job.eoptions.cost_fun = 'nmi'
    job.eoptions.tol = [repmat(0.02,1,3), repmat(0.01,1,6), repmat(0.001,1,3)];
    job.eoptions.sep = [4 2];
    job.eoptions.fwhm = [7 7];
    
    outpath = [preprocessedpathstem subjects{crun} '/'];
    job.ref = {[outpath 'wstructural.nii']};
    theseepis = find(strncmp(blocksout{crun},'Run',3));
    job.source = {[outpath 'atlas_all.nii']};
    
    filestocoregister = cell(1,length(theseepis));
    filestocoregister_list = [];
    for i = 1:length(search_labels)
        filestocoregister{i} = strcat(outpath,strrep(search_labels{i}, ' ', '_'),'.nii');
        filestocoregister_list = strvcat(filestocoregister_list, filestocoregister{i});
    end
    filestocoregister = cellstr(filestocoregister_list);
    
    job.other = filestocoregister
    
    try
        %spm_run_coreg(job)
        
        P = char(job.ref{:},job.source{:},job.other{:});
        %         %inflate the ROIs a bit to account for smaller brains than template
        %         for thisone=3:size(P,1)
        %             dilate_image_spm(P(thisone,:),5)
        % %             spm_imcalc(P(thisone,:), P(thisone,:), 'i1*10');
        % %             spm_smooth(P(thisone,:),P(thisone,:),10);
        % %             spm_imcalc(P(thisone,:),P(thisone,:),'i1>1');
        %         end
        flags=struct;
        flags.interp = 0;
        spm_reslice(P,flags)
        
        
        maskcoregisterworkedcorrectly(crun) = 1;
    catch
        maskcoregisterworkedcorrectly(crun) = 0;
    end
end

%% Now begin the MVPA proper! RSA within the mask first
nrun = size(subjects,2); % enter the number of runs here
data_smoo = 3; %Smoothing on MVPA data
mask_smoo = 3; %Smoothing on MVPA mask
mask_cond = {'sound' 'written'};
% 16M4 16M12 16MM4 16MM12 16WO 16R BP Null 6Mov

avgRDM = cell(size(subjects,2),length(mask_cond),length(conditions));
stats_p_r = cell(size(subjects,2),length(mask_cond),length(conditions));

for crun = 1:nrun
    
    data_path = [preprocessedpathstem subjects{crun} '/stats4_multi_' num2str(data_smoo) '/'];
    
    for mask_cond_num = 1:length(mask_cond)
        mask_path = [preprocessedpathstem subjects{crun} '/mask_' num2str(mask_smoo) '_' mask_cond{mask_cond_num} '_001.nii'];
        for cond_num = 1:length(conditions)
            tpattern_numbers = 9+[1:16]+(16*(cond_num-1));
            [avgRDM{crun,mask_cond_num,cond_num}, stats_p_r{crun,mask_cond_num,cond_num}] = module_rsa_job(tpattern_numbers,mask_path,data_path,cond_num,conditions{cond_num});
        end
    end
end
all_avgRDM{1} = avgRDM;
all_stats{1} = stats_p_r;

data_smoo = 3; %Smoothing on MVPA data
mask_smoo = 8; %Smoothing on MVPA mask
mask_cond = {'sound' 'written'};
% 16M4 16M12 16MM4 16MM12 16WO 16R BP Null 6Mov

avgRDM = cell(size(subjects,2),length(mask_cond),length(conditions));
stats_p_r = cell(size(subjects,2),length(mask_cond),length(conditions));

for crun = 1:nrun
    
    data_path = [preprocessedpathstem subjects{crun} '/stats4_multi_' num2str(data_smoo) '/'];
    
    for mask_cond_num = 1:length(mask_cond)
        mask_path = [preprocessedpathstem subjects{crun} '/mask_' num2str(mask_smoo) '_' mask_cond{mask_cond_num} '_001.nii'];
        for cond_num = 1:length(conditions)
            tpattern_numbers = 9+[1:16]+(16*(cond_num-1));
            [avgRDM{crun,mask_cond_num,cond_num}, stats_p_r{crun,mask_cond_num,cond_num}] = module_rsa_job(tpattern_numbers,mask_path,data_path,cond_num,conditions{cond_num});
        end
    end
end
all_avgRDM{2} = avgRDM;
all_stats{2} = stats_p_r;


data_smoo = 8; %Smoothing on MVPA data
mask_smoo = 8; %Smoothing on MVPA mask
mask_cond = {'sound' 'written'};
% 16M4 16M12 16MM4 16MM12 16WO 16R BP Null 6Mov

avgRDM = cell(size(subjects,2),length(mask_cond),length(conditions));
stats_p_r = cell(size(subjects,2),length(mask_cond),length(conditions));

for crun = 1:nrun
    
    data_path = [preprocessedpathstem subjects{crun} '/stats4_multi_' num2str(data_smoo) '/'];
    
    for mask_cond_num = 1:length(mask_cond)
        mask_path = [preprocessedpathstem subjects{crun} '/mask_' num2str(mask_smoo) '_' mask_cond{mask_cond_num} '_001.nii'];
        for cond_num = 1:length(conditions)
            tpattern_numbers = 9+[1:16]+(16*(cond_num-1));
            [avgRDM{crun,mask_cond_num,cond_num}, stats_p_r{crun,mask_cond_num,cond_num}] = module_rsa_job(tpattern_numbers,mask_path,data_path,cond_num,conditions{cond_num});
        end
    end
end
all_avgRDM{3} = avgRDM;
all_stats{3} = stats_p_r;

data_smoo = 3; %Smoothing on MVPA data
%mask_cond = {'rLeft_STG.nii' 'rLeft_PrG.nii' 'rLeft_FO.nii'};
mask_cond = {'rLeft_Superior_Temporal_Gyrus.nii'
    'rLeft_Angular_Gyrus.nii'
    'rLeft_Precentral_Gyrus.nii'
    'rLeft_Frontal_Operculum.nii'
    'rLeft_Inferior_Frontal_Angular_Gyrus.nii'
    };
% 16M4 16M12 16MM4 16MM12 16WO 16R BP Null 6Mov

avgRDM = cell(size(subjects,2),length(mask_cond),length(conditions));
stats_p_r = cell(size(subjects,2),length(mask_cond),length(conditions));

for crun = 1:nrun
    
    data_path = [preprocessedpathstem subjects{crun} '/stats4_multi_' num2str(data_smoo) '/'];
    
    for mask_cond_num = 1:length(mask_cond)
        mask_path = [preprocessedpathstem subjects{crun} '/' mask_cond{mask_cond_num}];
        for cond_num = 1:length(conditions)
            tpattern_numbers = 9+[1:16]+(16*(cond_num-1));
            [avgRDM{crun,mask_cond_num,cond_num}, stats_p_r{crun,mask_cond_num,cond_num}] = module_rsa_job(tpattern_numbers,mask_path,data_path,cond_num,conditions{cond_num});
        end
    end
end
all_avgRDM{4} = avgRDM;
all_stats{4} = stats_p_r;

%% Try again with parallelisation of different AR model orders
addpath(genpath('/imaging/mlr/users/tc02/toolboxes')); %Where is the RSA toolbox?

%data_smoo = 3; %Smoothing on MVPA data
all_aros = [1 3 6 12];
mask_cond = {'rLeft_Superior_Temporal_Gyrus.nii'
    'rLeft_Angular_Gyrus.nii'
    'rLeft_Precentral_Gyrus.nii'
    'rLeft_Frontal_Operculum.nii'
    'rLeft_Inferior_Frontal_Angular_Gyrus.nii'
    'rRight_Superior_Temporal_Gyrus.nii'
    'rRight_Angular_Gyrus.nii'
    'rRight_Precentral_Gyrus.nii'
    'rRight_Frontal_Operculum.nii'
    'rRight_Inferior_Frontal_Angular_Gyrus.nii'
    'rLeft_Cerebellar_Lobule_Cerebellar_Vermal_Lobules_VI-VII.nii'
    'rRight_Cerebellar_Lobule_Cerebellar_Vermal_Lobules_VI-VII.nii'
    };
mask_short_cond = {'lSTG'
    'lAG'
    'lPrG'
    'lFO'
    'lIFG'
    'rSTG'
    'rAG'
    'rPrG'
    'rFO'
    'rIFG'
    'lXXX'
    'rXXX'};
all_combs = combvec(1:size(subjects,2),1:length(mask_cond),1:length(conditions),1:length(all_aros))';
pat_aro_combs = combvec(1:size(subjects,2),1:length(all_aros))';

%type = 't-pat'; % Run based on the t-patterns
type = 'beta'; % Run based on the beta-patterns

switch type
    case 'beta'
        %First denan the beta images
        for thisone = 1:size(pat_aro_combs,1)
            crun = pat_aro_combs(thisone,1);
            aro = all_aros(pat_aro_combs(thisone,2));
            data_path = [preprocessedpathstem subjects{crun} '/stats4_multi_AR' num2str(aro) '/'];
            beta_files = dir([data_path '/Cbeta_0*']);
            
            parfor i = 1:size(beta_files,1)
                module_fslmaths_job([data_path 'Cbeta_' sprintf('%04d',i) '.nii'],'-nan',[data_path 'Cbeta_denan_' sprintf('%04d',i) '.nii']); %Account for the fact that spm_read_vols crashes with nan
            end
        end
end

parfor thisone = 1:size(all_combs,1)
    crun = all_combs(thisone,1);
    mask_cond_num = all_combs(thisone,2);
    cond_num = all_combs(thisone,3);
    aro = all_aros(all_combs(thisone,4));
    switch type
        case 't-pat'
            module_run_rsa_AR(crun,cond_num,mask_cond{mask_cond_num},conditions{cond_num},aro) % Run based on the t-patterns
        case 'beta'
            module_run_rsa_AR_beta(crun,cond_num,mask_cond{mask_cond_num},conditions{cond_num},aro) %Run based on the beta patterns
    end
    %module_run_rsa(crun,cond_num,mask_cond{mask_cond_num},conditions{cond_num},data_smoo)
    %module_run_rsa(crun,cond_num,mask_cond{mask_cond_num},['Subj_' num2str(crun) '_mask_' mask_cond{mask_cond_num} '_cond_' conditions{cond_num} '_smo_' num2str(data_smoo)],data_smoo)
    
end

for crun = 1:size(subjects,2)
    for mask_cond_num = 1:length(mask_cond)
        for cond_num = 1:length(conditions)
            for aro = 1:length(all_aros)
                mask_name = mask_cond{mask_cond_num};
                mask_short_name = mask_short_cond{mask_cond_num};
                switch type
                    case 't-pat'
                        thesedata = load(['./RSA_results/RSA_results_subj' num2str(crun) '_' conditions{cond_num} '_mask_' mask_name(1:end-4) '_AR' num2str(all_aros(aro))],'avgRDM','stats_p_r');
                    case 'beta'
                        thesedata = load(['./RSA_results/RSA_results_beta_subj' num2str(crun) '_' conditions{cond_num} '_mask_' mask_name(1:end-4) '_AR' num2str(all_aros(aro))],'avgRDM','stats_p_r');
                end
                %thesedata = load(['./RSA_results/RSA_results_subj' num2str(crun) '_' 'Subj_' num2str(crun) '_mask_' mask_cond{mask_cond_num} '_cond_' conditions{cond_num} '_smo_' num2str(data_smoo) '_mask_' mask_name(1:end-4) '_smooth_' num2str(data_smoo)],'avgRDM','stats_p_r');
                avgRDM{crun,mask_cond_num,cond_num,aro} = thesedata.avgRDM;
                this_cond_name = strrep(avgRDM{crun,mask_cond_num,cond_num,aro}.name,'Mismatch ','MM');
                this_cond_name = strrep(this_cond_name,'Match ','M');
                this_cond_name = strrep(this_cond_name,'RDM across sessions | condition ','');
                avgRDM{crun,mask_cond_num,cond_num,aro}.name = ['S' num2str(crun) this_cond_name '_' mask_short_name '_AR' num2str(all_aros(aro))];
                stats_p_r{crun,mask_cond_num,cond_num,aro} = thesedata.stats_p_r;
            end
        end
    end
end

userOptions.candRDMdifferencesTest='conditionRFXbootstrap';
userOptions.nBootstrap=100; % XXX CHange to 10000 when code finalised

judgmentRDM.RDM = zeros(16,16);
judgmentRDM.RDM(1:17:end) = 1;
judgmentRDM.RDM(2:68:end) = 1/3;
judgmentRDM.RDM(3:68:end) = 1/3;
judgmentRDM.RDM(4:68:end) = 1/3;
judgmentRDM.RDM(17:68:end) = 1/3;
judgmentRDM.RDM(19:68:end) = 1/3;
judgmentRDM.RDM(20:68:end) = 1/3;
judgmentRDM.RDM(33:68:end) = 1/3;
judgmentRDM.RDM(34:68:end) = 1/3;
judgmentRDM.RDM(36:68:end) = 1/3;
judgmentRDM.RDM(49:68:end) = 1/3;
judgmentRDM.RDM(50:68:end) = 1/3;
judgmentRDM.RDM(51:68:end) = 1/3;

judgmentRDM.RDM = 1-judgmentRDM.RDM;
judgmentRDM.name = 'vowels only';

base_figureindex = 250;
userOptions.figureIndex = [260, 360];

for crun = 1:size(subjects,2)
    for aro = 1:length(all_aros)
        userOptions.figureIndex = [base_figureindex+10*crun+aro, base_figureindex+200+10*crun+aro];
        
        subj_stats_p_r{crun,aro}=compareRefRDM2candRDMs(judgmentRDM, avgRDM(crun,:,:,aro), userOptions);
        
        
    end
end

%% Analyse by condition and brain region
addpath(genpath('/imaging/mlr/users/tc02/toolboxes')); %Where is the RSA toolbox?

if opennewanalysispool == 1
    %Re-open Parpool with larger worker pool
    
    cd('/group/language/data/thomascope/')
    try
        matlabpool 'close'
    catch
        delete(gcp)
    end
    workerpool = cbupool(24);
    workerpool.ResourceTemplate=['-l nodes=^N^,mem=192GB,walltime=168:00:00'];
    try
        matlabpool(workerpool)
    catch
        parpool(workerpool,workerpool.NumWorkers)
    end
    cd(currentdr)
end

all_smos = 3; %Smoothing on MVPA data

mask_cond = {'rLeft_Superior_Temporal_Gyrus.nii'
    'rLeft_Angular_Gyrus.nii'
    'rLeft_Precentral_Gyrus.nii'
    'rLeft_Frontal_Operculum.nii'
    'rLeft_Inferior_Frontal_Angular_Gyrus.nii'
    'rRight_Superior_Temporal_Gyrus.nii'
    'rRight_Angular_Gyrus.nii'
    'rRight_Precentral_Gyrus.nii'
    'rRight_Frontal_Operculum.nii'
    'rRight_Inferior_Frontal_Angular_Gyrus.nii'
    'rLeft_Cerebellar_Lobule_Cerebellar_Vermal_Lobules_VI-VII.nii'
    'rRight_Cerebellar_Lobule_Cerebellar_Vermal_Lobules_VI-VII.nii'
    };
mask_short_cond = {'lSTG'
    'lAG'
    'lPrG'
    'lFO'
    'lIFG'
    'rSTG'
    'rAG'
    'rPrG'
    'rFO'
    'rIFG'
    'lXXX'
    'rXXX'};
all_combs = combvec(1:size(subjects,2),1:length(mask_cond),1:length(conditions),1:length(all_smos))';
pat_smo_combs = combvec(1:size(subjects,2),1:length(all_smos))';

type = 't-pat'; % Run based on the t-patterns
%type = 'beta'; % Run based on the beta-patterns

switch type
    case 'beta'
        %First denan the beta images
        for thisone = 1:size(pat_smo_combs,1)
            crun = pat_smo_combs(thisone,1);
            smo = all_smos(pat_smo_combs(thisone,2));
            data_path = [preprocessedpathstem subjects{crun} '/stats_multi_' num2str(smo) '_nowritten/'];
            beta_files = dir([data_path '/Cbeta_0*']);
            
            parfor i = 1:size(beta_files,1)
                module_fslmaths_job([data_path 'Cbeta_' sprintf('%04d',i) '.nii'],'-nan',[data_path 'Cbeta_denan_' sprintf('%04d',i) '.nii']); %Account for the fact that spm_read_vols crashes with nan
            end
        end
end

parfor thisone = 1:size(all_combs,1)
    crun = all_combs(thisone,1);
    mask_cond_num = all_combs(thisone,2);
    cond_num = all_combs(thisone,3);
    smo = all_smos(all_combs(thisone,4));
    try
        switch type
            case 't-pat'
                module_run_rsa(crun,cond_num,mask_cond{mask_cond_num},conditions{cond_num},smo,setup_file) % Run based on the t-patterns
            case 'beta'
                module_run_rsa_beta(crun,cond_num,mask_cond{mask_cond_num},conditions{cond_num},smo,setup_file) %Run based on the beta patterns
        end
    end
    %module_run_rsa(crun,cond_num,mask_cond{mask_cond_num},conditions{cond_num},data_smoo)
    %module_run_rsa(crun,cond_num,mask_cond{mask_cond_num},['Subj_' num2str(crun) '_mask_' mask_cond{mask_cond_num} '_cond_' conditions{cond_num} '_smo_' num2str(data_smoo)],data_smoo)
    
end

for crun = 1:size(subjects,2)
    for mask_cond_num = 1:length(mask_cond)
        for cond_num = 1:length(conditions)
            for smo = 1:length(all_smos)
                try
                    mask_name = mask_cond{mask_cond_num};
                    mask_short_name = mask_short_cond{mask_cond_num};
                    switch type
                        case 't-pat'
                            thesedata = load(['./RSA_results/RSA_results_nowritten2_subj' num2str(crun) '_' conditions{cond_num} '_mask_' mask_name(1:end-4) '_smooth_' num2str(all_smos(smo))],'avgRDM','stats_p_r');
                        case 'beta'
                            thesedata = load(['./RSA_results/RSA_results_beta_nowritten2_subj' num2str(crun) '_' conditions{cond_num} '_mask_' mask_name(1:end-4) '_smooth_' num2str(all_smos(smo))],'avgRDM','stats_p_r');
                    end
                    %thesedata = load(['./RSA_results/RSA_results_subj' num2str(crun) '_' 'Subj_' num2str(crun) '_mask_' mask_cond{mask_cond_num} '_cond_' conditions{cond_num} '_smo_' num2str(data_smoo) '_mask_' mask_name(1:end-4) '_smooth_' num2str(data_smoo)],'avgRDM','stats_p_r');
                    avgRDM{crun,mask_cond_num,cond_num,smo} = thesedata.avgRDM;
                    this_cond_name = strrep(avgRDM{crun,mask_cond_num,cond_num,smo}.name,'Mismatch ','MM');
                    this_cond_name = strrep(this_cond_name,'Match ','M');
                    this_cond_name = strrep(this_cond_name,'RDM across sessions | condition ','');
                    avgRDM{crun,mask_cond_num,cond_num,smo}.name = ['S' num2str(crun) this_cond_name '_' mask_short_name '_sm' num2str(all_smos(smo))];
                    stats_p_r{crun,mask_cond_num,cond_num,smo} = thesedata.stats_p_r;
                end
            end
        end
    end
end

userOptions.candRDMdifferencesTest='conditionRFXbootstrap';
userOptions.nBootstrap=100; % XXX CHange to 10000 when code finalised

judgmentRDM.RDM = zeros(16,16);
judgmentRDM.RDM(1:17:end) = 1;
judgmentRDM.RDM(2:68:end) = 1/3;
judgmentRDM.RDM(3:68:end) = 1/3;
judgmentRDM.RDM(4:68:end) = 1/3;
judgmentRDM.RDM(17:68:end) = 1/3;
judgmentRDM.RDM(19:68:end) = 1/3;
judgmentRDM.RDM(20:68:end) = 1/3;
judgmentRDM.RDM(33:68:end) = 1/3;
judgmentRDM.RDM(34:68:end) = 1/3;
judgmentRDM.RDM(36:68:end) = 1/3;
judgmentRDM.RDM(49:68:end) = 1/3;
judgmentRDM.RDM(50:68:end) = 1/3;
judgmentRDM.RDM(51:68:end) = 1/3;

judgmentRDM.RDM = 1-judgmentRDM.RDM;
judgmentRDM.name = 'vowels only';

base_figureindex = 250;
userOptions.figureIndex = [260, 360];

for crun = 1:size(subjects,2)
    for smo = 1:length(all_smos)
        try
            userOptions.figureIndex = [base_figureindex+10*crun+smo, base_figureindex+200+10*crun+smo];
            
            subj_stats_p_r{crun,smo}=compareRefRDM2candRDMs(judgmentRDM, avgRDM(crun,:,:,smo), userOptions);
        end
        
    end
end

%% Analyse in scanner behaviour

cd('./behavioural_data')
control_response_averages = [];
control_rt_averages = [];
control_rt_medians = [];
control_AFCs = [];
patient_response_averages = [];
patient_rt_averages = [];
patient_rt_medians = [];
patient_AFCs = [];
graph_this = 0;
for crun = 1:length(subjects)
    
    [all_response_averages(crun,:), all_rt_averages(crun,:), all_rt_medians(crun,:), AFCs(crun,:)] = AFC_graph_this_subject_2021(subjects{crun}, dates{crun}, graph_this);
    
    if group(crun) == 1 % Controls
        control_response_averages(end+1,:) = all_response_averages(crun,:);
        control_rt_averages(end+1,:) = all_rt_averages(crun,:);
        control_rt_medians(end+1,:) = all_rt_medians(crun,:);
        control_AFCs(end+1,:) = AFCs(crun,:);
    elseif group(crun) == 2 % Patients
        patient_response_averages(end+1,:) = all_response_averages(crun,:);
        patient_rt_averages(end+1,:) = all_rt_averages(crun,:);
        patient_rt_medians(end+1,:) = all_rt_medians(crun,:);
        patient_AFCs(end+1,:) = AFCs(crun,:);
    end
end
cd('../')

figure
%cmap = colormap(parula(2));
colormap(jet)
set(gcf,'Position',[100 100 1600 800]);
subplot(3,1,1)
hold on
for this_x = 1:size(patient_response_averages,2)
    scatter(repmat(this_x+0.1,size(patient_response_averages,1),1)+rand(size(patient_response_averages,1),1)/30-(1/60),patient_response_averages(:,this_x),16,patient_AFCs/2)
    scatter(repmat(this_x-0.1,size(control_response_averages,1),1)+rand(size(control_response_averages,1),1)/30-(1/60),control_response_averages(:,this_x),16,control_AFCs/2)
end
xlim([0 6])
ylim([0 100])
title('Percent Correct')
set(gca,'XTick',[0:1:6])
set(gca,'XTickLabel',{'','Match 3','Mismatch 3','','Match 15','Mismatch 15',''},'XTickLabelRotation',15)

subplot(3,1,2)
hold on
for this_x = 1:size(patient_response_averages,2)
    scatter(repmat(this_x+0.1,size(patient_rt_averages,1),1)+rand(size(patient_rt_averages,1),1)/30-(1/60),patient_rt_averages(:,this_x),16,patient_AFCs/2+1)
    scatter(repmat(this_x-0.1,size(control_rt_averages,1),1)+rand(size(control_rt_averages,1),1)/30-(1/60),control_rt_averages(:,this_x),16,control_AFCs/2+1)
end
xlim([0 6])
title('Mean RT')
set(gca,'XTick',[0:1:6])
set(gca,'XTickLabel',{'','Match 3','Mismatch 3','','Match 15','Mismatch 15',''},'XTickLabelRotation',15)

subplot(3,1,3)
hold on
for this_x = 1:size(patient_rt_medians,2)
    scatter(repmat(this_x+0.1,size(patient_rt_medians,1),1)+rand(size(patient_rt_medians,1),1)/30-(1/60),patient_rt_medians(:,this_x),16,patient_AFCs/2)
    scatter(repmat(this_x-0.1,size(control_rt_medians,1),1)+rand(size(control_rt_medians,1),1)/30-(1/60),control_rt_medians(:,this_x),16,control_AFCs/2)
end
xlim([0 6])
title('Median RT')
set(gca,'XTick',[0:1:6])
set(gca,'XTickLabel',{'','Match 3','Mismatch 3','','Match 15','Mismatch 15',''},'XTickLabelRotation',15)

