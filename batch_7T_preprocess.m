% Batch script for preprocessing of pilot 7T data
% Written by TEC Feb 2018 and updated 2021
%
%% Setup environment
clear all
rmpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/'))
%addpath /imaging/local/software/spm_cbu_svn/releases/spm12_fil_r6906
addpath /group/language/data/thomascope/spm12_fil_r6906/
spm fmri

%% Define parameters
setup_file = 'PINFA_subjects_parameters';
eval(setup_file)
tr=2.5;
scriptdir = '/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/';

%% Options to skip steps
applytopup = 1;
opennewanalysispool = 1;

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

cat12normaliseworkedcorrectly = zeros(1,nrun);
jobs = repmat(jobfile, 1, 1);

parfor crun = 1:nrun
    spm('defaults', 'fMRI');
    spm_jobman('initcfg')
    try
        spm_jobman('run', jobs, inputs{:,crun});
        cat12normaliseworkedcorrectly(crun) = 1;
    catch
        cat12normaliseworkedcorrectly(crun) = 0;
    end
end

if ~all(cat12normaliseworkedcorrectly)
    error('failed at cat12 normalisation');
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
    jobfile{3} = {[scriptdir 'module_univariate_3runs_noneutral_lowthresh_job.m']};
    jobfile{4} = {[scriptdir 'module_univariate_4runs_noneutral_lowthresh_job.m']};
    inputs = cell(0, nrun);
    
    for crun = 1:nrun
        theseepis = find(strncmp(blocksout{crun},'Run',3));
        outpath = [preprocessedpathstem subjects{crun} '/'];
        filestoanalyse = cell(1,length(theseepis));
        
        tempDesign = module_get_event_times_AFC4(subjects{crun},dates{crun},length(theseepis),minvols(crun));
        
        inputs{1, crun} = cellstr([outpath 'stats6_' num2str(this_smooth)]);
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
        'con_0045.nii','Clarity Congruency Interaction Positive'
        'con_0050.nii','Clarity Congruency Interaction Negative'}; 
    expected_sessions = 4;
    
    visual_check = 0;
    nrun = size(all_conditions,1); % enter the number of runs here
    %jobfile = {'/group/language/data/thomascope/vespa/SPM12version/Standalone preprocessing pipeline/tc_source/batch_forwardmodel_job_noheadpoints.m'};
    
    this_scan = {};
    this_t_scan = {};
    firstlevel_folder = ['stats6_' num2str(this_smooth)];
    
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
%     inputs{1, crun} = cellstr([outpath 'stats5_multi_3']);
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
%         if any(cellfun(@isempty,inputs(1:(100*(sess-1))+101,crun))) % In case of trunkated run where an event did not occur, put it at the very end of the run so it isn't modelled but SPM doesn't crash
%             inputs{find(cellfun(@isempty,inputs(1:(100*(sess-1))+101,crun))),crun} = tr*(length(filestoanalyse{sess})-1);
%         end
%         %Now catch the cases where the task has continued beyond the end of a
%         %trunkated scan
%         for i = 3:size(inputs(1:(100*(sess-1))+101,crun),1)
%             if isnumeric(inputs{i,crun}) && all(inputs{i,crun}>tr*(length(filestoanalyse{sess})-1))
%                 inputs{i,crun} = tr*(length(filestoanalyse{sess})-1);
%             end
%         end
%     end
%     
%     jobs{crun} = jobfile{length(theseepis)};
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
% % Now repeat without the absent sound
% nrun = size(subjects,2); % enter the number of runs here
% jobfile = {};
% jobfile{3} = {[scriptdir 'module_univariate_3runs_noabsent_job.m']};
% jobfile{4} = {[scriptdir 'module_univariate_4runs_noabsent_job.m']};
% inputs = cell(0, nrun);
% 
% for crun = 1:nrun
%     theseepis = find(strncmp(blocksout{crun},'Run',3));
%     outpath = [preprocessedpathstem subjects{crun} '/'];
%     filestoanalyse = cell(1,length(theseepis));
%     
%     tempDesign = module_get_complex_event_times_AFC4(subjects{crun},dates{crun},length(theseepis),minvols(crun));
%     
%     inputs{1, crun} = cellstr([outpath 'stats5_multi_3_noabsent']);
%     for sess = 1:length(theseepis)
%         filestoanalyse{sess} = spm_select('ExtFPList',outpath,['^s3rtopup_' blocksin{crun}{theseepis(sess)}],1:minvols(crun));
%         inputs{(99*(sess-1))+2, crun} = cellstr(filestoanalyse{sess});
%         for cond_num = 1:80
%             inputs{(99*(sess-1))+2+cond_num, crun} = cat(2, tempDesign{sess}{cond_num})';
%         end
%         for cond_num = 81:96 %Response trials
%             inputs{(99*(sess-1))+2+cond_num, crun} = cat(2, tempDesign{sess}{cond_num+32})';
%         end
%         for cond_num = 97 %Button press
%             inputs{(99*(sess-1))+2+cond_num, crun} = cat(2, tempDesign{sess}{81})';
%         end
%         inputs{(99*(sess-1))+100, crun} = cellstr([outpath 'rp_topup_' blocksin{crun}{theseepis(sess)}(1:end-4) '.txt']);
%         if any(cellfun(@isempty,inputs(1:(99*(sess-1))+100,crun))) % In case of trunkated run where an event did not occur, put it at the very end of the run so it isn't modelled but SPM doesn't crash
%             inputs{find(cellfun(@isempty,inputs(1:(99*(sess-1))+100,crun))),crun} = tr*(length(filestoanalyse{sess})-1);
%         end
%         jobs{crun} = jobfile{length(theseepis)};
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
% %
% % %Now repeat with 8mm smoothing
% % 
% % for crun = 1:nrun
% %     theseepis = find(strncmp(blocksout{crun},'Run',3));
% %     outpath = [preprocessedpathstem subjects{crun} '/'];
% %     filestoanalyse = cell(1,length(theseepis));
% %     
% %     tempDesign = module_get_complex_event_times_AFC4(subjects{crun},dates{crun},length(theseepis),minvols(crun));
% %     
% %     inputs{1, crun} = cellstr([outpath 'stats5_multi_8']);
% %     for sess = 1:length(theseepis)
% %         filestoanalyse{sess} = spm_select('ExtFPList',outpath,['^s8rtopup_' blocksin{crun}{theseepis(sess)}],1:minvols(crun));
% %         inputs{(100*(sess-1))+2, crun} = cellstr(filestoanalyse{sess});
% %         for cond_num = 1:80
% %             inputs{(100*(sess-1))+2+cond_num, crun} = cat(2, tempDesign{sess}{cond_num})';
% %         end
% %         for cond_num = 81:96 %Response trials
% %             inputs{(100*(sess-1))+2+cond_num, crun} = cat(2, tempDesign{sess}{cond_num+32})';
% %         end
% %         for cond_num = 97 %Button press
% %             inputs{(100*(sess-1))+2+cond_num, crun} = cat(2, tempDesign{sess}{81})';
% %         end
% %         for cond_num = 98 %Absent sound (written only)
% %             inputs{(100*(sess-1))+2+cond_num, crun} = cat(2, tempDesign{sess}{129})';
% %         end
% %         inputs{(100*(sess-1))+101, crun} = cellstr([outpath 'rp_topup_' blocksin{crun}{theseepis(sess)}(1:end-4) '.txt']);
% %     end
% %     jobs{crun} = jobfile{length(theseepis)};
% %     if any(cellfun(@isempty,inputs(:,crun))) % In case of trunkated run where an event did not occur, put it at the very end of the run so it isn't modelled but SPM doesn't crash
% %         inputs{find(cellfun(@isempty,inputs(:,crun))),crun} = tr*(length(filestoanalyse{sess})-1);
% %     end
% % end
% % 
% % SPMworkedcorrectly = zeros(1,nrun);
% % parfor crun = 1:nrun
% %     spm('defaults', 'fMRI');
% %     spm_jobman('initcfg')
% %     try
% %         spm_jobman('run', jobs{crun}, inputs{:,crun});
% %         SPMworkedcorrectly(crun) = 1;
% %     catch
% %         SPMworkedcorrectly(crun) = 0;
% %     end
% % end
% 
% % %% Now create a more complex SPM with variable levels of AR whitening, with word omissions specified
% %
% % jobfile = {};
% % jobfile{3} = {[scriptdir 'module_univariate_3runs_complex_AR_job.m']};
% % jobfile{4} = {[scriptdir 'module_univariate_4runs_complex_AR_job.m']};
% %
% % all_aros = [1 3 6 12]; %Autoregressive model order
% % nrun = size(subjects,2)*length(all_aros); % enter the number of runs here
% % inputs = cell(0, size(subjects,2),length(all_aros));
% % for this_aro = 1:length(all_aros);
% % for crun = 1:size(subjects,2)
% %     aro = all_aros(this_aro);
% %
% %     theseepis = find(strncmp(blocksout{crun},'Run',3));
% %     outpath = [preprocessedpathstem subjects{crun} '/'];
% %     filestoanalyse = cell(1,length(theseepis));
% %
% %     tempDesign = module_get_complex_event_times(subjects{crun},dates{crun},length(theseepis),minvols(crun));
% %
% %     inputs{1, crun, this_aro} = cellstr([outpath 'stats5_multi_AR' num2str(aro)]);
% %     for sess = 1:length(theseepis)
% %         filestoanalyse{sess} = spm_select('ExtFPList',outpath,['^s3topup_' blocksin{crun}{theseepis(sess)}],1:minvols(crun));
% %         inputs{(100*(sess-1))+2, crun, this_aro} = cellstr(filestoanalyse{sess});
% %         for cond_num = 1:80
% %             inputs{(100*(sess-1))+2+cond_num, crun, this_aro} = cat(2, tempDesign{sess}{cond_num})';
% %         end
% %         for cond_num = 81:96 %Response trials
% %             inputs{(100*(sess-1))+2+cond_num, crun, this_aro} = cat(2, tempDesign{sess}{cond_num+32})';
% %         end
% %         for cond_num = 97 %Button press
% %             inputs{(100*(sess-1))+2+cond_num, crun, this_aro} = cat(2, tempDesign{sess}{81})';
% %         end
% %         for cond_num = 98 %Absent sound (written only)
% %             inputs{(100*(sess-1))+2+cond_num, crun, this_aro} = cat(2, tempDesign{sess}{129})';
% %         end
% %         inputs{(100*(sess-1))+101, crun, this_aro} = cellstr([outpath 'rp_topup_' blocksin{crun}{theseepis(sess)}(1:end-4) '.txt']);
% %     end
% %     %inputs{(100*(sess-1))+102, crun, this_aro} = 'AR(1)';
% %     inputs{(100*(sess-1))+102, crun, this_aro} = aro;
% %     jobs{crun} = jobfile{length(theseepis)};
% %
% % end
% % end
% %
% %
% % try
% %     matlabpool 'close'
% % catch
% %     delete(gcp)
% % end
% %
% %
% % workersrequested = 24;
% % workerpool = cbupool(workersrequested);
% % workerpool.ResourceTemplate=['-l nodes=^N^,mem=768GB,walltime=168:00:00'];
% % try
% %     matlabpool(workerpool)
% % catch
% %     parpool(workerpool,workerpool.NumWorkers)
% % end
% %
% % all_combs = combvec(1:size(subjects,2),1:length(all_aros))';
% % SPMworkedcorrectly = zeros(1,size(all_combs,1));
% % for thisone = 1:size(all_combs,1)
% %     crun = all_combs(thisone,1);
% %     this_aro = all_combs(thisone,2);
% %     spm('defaults', 'fMRI');
% %     spm_jobman('initcfg')
% %     try
% %         spm_jobman('run', jobs{crun}, inputs{:,crun, this_aro});
% %         SPMworkedcorrectly(thisone) = 1;
% %     catch
% %         SPMworkedcorrectly(thisone) = 0;
% %     end
% % end
 
%% Now create a ReML SPM without modelling the written word separately for future multivariate analysis (currently only implemented for 3 or 4 runs) - Native space (s3r)
nrun = size(subjects,2); % enter the number of runs here
jobfile = {};
jobfile{3} = {[scriptdir 'module_univariate_3runs_noabsent_lowthresh_job.m']};
jobfile{4} = {[scriptdir 'module_univariate_4runs_noabsent_lowthresh_job.m']};
inputs = cell(0, nrun);

for crun = 1:nrun
    theseepis = find(strncmp(blocksout{crun},'Run',3));
    outpath = [preprocessedpathstem subjects{crun} '/'];
    filestoanalyse = cell(1,length(theseepis));
    
    tempDesign = module_get_complex_event_times_nowritten_AFC4(subjects{crun},dates{crun},length(theseepis),minvols(crun));
    
    inputs{1, crun} = cellstr([outpath 'stats5_multi_3_nowritten2_midthresh']);
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
jobfile{3} = {[scriptdir 'module_univariate_3runs_noabsent_lowthresh_job.m']};
jobfile{4} = {[scriptdir 'module_univariate_4runs_noabsent_lowthresh_job.m']};
inputs = cell(0, nrun);

for crun = 1:nrun
    theseepis = find(strncmp(blocksout{crun},'Run',3));
    outpath = [preprocessedpathstem subjects{crun} '/'];
    filestoanalyse = cell(1,length(theseepis));
    
    tempDesign = module_get_complex_event_times_nowritten_AFC4(subjects{crun},dates{crun},length(theseepis),minvols(crun));
    
    inputs{1, crun} = cellstr([outpath 'stats5_multi_8_nowritten2']);
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

%% Repeat with an explicit grey matter mask - Now create a ReML SPM without modelling the written word separately for future multivariate analysis (currently only implemented for 3 or 4 runs) - Native space (s3r)
% This did not work very well, I think because of movement meaning that important voxels were not retained in the mask
nrun = size(subjects,2); % enter the number of runs here
jobfile = {};
jobfile{3} = {[scriptdir 'module_univariate_3runs_noabsent_lowthresh_masked_job.m']};
jobfile{4} = {[scriptdir 'module_univariate_4runs_noabsent_lowthresh_masked_job.m']};
inputs = cell(0, nrun);

redo = 1;

for crun = 1:nrun
    theseepis = find(strncmp(blocksout{crun},'Run',3));
    outpath = [preprocessedpathstem subjects{crun} '/'];
    filestoanalyse = cell(1,length(theseepis));
    
    tempDesign = module_get_complex_event_times_nowritten_AFC4(subjects{crun},dates{crun},length(theseepis),minvols(crun));
    
    inputs{1, crun} = cellstr([outpath 'stats5_multi_3_nowritten2_masked']);
    if redo && exist(char(inputs{1, crun}));
        delete([char(inputs{1, crun}) '/SPM.mat'])
    end
    
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
    inputs{(99*(sess-1))+101,crun} = cellstr([outpath 's5_structural_mask_0.05.nii']);
    if ~exist(inputs{(99*(sess-1))+101,crun}{1},'file') %Make a grey matter native space mask thresholded at 5%, mainly to exclude muscle
        spm_imcalc([outpath 'structural.nii'],[outpath 'temp.nii'],'i1>0.05')
        spm_smooth([outpath 'temp.nii'],[outpath 'temp2.nii'],[5 5 5])
        spm_imcalc([outpath 'temp2.nii'],[outpath 's5_structural_mask_0.05.nii'],'i1>0.05')
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


%% Now run the cross validated Mahalanobis distance and RSM on each subject on the whole brain downsampled at 2 (quick)
nrun = size(subjects,2); % enter the number of runs here
mahalanobisworkedcorrectly = zeros(1,nrun);
downsamp_ratio = 2; %Downsampling in each dimension, must be an integer, 2 is 8 times faster than 1 (2 cubed). 
parfor crun = 1:nrun
    addpath(genpath('./RSA_scripts'))
    GLMDir = [preprocessedpathstem subjects{crun} '/stats5_multi_3_nowritten2'];
    try
        TDTCrossnobisAnalysis_1Subj(GLMDir,downsamp_ratio)
        mahalanobisworkedcorrectly(crun) = 1;
    catch
        mahalanobisworkedcorrectly(crun) = 0;
    end
end

%% Or run the cross validated Mahalanobis distance and RSM on each subject on the whole brain not downsampled, but in parallel over voxels (slow, and produces around 12Gb output data per subject, but more powerful statistics) - the bigger the worker pool the better.
run_not_downsampled = 0; % NB: MAKES HUGE FILES!
if run_not_downsampled
    nrun = size(subjects,2); % enter the number of runs here
    mahalanobisparallelworkedcorrectly = zeros(1,nrun);
    if opennewanalysispool == 1
        delete(gcp) % Make a bigger pool for this step.
        Poolinfo = cbupool(120,'--mem-per-cpu=1G --time=167:00:00 --exclude=node-i[01-15]');
        parpool(Poolinfo,Poolinfo.NumWorkers);
    end
    for crun = 1:nrun
        if numel(gcp('nocreate')) == 0 % If parallel pool crashes, this should allow the loop to simply resume at the next subject
            Poolinfo = cbupool(60,'--mem-per-cpu=1G --time=167:00:00 --exclude=node-i[01-15]');
            parpool(Poolinfo,Poolinfo.NumWorkers);
        end
        addpath(genpath('./RSA_scripts'))
        GLMDir = [preprocessedpathstem subjects{crun} '/stats5_multi_3_nowritten2'];
        try
            TDTCrossnobisAnalysis_parallelsearch(GLMDir)
            mahalanobisparallelworkedcorrectly(crun) = 1;
        catch
            mahalanobisparallelworkedcorrectly(crun) = 0;
        end
    end
    if opennewanalysispool == 1
        delete(gcp)
        if size(subjects,2) > 64
            workersrequested = 64;
            fprintf([ '\n\nUnable to ask for a worker per run; asking for 64 instead\n\n' ]);
        else
            workersrequested = size(subjects,2);
        end
        Poolinfo = cbupool(workersrequested,'--mem-per-cpu=12G --time=167:00:00 --exclude=node-i[01-15]');
        parpool(Poolinfo,Poolinfo.NumWorkers);
    end
end


%% Do an RSA analysis separately if you want (already integrated into previous step for vowels, but now can compare new models etc without repeating the time consuming cross-nobis)
nrun = size(subjects,2); % enter the number of runs here
RSAnobisworkedcorrectly = zeros(1,nrun);
downsamp_ratio = 2; %Downsampling in each dimension, must be an integer, 2 is 8 times faster than 1 (2 cubed). 
parfor crun = 1:nrun
    addpath(genpath('./RSA_scripts'))
    GLMDir = [preprocessedpathstem subjects{crun} '/stats4_multi_3_nowritten2'];
    try
        module_make_effect_maps(GLMDir,downsamp_ratio)
        RSAnobisworkedcorrectly(crun) = 1;
    catch
        RSAnobisworkedcorrectly(crun) = 0;
    end
end

%% Do a partial-correlation based RSA analysis to tease apart the written and spoken word representations
nrun = size(subjects,2); % enter the number of runs here
partialRSAnobisworkedcorrectly = zeros(1,nrun);
downsamp_ratio = 2; %Downsampling in each dimension, must be an integer, 2 is 8 times faster than 1 (2 cubed). 
parfor crun = 1:nrun
    addpath(genpath('./RSA_scripts'))
    GLMDir = [preprocessedpathstem subjects{crun} '/stats4_multi_3_nowritten2'];
    try
        module_make_partial_effect_maps(GLMDir,downsamp_ratio)
        partialRSAnobisworkedcorrectly(crun) = 1;
    catch
        partialRSAnobisworkedcorrectly(crun) = 0;
    end
end

%% Now normalise the native space RSA maps into template space with CAT12 deformation fields calculated earlier
nrun = size(subjects,2); % enter the number of runs here
native2templateworkedcorrectly = zeros(1,nrun);
downsamp_ratio = 2; %Downsampling in each dimension, must be an integer, 2 is 8 times faster than 1 (2 cubed). 
parfor crun = 1:nrun
    addpath(genpath('./RSA_scripts'))
    GLMDir = [preprocessedpathstem subjects{crun} '/stats4_multi_3_nowritten2'];
    outpath = [preprocessedpathstem subjects{crun} '/'];
    try
        module_nativemap_2_template(GLMDir,downsamp_ratio,outpath)
        native2templateworkedcorrectly(crun) = 1;
    catch
        native2templateworkedcorrectly(crun) = 0;
    end
end

%% Now do a second level analysis on the searchlights
crun = 1;
age_lookup = readtable('Pinfa_ages.csv');
downsamp_ratio = 2; %Downsampling in each dimension, must be an integer, 2 is 8 times faster than 1 (2 cubed). 
rmpath('/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/RSA_scripts/es_scripts_fMRI') %Stops SPM getting defaults for second level if on path

GLMDir = [preprocessedpathstem subjects{crun} '/stats4_multi_3_nowritten2']; %Template, first subject
outpath = [preprocessedpathstem '/stats4_multi_3_nowritten2/searchlight/downsamp_' num2str(downsamp_ratio) filesep 'second_level']; %Results directory

searchlightsecondlevel = [];
searchlightsecondlevel = module_searchlight_secondlevel(GLMDir,subjects,group,age_lookup,outpath,downsamp_ratio);

%% Now do a correlation analysis with the Bayesian perceptual model parameters against selected models
model_run_date = '13-May-2021';
try
    load(['./modelparameters/modelparameters_' model_run_date '.mat'])
catch
    [all_sigma_pred,all_thresholds,controls_sigma_pred,controls_threshold,patients_sigma_pred,patients_threshold] = module_bayesian_behaviour(subjects,group,dates);
    save(['./modelparameters/modelparameters_' date '.mat'],'all_sigma_pred','all_thresholds','controls_sigma_pred','controls_threshold','patients_sigma_pred','patients_threshold');
end

downsamp_ratio = 1;
age_lookup = readtable('Pinfa_ages.csv');

Basefilepath = [preprocessedpathstem subjects{1} '/stats4_multi_3_nowritten2/TDTcrossnobis/spearman/weffect-map_']; %Template, first subject
outpath = [preprocessedpathstem '/stats4_multi_3_nowritten2/searchlight/downsamp_' num2str(downsamp_ratio) filesep 'correlations' filesep 'prediction_precision' filesep]; %Results directory

condition_names = {
  'M to MM Shared Segments:  Cross Negative partialling '
  'Match Clear to Mismatch Unclear only cross.nii'
};

covariatesecondlevelworkedcorrectly = zeros(1,size(condition_names,1));
for crun = 1:length(condition_names)
    thisfilepath = [Basefilepath condition_names{crun} '.nii'];
    thisoutpath = [outpath condition_names{crun}];
    covariatesecondlevelworkedcorrectly(crun) = module_multivariate_correlation(nanmean(all_sigma_pred),thisfilepath,subjects,age_lookup,thisoutpath);
    % cd(deblank(thisoutpath))
    % pause % If you want to view the outputs.
end

XXX WIP

%% Now normalise the template space masks into native space

images2normalise = {%'/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/RSA_scripts/Blank_ROI/blank_mask.nii'
    %'/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/atlas_Neuromorphometrics/Blank_2016_inflated.nii' %Blank and Davis 2018 mask
    %     '/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/atlas_Neuromorphometrics/Left_IFG_Written_Cluster.nii' %Cross-decoding Match unclear to Mismatch unclear
    %     '/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/atlas_Neuromorphometrics/Left_Precentral_Written_Cluster.nii' %Cross-decoding Match unclear to Mismatch unclear
    %     '/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/atlas_Neuromorphometrics/Left_IFG_cross_group_cluster.nii' %M to MM Shared Segments:  Cross Negative partialling
%     '/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/atlas_Neuromorphometrics/Left_Frontal_Univariate_MM>M.nii'
%     '/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/atlas_Neuromorphometrics/Left_Temporal_Univariate_MM>M.nii'
%     '/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/atlas_Neuromorphometrics/Left_PostSTG_Univariate_Interaction.nii'
%     '/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/atlas_Neuromorphometrics/Left_Precentral_Univariate_Interaction1.nii'
%     '/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/atlas_Neuromorphometrics/Left_Precentral_Univariate_Interaction2.nii'
%     '/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/atlas_Neuromorphometrics/Left_Precentral_Univariate_Interaction3.nii'
%     '/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/atlas_Neuromorphometrics/Left_Angular_Univariate_Interaction1.nii'
%     '/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/atlas_Neuromorphometrics/Left_Angular_Univariate_Interaction2.nii'
%     '/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/atlas_Neuromorphometrics/Left_Precentral_Univariate_Interaction_combined.nii'
%     '/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/atlas_Neuromorphometrics/Left_Angular_Univariate_Interaction_combined.nii'
%     '/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/atlas_Neuromorphometrics/Left_STG_Univariate8mm_15>3.nii'
%     '/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/atlas_Neuromorphometrics/Left_STG_Univariate3mm_15>3.nii'
%     '/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/atlas_Neuromorphometrics/Left_PrG_SSMatchnoself_combined.nii'
    '/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/atlas_Neuromorphometrics/Left_PrG_All_Shared_Segments.nii'
    };

% search_labels = {
%     'Left STG'
%     'Left PT'
%     'Left PrG'
%     'Left FO'
%     'Left TrIFG'
%     };

% xA=spm_atlas('load','Neuromorphometrics');

search_labels = {
%     'Left Superior Temporal Gyrus'
%     'Left Angular Gyrus'
%     'Left Precentral Gyrus'
%     'Left Frontal Operculum'
%     'Left Inferior Frontal Angular Gyrus'
%     'Right Superior Temporal Gyrus'
%     'Right Angular Gyrus'
%     'Right Precentral Gyrus'
%     'Right Frontal Operculum'
%     'Right Inferior Frontal Angular Gyrus'
%     'Left Cerebellar Lobule Cerebellar Vermal Lobules VI-VII'
%     'Right Cerebellar Lobule Cerebellar Vermal Lobules VI-VII'
    };

% cat_install_atlases
% xA=spm_atlas('load','dartel_neuromorphometrics');
% for i = 1:size(xA.labels,2)
%     all_labels{i} = xA.labels(i).name;
% end

% S = cell(1,length(search_labels));
% for i = 1:length(search_labels)
%     S{i} = find(strncmp(all_labels,search_labels{i},size(search_labels{i},2)));
% end
% if ~exist('./atlas_Neuromorphometrics/','dir')
%     mkdir('./atlas_Neuromorphometrics/');
% end
% for i = 1:size(S,2)
%     fname=strcat(strrep(search_labels{i}, ' ', '_'),'.nii');
%     VM=spm_atlas('mask',xA,xA.labels(S{i}).name);
%     VM.fname=['./atlas_Neuromorphometrics/' fname];
%     spm_write_vol(VM,spm_read_vols(VM));
%     images2normalise{end+1} = [pwd '/atlas_Neuromorphometrics/' fname];
% end

nrun = size(subjects,2); % enter the number of runs here
template2nativeworkedcorrectly = zeros(1,nrun);
parfor crun = 1:nrun
    addpath(genpath('./RSA_scripts'))
    outpath = [preprocessedpathstem subjects{crun} '/'];
    reslice_template = [preprocessedpathstem subjects{crun} '/stats4_multi_3_nowritten2/mask.nii']; %Template for reslicing
    try
        module_template_2_nativemap(images2normalise,outpath,1,reslice_template);
        template2nativeworkedcorrectly(crun) = 1;
    catch
        template2nativeworkedcorrectly(crun) = 0;
    end
end

%% Analyse by condition and brain region
addpath(genpath('/imaging/mlr/users/tc02/toolboxes')); %Where is the RSA toolbox?

% masks = { %rw for re-sliced after warping into native space (I have also re-binarised). Underscores for spaces in the atlas search above.
%     'rwblank_mask'
%     'rwLeft_STG'
%     'rwLeft_PT'
%     'rwLeft_PrG'
%     'rwLeft_FO'
%     'rwLeft_TrIFG'
%     };

masks = {
    'rwLeft_PrG_All_Shared_Segments'
%'rwLeft_PrG_SSMatchnoself_combined'
%     'rwLeft_STG_Univariate8mm_15>3'
%     'rwLeft_STG_Univariate3mm_15>3'
%     'rwLeft_Precentral_Univariate_Interaction_combined'
%     'rwLeft_Angular_Univariate_Interaction_combined'
%         'rwLeft_PostSTG_Univariate_Interaction'
%     'rwLeft_Precentral_Univariate_Interaction1'
%     'rwLeft_Precentral_Univariate_Interaction2'
%     'rwLeft_Precentral_Univariate_Interaction3'
%     'rwLeft_Angular_Univariate_Interaction1'
%     'rwLeft_Angular_Univariate_Interaction2'
%     'rwLeft_Frontal_Univariate_MM>M'
%     'rwLeft_Temporal_Univariate_MM>M'
%     'rwLeft_IFG_cross_group_cluster'
%     'rwBlank_2016_inflated'
%     'rwL_STG_cross-segment_cluster'
%     'rwLeft_Superior_Temporal_Gyrus'
%     'rwLeft_Angular_Gyrus'
%     'rwLeft_Precentral_Gyrus'
%     'rwLeft_Frontal_Operculum'
%     'rwLeft_Inferior_Frontal_Angular_Gyrus'
%     'rwRight_Superior_Temporal_Gyrus'
%     'rwRight_Angular_Gyrus'
%     'rwRight_Precentral_Gyrus'
%     'rwRight_Frontal_Operculum'
%     'rwRight_Inferior_Frontal_Angular_Gyrus'
%     'rwLeft_IFG_Written_Cluster'
%     'rwLeft_Precentral_Written_Cluster'
    };

GLMDir = [preprocessedpathstem subjects{1} '/stats4_multi_3_nowritten2']; %Template, first subject
temp = load([GLMDir filesep 'SPM.mat']);
labelnames = {};
for i = 1:length(temp.SPM.Sess(1).U)
    if ~strncmp(temp.SPM.Sess(1).U(i).name,{'Match','Mismatch','Written'},5)
        continue
    else
        labelnames(end+1) = temp.SPM.Sess(1).U(i).name;
    end
end        
labelnames_denumbered = {};
for i = 1:length(labelnames)
    labelnames_denumbered{i} = labelnames{i}(isletter(labelnames{i})|isspace(labelnames{i}));
end
conditionnames = unique(labelnames_denumbered,'stable');
clear temp labelnames_denumbered labelnames

nrun = size(subjects,2); % enter the number of runs here
mahalanobisroiworkedcorrectly = zeros(1,nrun);
parfor crun = 1:nrun
    addpath(genpath('./RSA_scripts'))
    GLMDir = [preprocessedpathstem subjects{crun} '/stats4_multi_3_nowritten2']; %Where is the SPM model?
    mask_dir = [preprocessedpathstem subjects{crun}]; %Where are the native space ROI masks?
    try
        TDTCrossnobisAnalysis_roi(GLMDir,mask_dir,masks);
        mahalanobisroiworkedcorrectly(crun) = 1;
    catch
        mahalanobisroiworkedcorrectly(crun) = 0;
    end
end

%% Now do RSA on ROI data
nrun = size(subjects,2); % enter the number of runs here
RSAroiworkedcorrectly = zeros(1,nrun);
partialRSAroiworkedcorrectly = zeros(1,nrun);
masks = {
'rwLeft_PrG_All_Shared_Segments'
%     'rwLeft_PrG_SSMatchnoself_combined'
%     'rwLeft_STG_Univariate8mm_15>3'
%     'rwLeft_STG_Univariate3mm_15>3'
    %     'rwLeft_Precentral_Univariate_Interaction_combined'
    %     'rwLeft_Angular_Univariate_Interaction_combined'
    %     'rwLeft_PostSTG_Univariate_Interaction'
    %     'rwLeft_Precentral_Univariate_Interaction1'
    %     'rwLeft_Precentral_Univariate_Interaction2'
    %     'rwLeft_Precentral_Univariate_Interaction3'
%     'rwLeft_Angular_Univariate_Interaction1'
%     'rwLeft_Angular_Univariate_Interaction2'
%     'rwLeft_Frontal_Univariate_MM>M'
%     'rwLeft_Temporal_Univariate_MM>M'
%     'rwLeft_IFG_cross_group_cluster'
%     'rwBlank_2016_inflated'
%     'rwL_STG_cross-segment_cluster'
%     'rwLeft_Superior_Temporal_Gyrus'
%     'rwLeft_Angular_Gyrus'
%     'rwLeft_Precentral_Gyrus'
%     'rwLeft_Frontal_Operculum'
%     'rwLeft_Inferior_Frontal_Angular_Gyrus'
%     'rwRight_Superior_Temporal_Gyrus'
%     'rwRight_Angular_Gyrus'
%     'rwRight_Precentral_Gyrus'
%     'rwRight_Frontal_Operculum'
%     'rwRight_Inferior_Frontal_Angular_Gyrus'
%     'rwLeft_IFG_Written_Cluster'
%     'rwLeft_Precentral_Written_Cluster'
    };

parfor crun = 1:nrun
    addpath(genpath('./RSA_scripts'))
    GLMDir = [preprocessedpathstem subjects{crun} '/stats4_multi_3_nowritten2']; %Where is the SPM model?
    try
        module_roi_RSA(GLMDir,masks)
        RSAroiworkedcorrectly(crun) = 1;
    catch
        RSAroiworkedcorrectly(crun) = 0;
    end
end
parfor crun = 1:nrun
    addpath(genpath('./RSA_scripts'))
    GLMDir = [preprocessedpathstem subjects{crun} '/stats4_multi_3_nowritten2']; %Where is the SPM model?
    try
        module_make_partial_roi_RSA(GLMDir,masks)
        partialRSAroiworkedcorrectly(crun) = 1;
    catch
        partialRSAroiworkedcorrectly(crun) = 0;
    end
end


%% Compare across conditions in STG as a sanity check then go on to do all ROIs
GLMDir = [preprocessedpathstem subjects{1} '/stats4_multi_3_nowritten2']; %Template, first subject
temp = load([GLMDir filesep 'SPM.mat']);
labelnames = {};
for i = 1:length(temp.SPM.Sess(1).U)
    if ~strncmp(temp.SPM.Sess(1).U(i).name,{'Match','Mismatch','Written'},5)
        continue
    else
        labelnames(end+1) = temp.SPM.Sess(1).U(i).name;
    end
end        
labelnames_denumbered = {};
for i = 1:length(labelnames)
    labelnames_denumbered{i} = labelnames{i}(isletter(labelnames{i})|isspace(labelnames{i}));
end
conditionnames = unique(labelnames_denumbered,'stable');

clear this_model_name mask_names
this_model_name{1} = {
%     'Match Unclear vowels'
%     'Match Clear vowels'
%     'Mismatch Unclear vowels'
%     'Mismatch Clear vowels'
    'Match Unclear shared_segments'
    'Match Clear shared_segments'
    'Mismatch Unclear shared_segments'
    'Mismatch Clear shared_segments'
%     'Written vowels'
    'Written shared_segments'
    };

this_model_name{2} = {
    'Match Unclear to Mismatch Unclear Cross-decode_Match'
    'Match Unclear to Mismatch Unclear SS_Match'
    'Match Unclear to Mismatch Unclear SS_Match - no self'
    'Match Unclear to Mismatch Clear Cross-decode_Match'
    'Match Unclear to Mismatch Clear SS_Match'
    'Match Unclear to Mismatch Clear SS_Match - no self'
    'Match Clear to Mismatch Unclear Cross-decode_Match'
    'Match Clear to Mismatch Unclear SS_Match'
    'Match Clear to Mismatch Unclear SS_Match - no self'
    'Match Clear to Mismatch Clear Cross-decode_Match'
    'Match Clear to Mismatch Clear SS_Match'
    'Match Clear to Mismatch Clear SS_Match - no self'
    'Match Unclear to Written Cross-decode_Match'
    'Match Clear to Written Cross-decode_Match'
    'Mismatch Unclear to Written Cross-decode_Match'
    'Mismatch Clear to Written Cross-decode_Match'
    };

this_model_name{3} = {'All spoken Cross-decode_Match'
    'All spoken SS_Match'
    'All spoken SS_Match - no self'
    'Spoken to Written Cross-decode_Match'
    'Spoken to Written SS_Match - no self'
    'Match to Mismatch Shared Segments - no self'
    'Match to Mismatch SS_Match - no self'
    'Match to Mismatch combined_SS - no self - rescaled'
    'Match to Mismatch only cross'
    'Match to Mismatch only not cross'
    };

this_model_name{4} = {
    'Match Unclear to Match Clear Cross-decode_Match';
    'Mismatch Unclear to Mismatch Clear Cross-decode';
    'Match Unclear to Match Clear Cross-decode';
    'Mismatch Unclear to Mismatch Clear Cross-decode_Match';
    };

this_model_name{5} = {
    'M to MM Shared Segments:  Cross Negative partialling '
    'M to MM Shared Segments:  partialling  Cross Negative';
    };


% this_model_name{6} = {
%     'Match Unclear to Mismatch Unclear Cross-decode'
%     'Match Unclear to Mismatch Unclear Shared Segments - cross'
%     'Match Unclear to Mismatch Unclear Shared Segments - no self'
%     'Match Unclear to Mismatch Clear Cross-decode'
%     'Match Unclear to Mismatch Clear Shared Segments - cross'
%     'Match Unclear to Mismatch Clear Shared Segments - no self'
%     'Match Clear to Mismatch Unclear Cross-decode'
%     'Match Clear to Mismatch Unclear Shared Segments - cross'
%     'Match Clear to Mismatch Unclear Shared Segments - no self'
%     'Match Clear to Mismatch Clear Cross-decode'
%     'Match Clear to Mismatch Clear Shared Segments - cross'
%     'Match Clear to Mismatch Clear Shared Segments - no self'
%     'Match Unclear to Written Cross-decode'
%     'Match Clear to Written Cross-decode'
%     'Mismatch Unclear to Written Cross-decode'
%     'Mismatch Clear to Written Cross-decode'
%     };

nrun = size(subjects,2); % enter the number of runs here
% First load in the similarities
RSA_ROI_data_exist = zeros(1,nrun);
all_data = [];
mask_names{1} = {
         'rwLeft_IFG_cross_group_cluster'
    %     'rwLeft_Superior_Temporal_Gyrus';
    'rwL_STG_cross-segment_cluster'
    %     'rwBlank_2016_inflated'
    'rwLeft_Frontal_Univariate_MM>M'
    'rwLeft_Temporal_Univariate_MM>M'
    %     'rwLeft_PostSTG_Univariate_Interaction'
    %     %     'rwLeft_Precentral_Univariate_Interaction1'
    %     %     'rwLeft_Precentral_Univariate_Interaction2'
    %     %     'rwLeft_Precentral_Univariate_Interaction3'
    %     'rwLeft_Angular_Univariate_Interaction1'
    %     'rwLeft_Angular_Univariate_Interaction2'
    'rwLeft_Precentral_Univariate_Interaction_combined'
    'rwLeft_Angular_Univariate_Interaction_combined'
    %     'rwLeft_STG_Univariate8mm_15>3'
    'rwLeft_STG_Univariate3mm_15>3'
    'rwLeft_PrG_SSMatchnoself_combined'
    'rwLeft_PrG_All_Shared_Segments'
    };

% mask_names{2} = {
%     %'rwL_STG_cross-segment_cluster'
%     'rwLeft_Angular_Gyrus'
%     'rwLeft_Precentral_Gyrus'
%     'rwLeft_Frontal_Operculum'
%     'rwLeft_Inferior_Frontal_Angular_Gyrus'
%     'rwRight_Superior_Temporal_Gyrus'
%     'rwRight_Angular_Gyrus'
%     'rwRight_Precentral_Gyrus'
%     'rwRight_Frontal_Operculum'
%     'rwRight_Inferior_Frontal_Angular_Gyrus'
%     };

for j = 1:length(this_model_name)
    for k = 1:length(mask_names)
        for i = 1:length(mask_names{k})
            all_data = [];
            for crun = 1:nrun
                %ROI_RSA_dir = [preprocessedpathstem subjects{crun} '/stats4_multi_3_nowritten2/TDTcrossnobis_ROI/RSA/spearman']; %Where are the results>
                ROI_RSA_dir = [preprocessedpathstem subjects{crun} '/stats4_multi_3_nowritten2/TDTcrossnobis_ROI/' mask_names{k}{i} '/RSA/spearman'];
                if ~exist(fullfile(ROI_RSA_dir,['roi_effects_' this_model_name{j}{1} '.mat']),'file')
                    ROI_RSA_dir = [preprocessedpathstem subjects{crun} '/stats4_multi_3_nowritten2/TDTcrossnobis_ROI' mask_names{k}{i} '/RSA/spearman']; % Stupid coding error earlier in analysis led to misnamed directories
                end
                for m = 1:length(this_model_name{j})
                    try
                        temp_data = load(fullfile(ROI_RSA_dir,['roi_effects_' this_model_name{j}{m} '.mat']));
                        all_data(m,:,crun) = temp_data.roi_effect; %Create a matrix of condition by ROI by subject
                        RSA_ROI_data_exist(crun) = 1;
                    catch
                        warning(['No data for ' subjects{crun} ' probably because of SPM dropout, ignoring them'])
                        %error
                        RSA_ROI_data_exist(crun) = 0;
                        continue
                    end
                end
            end
            roi_names = temp_data.roi_names;
            clear temp_data
            disp(['Excluding subjects ' num2str(find(RSA_ROI_data_exist==0)) ' belonging to groups ' num2str(group(RSA_ROI_data_exist==0)) ' maybe check them'])
            all_data(:,:,RSA_ROI_data_exist==0) = NaN;
            
            %this_ROI = find(strcmp('rwLeft_Superior_Temporal_Gyrus',roi_names));
            this_ROI = find(strcmp(mask_names{k}{i},roi_names));
            figure
            set(gcf,'Position',[100 100 1600 800]);
            set(gcf, 'PaperPositionMode', 'auto');
            hold on
            errorbar([1:length(this_model_name{j})]-0.1,nanmean(squeeze(all_data(:,this_ROI,group==1&RSA_ROI_data_exist)),2),nanstd(squeeze(all_data(:,this_ROI,group==1&RSA_ROI_data_exist))')/sqrt(sum(group==1&RSA_ROI_data_exist)),'kx')
            errorbar([1:length(this_model_name{j})]+0.1,nanmean(squeeze(all_data(:,this_ROI,group==2&RSA_ROI_data_exist)),2),nanstd(squeeze(all_data(:,this_ROI,group==2&RSA_ROI_data_exist))')/sqrt(sum(group==2&RSA_ROI_data_exist)),'rx')
            xlim([0 length(this_model_name{j})+1])
            set(gca,'xtick',[1:length(this_model_name{j})],'xticklabels',this_model_name{j},'XTickLabelRotation',45,'TickLabelInterpreter','none')
            plot([0 length(this_model_name{j})+1],[0,0],'k--')
            title([mask_names{k}{i}(3:end) ' RSA'],'Interpreter','none')
            if verLessThan('matlab', '9.2')
                legend('Controls','Patients','location','southeast')
            else
                legend('Controls','Patients','location','southeast','AutoUpdate','off')
            end
            [h,p] = ttest(squeeze(all_data(:,this_ROI,logical(RSA_ROI_data_exist)))');
            these_y_lims = ylim;
            if sum(h)~=0
                plot(find(h),these_y_lims(2)-diff(these_y_lims/10),'k*')
            end
            [h,p] = ttest2(squeeze(all_data(:,this_ROI,group==1&logical(RSA_ROI_data_exist)))',squeeze(all_data(:,this_ROI,group==2&logical(RSA_ROI_data_exist)))');
            if sum(h)~=0
                plot(find(h),these_y_lims(2)-diff(these_y_lims/20),'kx')
            end
            drawnow
        end
    end
end


%% Now compare across ROI for each condition - WORK IN PROGRESS
GLMDir = [preprocessedpathstem subjects{1} '/stats4_multi_3_nowritten2']; %Template, first subject
temp = load([GLMDir filesep 'SPM.mat']);
labelnames = {};
for i = 1:length(temp.SPM.Sess(1).U)
    if ~strncmp(temp.SPM.Sess(1).U(i).name,{'Match','Mismatch','Written'},5)
        continue
    else
        labelnames(end+1) = temp.SPM.Sess(1).U(i).name;
    end
end        
labelnames_denumbered = {};
for i = 1:length(labelnames)
    labelnames_denumbered{i} = labelnames{i}(isletter(labelnames{i})|isspace(labelnames{i}));
end
conditionnames = unique(labelnames_denumbered,'stable');

clear temp labelnames_denumbered labelnames

nrun = size(subjects,2); % enter the number of runs here
% First load in the similarities
RSA_ROI_data_exist = zeros(1,nrun);
all_data = [];
for crun = 1:nrun
    ROI_RSA_dir = [preprocessedpathstem subjects{crun} '/stats4_multi_3_nowritten2/TDTcrossnobis_ROI/RSA/spearman']; %Where are the results>
    for m = 1:length(conditionnames)
        try
            temp_data = load(fullfile(ROI_RSA_dir,['roi_effects_' conditionnames{m} '.mat']));
            all_data(m,:,crun) = temp_data.roi_effect; %Create a matrix of condition by ROI by subject
            RSA_ROI_data_exist(crun) = 1;
        catch
            warning(['No data for ' subjects{crun} ' probably because of SPM dropout, ignoring them'])
            RSA_ROI_data_exist(crun) = 0;
            continue
        end
    end
end
roi_names = temp_data.roi_names;
clear temp_data
disp(['Excluding subjects ' num2str(find(RSA_ROI_data_exist==0)) ' belonging to groups ' num2str(group(RSA_ROI_data_exist==0)) ' maybe check them'])
all_data(:,:,RSA_ROI_data_exist==0) = NaN;

addpath('./plotting')
for m = 1:length(conditionnames)
    figure
    set(gcf,'Position',[100 100 1600 800]);
    set(gcf, 'PaperPositionMode', 'auto');
    hold on
    %violin(squeeze(all_data(m,:,group==1&RSA_ROI_data_exist))','x',[1:length(roi_names)]-0.1)
    errorbar([1:length(roi_names)]-0.1,nanmean(squeeze(all_data(m,:,group==1&RSA_ROI_data_exist)),2),nanstd(squeeze(all_data(m,:,group==1&RSA_ROI_data_exist))')/sqrt(sum(group==1&RSA_ROI_data_exist)),'kx')
    errorbar([1:length(roi_names)]+0.1,nanmean(squeeze(all_data(m,:,group==2&RSA_ROI_data_exist)),2),nanstd(squeeze(all_data(m,:,group==2&RSA_ROI_data_exist))')/sqrt(sum(group==2&RSA_ROI_data_exist)),'rx')
    %violin(squeeze(all_data(m,:,group==2&RSA_ROI_data_exist))','x',[1:length(roi_names)]+0.1)
    xlim([0 length(roi_names)+1])
    set(gca,'xtick',[1:length(roi_names)],'xticklabels',roi_names,'XTickLabelRotation',45,'TickLabelInterpreter','none')
    data_bounds = [min(min(squeeze(all_data(m,:,:)))), max(max(squeeze(all_data(m,:,:))))];
    %ylim([data_bounds(1) - (diff(data_bounds)/10) data_bounds(2) + (diff(data_bounds)/6)])
    plot([0 length(roi_names)+1],[0,0],'k--')
    title(conditionnames{m},'Interpreter','none')
end

%% Analyse in scanner behaviour
graph_individuals = 0;
plot_behavioural_data(subjects, dates, group, graph_individuals) % Requires matlab2016a or newer for movmean



% %% In progress - extract the neural RDM from a given mask
% nrun = size(subjects,2); % enter the number of runs here
% avneuralRDM_clusterworkedcorrectly = zeros(1,nrun);
% downsamp_ratio = 2; %Downsampling in each dimension, must be an integer, 2 is 8 times faster than 1 (2 cubed). 
% 
% extraction_masks = {'rwL_STG_cross-segment_cluster.nii'};
% parfor crun = 1:nrun
%     addpath(genpath('./RSA_scripts'))
%     GLMDir = [preprocessedpathstem subjects{crun} '/stats4_multi_3_nowritten2'];
%     try
%         module_extract_avneuralRDM_cluster(GLMDir,downsamp_ratio,extraction_masks)
%         avneuralRDM_clusterworkedcorrectly(crun) = 1;
%     catch
%         avneuralRDM_clusterworkedcorrectly(crun) = 0;
%     end
% end
% crun = 1;
% GLMDir = [preprocessedpathstem subjects{crun} '/stats4_multi_3_nowritten2'];
% addpath(genpath('./RSA_scripts'))
% module_across_subj_avneuralRDM(GLMDir,downsamp_ratio,extraction_masks,subjects,group);

%% Explore univariate contrast for subject quality check
% On the basis of this and his inability to do the task, P7P14 excluded,
% all other subjects accepted as having at least a large Normal > Written response

this_smooth = 8;
all_conditions = {
        'con_0010.nii','Mismatch > Match';
        'con_0015.nii','Normal > Written';
        'con_0020.nii','Written > Normal';
        'con_0030.nii','Clear > Unclear';
        'con_0035.nii','Unclear > Clear';
        'con_0040.nii','Clarity Congruency Interaction'
        };
explore_univariate_contrast(subjects,preprocessedpathstem,this_smooth,all_conditions)

%% Now run the freesurfer (Must be done after traditional segmentation for masking)
skullstripped = 2; %Use segmentation with imfill to mask the raw image - works best.
module_run_freesurfer;

%% Now do stats on the freesurfer
datafolder = [preprocessedpathstem '/freesurfer_masked/'];
group_names = {'Control','nfvPPA'};
age_lookup = readtable('Pinfa_ages.csv');
studyname = 'PINFA';
precached = 0;
view_data = 0;
smoothing_kernel = 15; %In mm
module_make_FSGD(subjects, group, group_names, age_lookup, datafolder, studyname)
module_run_FSGLM(datafolder, group_names, studyname, precached, view_data, smoothing_kernel)
