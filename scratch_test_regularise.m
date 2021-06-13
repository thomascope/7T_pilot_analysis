nrun = size(subjects,2); % enter the number of runs here

%%First copy targets to a new folder, coregister and crop necks
parfor crun = 1:nrun
    Out_file = [preprocessedpathstem subjects{crun} '/Regularised_structural.nii'];
    if crun>1 % Subject 1 was a pilot and annoyingly did not have INV1 saved
        mkdir([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '_regularised'])
        mkdir([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'INV2'))} '_regularised'])
        % Find INV1 path - series number expected to be 1 or 2 less than INV2
        inv2folder = [rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'INV2'))}];
        Series_number_location = findstr(inv2folder,'Series_')+8;
        Series_number = inv2folder(Series_number_location:Series_number_location+1);
        inv1folder = strrep(inv2folder,'INV2','INV1');
        flag = 2;
        inv1folder = strrep(inv1folder,['Series_0' Series_number],['Series_0' num2str(str2num(Series_number)-flag)]);
        if ~exist(inv1folder,'dir')
            flag = 1;
            inv1folder = strrep(inv2folder,'INV2','INV1');
            inv1folder = strrep(inv1folder,['Series_0' Series_number],['Series_0' num2str(str2num(Series_number)-flag)]);
            if ~exist(inv1folder,'dir')
                error(['Cannot find INV1 folder in the expected location, please check'])
            end
        end
        inv1rawfile = strrep(blocksin{crun}{find(strcmp(blocksout{crun},'INV2'))},[Series_number],[num2str(str2num(Series_number)-flag)]);
        if ~exist([inv1folder '/' inv1rawfile],'file')
            error(['Cannot find INV1 file in the expected location, please check'])
        end
        mkdir([inv1folder '_regularised'])
        copyfile([inv1folder '/' inv1rawfile],[inv1folder '_regularised/' inv1rawfile]);
        copyfile([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '/' blocksin{crun}{find(strcmp(blocksout{crun},'structural'))}],[rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '_regularised/' blocksin{crun}{find(strcmp(blocksout{crun},'structural'))}]);
        copyfile([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'INV2'))} '/' blocksin{crun}{find(strcmp(blocksout{crun},'INV2'))}],[rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'INV2'))} '_regularised/' blocksin{crun}{find(strcmp(blocksout{crun},'INV2'))}]);
        INV1_file = [inv1folder '_regularised/' inv1rawfile];
        UNI_file = [rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '_regularised/' blocksin{crun}{find(strcmp(blocksout{crun},'structural'))}];
        INV2_file = [rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'INV2'))} '_regularised/' blocksin{crun}{find(strcmp(blocksout{crun},'INV2'))}];
        module_regularise_MP2RAGE(UNI_file,INV1_file,INV2_file,Out_file)
        crop_images(Out_file, 1)
    else
        %Can't regularise subject 1 - will need to manually check every segmentation very carefully.
        copyfile([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '/' blocksin{crun}{find(strcmp(blocksout{crun},'structural'))}],Out_file)
        crop_images(Out_file, 1)       
    end
end
Regularised_structurals = {};
for crun = 1:nrun
    Regularised_structurals{crun} = [preprocessedpathstem subjects{crun} '/Regularised_structural.nii'];
end

if visualise_data
    spm_check_registration(char(['/group/language/data/thomascope/spm12_fil_r7771/tpm/TPM.nii';Regularised_structurals']))
end

%% Now do cat12 segmentation
rmpath(genpath('/group/language/data/thomascope/spm12_fil_r6906/'))
addpath /group/language/data/thomascope/spm12_fil_r7771/ % Newset version of cat12
spm fmri

nrun = size(subjects,2); % enter the number of runs here
%jobfile = {'/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/module_cat12_segment_job.m'};
jobfile = {'/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/module_cat12_segment_expert_job.m'};
inputs = cell(1, nrun);
for crun = 1:nrun
    inputs{1, crun} = cellstr(Regularised_structurals{crun}); % CAT12: Segmentation: Volumes - cfg_files
    copyfile(char(inputs{1,crun}),strrep(char(inputs{1,crun}),'.nii','_expert_heavy.nii'))
    inputs{1, crun} = cellstr(strrep(char(inputs{1,crun}),'.nii','_expert_heavy.nii'));
end

cat12segmentworkedcorrectly = zeros(1,nrun);
jobs = repmat(jobfile, 1, 1);

parfor crun = 1:nrun
    rmpath(genpath('/group/language/data/thomascope/spm12_fil_r6906/'))
    addpath /group/language/data/thomascope/spm12_fil_r7771/ % Newset version of cat12
    spm fmri
    spm('defaults', 'fMRI');
    spm_jobman('initcfg')
    cat12('expert')
    try
        spm_jobman('run', jobs, inputs{:,crun});
        cat12segmentworkedcorrectly(crun) = 1;
    catch
        cat12segmentworkedcorrectly(crun) = 0;
    end
end

if ~all(cat12segmentworkedcorrectly)
    error('failed at cat12');
end

% Put old SPM back in action
rmpath(genpath('/group/language/data/thomascope/spm12_fil_r7771/'))
%addpath /imaging/local/software/spm_cbu_svn/releases/spm12_fil_r6906
addpath /group/language/data/thomascope/spm12_fil_r6906/
spm fmri

%Check scans
if visualise_data
    for crun = 1:nrun
        cmd = ['freeview -v ' preprocessedpathstem subjects{crun} '/Regularised_structural.nii ' preprocessedpathstem subjects{crun} '/mri/p1Regularised_structural_expert_heavy.nii:colormap=Heat:opacity=0.2 ' preprocessedpathstem subjects{crun} '/mri/p2Regularised_structural_expert_heavy.nii:colormap=Heat:opacity=0.2 ']
        system(cmd)
    end
end

%% Now calculate TIV

% List of open inputs
% Estimate TIV and global tissue volumes: XML files - cfg_files
nrun = 1; % enter the number of runs here
jobfile = {'/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/module_cat12_tiv_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(2, nrun);
postfix = '_expert_heavy';
all_cat12_xmls = {};
for crun = 1:size(subjects,2)
    all_cat12_xmls(crun) = cellstr(strrep(char(Regularised_structurals{crun}),'Regularised_structural.nii',['/report/cat_Regularised_structural' postfix '.xml']));
end

inputs{1, 1} = cellstr(char(all_cat12_xmls)); % CAT12: Segmentation: Volumes - cfg_files
inputs{2, 1} = [preprocessedpathstem '/cat12VBM/TIV' postfix '.txt']; % Output file
mkdir([preprocessedpathstem '/cat12VBM/TIV'])
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});

%% Now run a VBM
age_lookup = readtable('Pinfa_ages.csv');
tiv_lookup = csvread([preprocessedpathstem '/cat12VBM/TIV' postfix '.txt']);
visual_check = 0;
nrun = size(subjects,2); % enter the number of runs here
%jobfile = {'/group/language/data/thomascope/vespa/SPM12version/Standalone preprocessing pipeline/tc_source/batch_forwardmodel_job_noheadpoints.m'};
group1_mrilist = {}; %NB: Patient MRIs, so here group 2 (sorry)
group1_ages = [];
group1_tivs = [];
group2_mrilist = {};
group2_ages = [];
group2_tivs = [];
this_scan = {};

for crun = 1:nrun
    this_age = age_lookup.Age(strcmp(age_lookup.Study_ID,subjects{crun}));
    this_scan(crun) = cellstr(strrep(char(Regularised_structurals{crun}),'Regularised_structural.nii',['Regularised_structural' postfix '.nii']));
    if group(crun) == 1 % Controls
        group2_mrilist(end+1) = this_scan(crun);
        group2_ages(end+1) = this_age;
        group2_tivs(end+1) = tiv_lookup(crun);
    elseif group(crun) == 2 % Patients
        group1_mrilist(end+1) = this_scan(crun);
        group1_ages(end+1) = this_age;
        group1_tivs(end+1) = tiv_lookup(crun);
    end
end

module_cat12_vbm_job(group1_mrilist, group1_ages, group1_tivs, group2_mrilist, group2_ages, group2_tivs, [preprocessedpathstem '/cat12' postfix],postfix)