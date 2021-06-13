function module_cat12_vbm_job(group1_mrilist, group1_ages, group1_tivs, group2_mrilist, group2_ages, group2_tivs, preprocessfolder, postfix)
% A generic function for taking images from two groups (1=patients, 2=controls) and performing a VBM
% Group 1 should be the same size as or smaller than group 2, and be the
% disease group (i.e. you must have at least as many controls as patients)
% Create a parallel pool before calling the module, preferably with the
% number of workers equalling the number of subjects.
% There is the facility to pre-segment your images, for example you may
% want to use both the uni and INV2 images from an MP2RAGE, which is not
% supported here.

visual_check = 0;

if size(group1_mrilist,1) < size(group1_mrilist,2)
    group1_mrilist = group1_mrilist';
end
if size(group2_mrilist,1) < size(group2_mrilist,2)
    group2_mrilist = group2_mrilist';
end
if size(group1_ages,1) < size(group1_ages,2)
    group1_ages = group1_ages';
end
if size(group2_ages,1) < size(group2_ages,2)
    group2_ages = group2_ages';
end
if size(group1_tivs,1) < size(group1_tivs,2)
    group1_tivs = group1_tivs';
end
if size(group2_tivs,1) < size(group2_tivs,2)
    group2_tivs = group2_tivs';
end

mrilist = [group1_mrilist; group2_mrilist];

split_stem = regexp(mrilist, '/', 'split');

old_imagepaths = cell(length(mrilist),3);
new_imagepaths = cell(length(mrilist),3);
for i = 1:length(mrilist)
    for j = 1:2
        old_imagepaths(i,j) = cellstr(['/' fullfile(split_stem{i}{1:end-1}) '/mri/mwp' num2str(j) split_stem{i}{end}]);
        new_imagepaths(i,j) = cellstr(['/' fullfile(split_stem{i}{1:end-1}) '/mri/mwp' num2str(j) '_ns_' split_stem{i}{end}]);
        try
            movefile(char(old_imagepaths(i,j)),char(new_imagepaths(i,j)))
        catch
            assert(~~exist(['/' fullfile(split_stem{i}{1:end-1}) '/mri/mwp' num2str(j) '_ns_' split_stem{i}{end}]), ['Something went wrong with file ' num2str(i)])
        end
    end
end


%% Then make a DARTEL template based on only the core images (n patients and n controls)
core_imagepaths = [group1_mrilist; group2_mrilist(1:length(group1_mrilist))];

nrun = 1; % enter the number of runs here
jobfile = {'./vbm_scripts/VBM_batch_dartel_namedout_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(2,1);
inputs{1,1} = cell(length(core_imagepaths),1);
inputs{2,1} = cell(length(core_imagepaths),1);
split_stem = regexp(core_imagepaths, '/', 'split');

for i = 1:length(core_imagepaths)
    inputs{1,1}(i) = cellstr(['/' fullfile(split_stem{i}{1:end-1}) '/mri/rp1' split_stem{i}{end}]);
    inputs{2,1}(i) = cellstr(['/' fullfile(split_stem{i}{1:end-1}) '/mri/rp2' split_stem{i}{end}]);
end

inputs{3,1} = ['Template ' postfix];
dartelworkedcorrectly = zeros(1,nrun);
jobs = repmat(jobfile, 1, 1);

parfor crun = 1:nrun %Still submit as a parfor to avoid overloading a login node
    spm('defaults', 'PET');
    spm_jobman('initcfg')
    try
        spm_jobman('run', jobs, inputs{:,crun});
        dartelworkedcorrectly(crun) = 1;
    catch
        dartelworkedcorrectly(crun) = 0;
    end
end

%% Then apply the DARTEL template to the remaining controls

%Find the controls not processed above

[C,ia,ib] = intersect(group2_mrilist,group2_mrilist(1:length(group1_mrilist)),'rows');
total_control_numbers = 1:length(group2_mrilist);
indexes_to_process = setdiff(total_control_numbers,ia);
remaining_control_imagepaths = group2_mrilist(indexes_to_process);

nrun = 1; % enter the number of runs here
jobfile = {'./vbm_scripts/VBM_batch_templated_dartel.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(8, nrun);
inputs{1,1} = cell(length(remaining_control_imagepaths),1);
inputs{2,1} = cell(length(remaining_control_imagepaths),1);

split_stem = regexp(remaining_control_imagepaths, '/', 'split');

for i = 1:length(remaining_control_imagepaths)
    inputs{1,1}(i) = cellstr(['/' fullfile(split_stem{i}{1:end-1}) '/mri/rp1' split_stem{i}{end}]);
    inputs{2,1}(i) = cellstr(['/' fullfile(split_stem{i}{1:end-1}) '/mri/rp2' split_stem{i}{end}]);
end

split_stem = regexp(core_imagepaths, '/', 'split');
inputs{3,1} = cellstr(['/' fullfile(split_stem{1}{1:end-1}) '/mri/Template ' postfix ' _1.nii']);
inputs{4,1} = cellstr(['/' fullfile(split_stem{1}{1:end-1}) '/mri/Template ' postfix ' _2.nii']);
inputs{5,1} = cellstr(['/' fullfile(split_stem{1}{1:end-1}) '/mri/Template ' postfix ' _3.nii']);
inputs{6,1} = cellstr(['/' fullfile(split_stem{1}{1:end-1}) '/mri/Template ' postfix ' _4.nii']);
inputs{7,1} = cellstr(['/' fullfile(split_stem{1}{1:end-1}) '/mri/Template ' postfix ' _5.nii']);
inputs{8,1} = cellstr(['/' fullfile(split_stem{1}{1:end-1}) '/mri/Template ' postfix ' _6.nii']);

path_to_template_6 = cellstr(['/' fullfile(split_stem{1}{1:end-1}) '/mri/Template ' postfix ' _6.nii']);

templateddartelworkedcorrectly = zeros(1,nrun);
jobs = repmat(jobfile, 1, 1);

parfor crun = 1:nrun %Still submit as a parfor to avoid overloading a login node
    spm('defaults', 'PET');
    spm_jobman('initcfg')
    try
        spm_jobman('run', jobs, inputs{:,crun});
        templateddartelworkedcorrectly(crun) = 1;
    catch
        templateddartelworkedcorrectly(crun) = 0;
    end
end

%% Now normalise all scans

nrun = length(mrilist);
jobfile = {'./vbm_scripts/VBM_batch_normalise.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(3, nrun);

split_stem = regexp(mrilist, '/', 'split');
split_stem_template = regexp(core_imagepaths, '/', 'split');
path_to_template_6 = cellstr(['/' fullfile(split_stem_template{1}{1:end-1}) '/mri/Template ' postfix ' _6.nii']);

for crun = 1:nrun
    inputs{1, crun} = path_to_template_6; % Normalise to MNI Space: Dartel Template - cfg_files
    inputs{2, crun} = cellstr(['/' fullfile(split_stem{crun}{1:end-1}) '/mri/u_rp1' split_stem{crun}{end}(1:end-4) '_Template ' postfix '.nii']);
    if ~exist(char(inputs{2,crun}),'file')
        inputs{2, crun} = cellstr(['/' fullfile(split_stem{crun}{1:end-1}) '/mri/u_rp1' split_stem{crun}{end}]); %Deal with the difference in naming convention depending if part of the dartel templating or not
    end
    inputs{3, crun} = cellstr(['/' fullfile(split_stem{crun}{1:end-1}) '/mri/rp1' split_stem{crun}{end}]);
end

normaliseworkedcorrectly = zeros(1,nrun);
jobs = repmat(jobfile, 1, 1);

parfor crun = 1:nrun
    spm('defaults', 'PET');
    spm_jobman('initcfg')
    try
        spm_jobman('run', jobs, inputs{:,crun});
        normaliseworkedcorrectly(crun) = 1;
    catch
        normaliseworkedcorrectly(crun) = 0;
    end
end

%% Now do group stats with TIV and age file as covariates in the ANOVA
nrun = 1;
jobfile = {'./vbm_scripts/VBM_batch_factorial_TIV_age.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(6, nrun);

stats_folder = {[preprocessedpathstem '/cat12VBM/' postfix '/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask']};
split_stem_group2 = regexp(group2_mrilist, '/', 'split');
split_stem_group1 = regexp(group1_mrilist, '/', 'split');

inputs{1, 1} = stats_folder;

for crun = 1:nrun
    inputs{2, 1} = cell(length(group1_mrilist),1);
    for i = 1:length(group1_mrilist)
        inputs{2,crun}(i) = cellstr(['/' fullfile(split_stem_group1{i}{1:end-1}) '/mri/smwp1' split_stem_group1{i}{end}]);
    end
    inputs{3, 1} = cell(length(group2_mrilist),1);
    for i = 1:length(group2_mrilist)
        inputs{3,crun}(i) = cellstr(['/' fullfile(split_stem_group2{i}{1:end-1}) '/mri/smwp1' split_stem_group2{i}{end}]);
    end
end

inputs{4, 1} = [group1_tivs;group2_tivs];

inputs{5, 1} = [group1_ages; group2_ages];

inputs{6, 1} = {'control_majority_unsmoothed_mask_c1_thr0.05_cons0.8.img'};

if ~exist(char(inputs{6, 1}),'file')
    core_imagepaths = [group1_mrilist; group2_mrilist(1:length(group1_mrilist))];
    split_stem_template = regexp(core_imagepaths, '/', 'split');
    path_to_template_6 = cellstr(['/' fullfile(split_stem_template{1}{1:end-1}) '/mri/Template ' postfix ' _6.nii']);
    make_VBM_explicit_mask(group2_mrilist, path_to_template_6, 'control')
end

spm_jobman('run', jobs, inputs{:});

inputs = cell(1, nrun);
inputs{1, 1} =  {[char(stats_folder) '/SPM.mat']};

jobfile = {'./vbm_scripts/VBM_batch_estimate.m'};
jobs = repmat(jobfile, 1, nrun);

spm_jobman('run', jobs, inputs{:});

jobfile = {'./vbm_scripts/VBM_batch_contrast.m'};
jobs = repmat(jobfile, 1, nrun);

spm_jobman('run', jobs, inputs{:});

jobfile = {'./vbm_scripts/VBM_batch_results.m'};
jobs = repmat(jobfile, 1, nrun);

spm_jobman('run', jobs, inputs{:});

%% Make average brain by normalising MPRAGES to template and then averaging them

nrun = length(mrilist);
jobfile = {'./vbm_scripts/VBM_batch_normalise_unmodulated_unsmoothed.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(3, nrun);

core_imagepaths = [group1_mrilist; group2_mrilist(1:length(group1_mrilist))];
split_stem = regexp(mrilist, '/', 'split');
split_stem_template = regexp(core_imagepaths, '/', 'split');
path_to_template_6 = cellstr(['/' fullfile(split_stem_template{1}{1:end-1}) '/Template ' postfix '_6.nii']);



for crun = 1:nrun
    inputs{1, crun} = path_to_template_6; % Normalise to MNI Space: Dartel Template - cfg_files
    inputs{2, crun} = cellstr(['/' fullfile(split_stem{crun}{1:end-1}) '/mri/u_rp1' split_stem{crun}{end}(1:end-4) '_Template ' postfix '.nii']);
    if ~exist(char(inputs{2,crun}),'file')
        inputs{2, crun} = cellstr(['/' fullfile(split_stem{crun}{1:end-1}) '/mri/u_rp1' split_stem{crun}{end}]); %Deal with the difference in naming convention depending if part of the dartel templating or not
    end
    inputs{3, crun} = cellstr(['/' fullfile(split_stem{crun}{1:end-1}) '/' split_stem{crun}{end}]);
end

normaliseforaverageworkedcorrectly = zeros(1,nrun);
jobs = repmat(jobfile, 1, 1);

parfor crun = 1:nrun
    spm('defaults', 'PET');
    spm_jobman('initcfg')
    try
        spm_jobman('run', jobs, inputs{:,crun});
        normaliseforaverageworkedcorrectly(crun) = 1;
    catch
        normaliseforaverageworkedcorrectly(crun) = 0;
    end
end

nrun = 1;
jobfile = {'./vbm_scripts/VBM_batch_imcalc_average.m'};
inputs = cell(2, nrun);

split_stem = regexp(core_imagepaths, '/', 'split');
inputs{1,1} = cell(length(core_imagepaths),1);

for i = 1:length(core_imagepaths)
    inputs{1,1}(i) = cellstr(['/' fullfile(split_stem{i}{1:end-1}) '//w' split_stem{i}{end}]);
end

inputs{2,1} = ['./cat12/average_matched_control_patient_T1head'];

imcalcworkedcorrectly = zeros(1,nrun);
jobs = repmat(jobfile, 1, 1);

parfor crun = 1:nrun %Still submit as a parfor to avoid overloading a login node
    spm('defaults', 'PET');
    spm_jobman('initcfg')
    try
        spm_jobman('run', jobs, inputs{:,crun});
        imcalcworkedcorrectly(crun) = 1;
    catch
        imcalcworkedcorrectly(crun) = 0;
    end
end

%% Now repeat for white matter

%% Now normalise all white matter scans

nrun = length(mrilist);
jobfile = {'./vbm_scripts/VBM_batch_normalise.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(3, nrun);

split_stem = regexp(mrilist, '/', 'split');
split_stem_template = regexp(core_imagepaths, '/', 'split');
path_to_template_6 = cellstr(['/' fullfile(split_stem_template{1}{1:end-1}) '/mri/Template ' postfix ' _6.nii']);

for crun = 1:nrun
    inputs{1, crun} = path_to_template_6; % Normalise to MNI Space: Dartel Template - cfg_files
    inputs{2, crun} = cellstr(['/' fullfile(split_stem{crun}{1:end-1}) '/mri/u_rp2' split_stem{crun}{end}(1:end-4) '_Template ' postfix '.nii']);
    if ~exist(char(inputs{2,crun}),'file')
        inputs{2, crun} = cellstr(['/' fullfile(split_stem{crun}{1:end-1}) '/mri/u_rp2' split_stem{crun}{end}]); %Deal with the difference in naming convention depending if part of the dartel templating or not
    end
    inputs{3, crun} = cellstr(['/' fullfile(split_stem{crun}{1:end-1}) '/mri/rp2' split_stem{crun}{end}]);
end

normaliseworkedcorrectly = zeros(1,nrun);
jobs = repmat(jobfile, 1, 1);

parfor crun = 1:nrun
    spm('defaults', 'PET');
    spm_jobman('initcfg')
    try
        spm_jobman('run', jobs, inputs{:,crun});
        normaliseworkedcorrectly(crun) = 1;
    catch
        normaliseworkedcorrectly(crun) = 0;
    end
end

%% Now do group stats with TIV and age file as covariates in the ANOVA
nrun = 1;
jobfile = {'./vbm_scripts/VBM_batch_factorial_TIV_age.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(6, nrun);

stats_folder = {[preprocessedpathstem '/cat12VBM/' postfix '/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask']};
split_stem_group2 = regexp(group2_mrilist, '/', 'split');
split_stem_group1 = regexp(group1_mrilist, '/', 'split');

inputs{1, 1} = stats_folder;

for crun = 1:nrun
    inputs{2, 1} = cell(length(group1_mrilist),1);
    for i = 1:length(group1_mrilist)
        inputs{2,crun}(i) = cellstr(['/' fullfile(split_stem_group1{i}{1:end-1}) '/mri/smwp2' split_stem_group1{i}{end}]);
    end
    inputs{3, 1} = cell(length(group2_mrilist),1);
    for i = 1:length(group2_mrilist)
        inputs{3,crun}(i) = cellstr(['/' fullfile(split_stem_group2{i}{1:end-1}) '/mri/smwp2' split_stem_group2{i}{end}]);
    end
end

inputs{4, 1} = [group1_tivs;group2_tivs];

inputs{5, 1} = [group1_ages; group2_ages];

inputs{6, 1} = {'control_majority_unsmoothed_mask_c1_thr0.05_cons0.8.img'};

if ~exist(char(inputs{6, 1}),'file')
    core_imagepaths = [group1_mrilist; group2_mrilist(1:length(group1_mrilist))];
    split_stem_template = regexp(core_imagepaths, '/', 'split');
    path_to_template_6 = cellstr(['/' fullfile(split_stem_template{1}{1:end-1}) '/mri/Template ' postfix ' _6.nii']);
    make_VBM_explicit_mask(group2_mrilist, path_to_template_6, 'control')
end

spm_jobman('run', jobs, inputs{:});

inputs = cell(1, nrun);
inputs{1, 1} =  {[char(stats_folder) '/SPM.mat']};

jobfile = {'./vbm_scripts/VBM_batch_estimate.m'};
jobs = repmat(jobfile, 1, nrun);

spm_jobman('run', jobs, inputs{:});

jobfile = {'./vbm_scripts/VBM_batch_contrast.m'};
jobs = repmat(jobfile, 1, nrun);

spm_jobman('run', jobs, inputs{:});

jobfile = {'./vbm_scripts/VBM_batch_results.m'};
jobs = repmat(jobfile, 1, nrun);

spm_jobman('run', jobs, inputs{:});


%% END
disp('VBM done!')