function module_vbm_job(group1_mrilist, group1_ages, group2_mrilist, group2_ages, preprocessfolder, segmented)
% A generic function for taking images from two groups (1=patients, 2=controls) and performing a VBM
% Group 1 should be the same size as or smaller than group 2, and be the
% disease group (i.e. you must have at least as many controls as patients)
% Create a parallel pool before calling the module, preferably with the
% number of workers equalling the number of subjects.
% There is the facility to pre-segment your images, for example you may
% want to use both the uni and INV2 images from an MP2RAGE, which is not
% supported here.

visual_check = 1

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

mrilist = [group1_mrilist; group2_mrilist];

split_stem = regexp(mrilist, '/', 'split');
if segmented
    for i = 1:size(mrilist,1)
        if ~exist(['/' fullfile(split_stem{i}{1:end-1}) '/c1' split_stem{i}{end}],'file') || ~exist(['/' fullfile(split_stem{i}{1:end-1}) '/c2' split_stem{i}{end}],'file')
            error(['Use existing segmentation specified, but this does not exist for ' fullfile(split_stem{i})])
        end
    end
end

%% First segment all the images
if ~segmented
nrun = length(mrilist);
jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_segment.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(1, nrun);

for crun = 1:nrun
    inputs{1, crun} = cellstr(mrilist{crun}); % for dartel templating
end

segmentworkedcorrectly = zeros(1,nrun);
jobs = repmat(jobfile, 1, 1);

if isempty(gcp('nocreate'))
    cbupool(nrun)
end

parfor crun = 1:nrun
    spm('defaults', 'PET');
    spm_jobman('initcfg')
    try
        spm_jobman('run', jobs, inputs{:,crun});
        segmentworkedcorrectly(crun) = 1;
    catch
        segmentworkedcorrectly(crun) = 0;
    end
end
end
%% Now calculate the TIV and then rename all of the mc files as will be overwritten by DARTEL (if not smoothed)
nrun = 1;
jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_TIV.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(2, nrun);

for crun = 1:nrun
    inputs{1, crun} = cell(length(mrilist),1);
    for i = 1:length(mrilist)
        inputs{1, crun}(i) = cellstr([mrilist{i}(1:end-4) '_seg8.mat']); % for dartel templating
    end
end
mkdir([preprocessfolder filesep ''])
inputs{2,1} = [preprocessfolder filesep 'volumes_VBM.csv'];
tiv_filename = [inputs{2,1}];

TIVworkedcorrectly = zeros(1,nrun);
jobs = repmat(jobfile, 1, 1);

parfor crun = 1:nrun %Even though single shot, submit to parfor to avoid overloading login node
    spm('defaults', 'PET');
    spm_jobman('initcfg')
    try
        spm_jobman('run', jobs, inputs{:,crun});
        TIVworkedcorrectly(crun) = 1;
    catch
        TIVworkedcorrectly(crun) = 0;
    end
end

split_stem = regexp(mrilist, '/', 'split');
old_imagepaths = cell(length(mrilist),3);
new_imagepaths = cell(length(mrilist),3);
for i = 1:length(mrilist)
    for j = 1:3
        old_imagepaths(i,j) = cellstr(['/' fullfile(split_stem{i}{1:end-1}) '/mwc' num2str(j) split_stem{i}{end}]);
        new_imagepaths(i,j) = cellstr(['/' fullfile(split_stem{i}{1:end-1}) '/mwc' num2str(j) '_ns_' split_stem{i}{end}]);
        try
            movefile(char(old_imagepaths(i,j)),char(new_imagepaths(i,j)))
        catch
            assert(~~exist(['/' fullfile(split_stem{i}{1:end-1}) '/mwc' num2str(j) '_ns_' split_stem{i}{end}]), ['Something went wrong with file ' num2str(i)])
        end
    end
end


%% Then make a DARTEL template based on only the core images (n patients and n controls)
core_imagepaths = [group1_mrilist; group2_mrilist(1:length(group1_mrilist))];

nrun = 1; % enter the number of runs here
jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_dartel.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(2,1);
inputs{1,1} = cell(length(core_imagepaths),1);
inputs{2,1} = cell(length(core_imagepaths),1);
split_stem = regexp(core_imagepaths, '/', 'split');

for i = 1:length(core_imagepaths)
    inputs{1,1}(i) = cellstr(['/' fullfile(split_stem{i}{1:end-1}) '/rc1' split_stem{i}{end}]);
    inputs{2,1}(i) = cellstr(['/' fullfile(split_stem{i}{1:end-1}) '/rc2' split_stem{i}{end}]);
end

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
jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_templated_dartel.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(8, nrun);
inputs{1,1} = cell(length(remaining_control_imagepaths),1);
inputs{2,1} = cell(length(remaining_control_imagepaths),1);

split_stem = regexp(remaining_control_imagepaths, '/', 'split');

for i = 1:length(remaining_control_imagepaths)
    inputs{1,1}(i) = cellstr(['/' fullfile(split_stem{i}{1:end-1}) '/rc1' split_stem{i}{end}]);
    inputs{2,1}(i) = cellstr(['/' fullfile(split_stem{i}{1:end-1}) '/rc2' split_stem{i}{end}]);
end

split_stem = regexp(core_imagepaths, '/', 'split');
inputs{3,1} = cellstr(['/' fullfile(split_stem{1}{1:end-1}) '/Template_1.nii']);
inputs{4,1} = cellstr(['/' fullfile(split_stem{1}{1:end-1}) '/Template_2.nii']);
inputs{5,1} = cellstr(['/' fullfile(split_stem{1}{1:end-1}) '/Template_3.nii']);
inputs{6,1} = cellstr(['/' fullfile(split_stem{1}{1:end-1}) '/Template_4.nii']);
inputs{7,1} = cellstr(['/' fullfile(split_stem{1}{1:end-1}) '/Template_5.nii']);
inputs{8,1} = cellstr(['/' fullfile(split_stem{1}{1:end-1}) '/Template_6.nii']);

path_to_template_6 = cellstr(['/' fullfile(split_stem{1}{1:end-1}) '/Template_6.nii']);

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
jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_normalise.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(3, nrun);

split_stem = regexp(mrilist, '/', 'split');
split_stem_template = regexp(core_imagepaths, '/', 'split');
path_to_template_6 = cellstr(['/' fullfile(split_stem_template{1}{1:end-1}) '/Template_6.nii']);

for crun = 1:nrun
    inputs{1, crun} = path_to_template_6; % Normalise to MNI Space: Dartel Template - cfg_files
    inputs{2, crun} = cellstr(['/' fullfile(split_stem{crun}{1:end-1}) '/u_rc1' split_stem{crun}{end}(1:end-4) '_Template.nii']);
    if ~exist(char(inputs{2,crun}),'file')
        inputs{2, crun} = cellstr(['/' fullfile(split_stem{crun}{1:end-1}) '/u_rc1' split_stem{crun}{end}]); %Deal with the difference in naming convention depending if part of the dartel templating or not
    end
    inputs{3, crun} = cellstr(['/' fullfile(split_stem{crun}{1:end-1}) '/c1' split_stem{crun}{end}]);
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

%% Now read in TIV file
if exist('tiv_filename','var')
    filename = tiv_filename;
else
    filename = [preprocessfolder filesep 'volumes_VBM.csv'];
end
delimiter = ',';
startRow = 2;
endRow = inf;
formatSpec = '%s%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
fclose(fileID);
tiv= dataArray{2}+dataArray{3}+dataArray{4};

%% Now do group stats with TIV and age file as covariates in the ANOVA
nrun = 1;
jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_factorial_TIV_age.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(6, nrun);

stats_folder = {[preprocessfolder filesep 'VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask']};
split_stem_group2 = regexp(group2_mrilist, '/', 'split');
split_stem_group1 = regexp(group1_mrilist, '/', 'split');

inputs{1, 1} = stats_folder;

for crun = 1:nrun
    inputs{2, 1} = cell(length(group1_mrilist),1);
    for i = 1:length(group1_mrilist)
        inputs{2,crun}(i) = cellstr(['/' fullfile(split_stem_group1{i}{1:end-1}) '/smwc1' split_stem_group1{i}{end}]);
    end
    inputs{3, 1} = cell(length(group2_mrilist),1);
    for i = 1:length(group2_mrilist)
        inputs{3,crun}(i) = cellstr(['/' fullfile(split_stem_group2{i}{1:end-1}) '/smwc1' split_stem_group2{i}{end}]);
    end
end

try
    inputs{4, 1} = tiv;
catch
    tivstem = [preprocessfolder filesep 'volumes_VBM'];
    filename = [tivstem '.csv'];
    delimiter = ',';
    startRow = 2;
    endRow = inf;
    formatSpec = '%s%f%f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
    fclose(fileID);
    tiv= dataArray{2}+dataArray{3}+dataArray{4};
    inputs{4, 1} = tiv;
end
inputs{5, 1} = [group1_ages; group2_ages];
inputs{6, 1} = {'control_majority_unsmoothed_mask_c1_thr0.05_cons0.8.img'};

if ~exist(char(inputs{6, 1}),'file')
    core_imagepaths = [group1_mrilist; group2_mrilist(1:length(group1_mrilist))];
    split_stem_template = regexp(core_imagepaths, '/', 'split');
    path_to_template_6 = cellstr(['/' fullfile(split_stem_template{1}{1:end-1}) '/Template_6.nii']);
    make_VBM_explicit_mask(group2_mrilist, path_to_template_6, 'control')
end


if visual_check
    all_input_images = [inputs{2};inputs{3}];
    spm_check_registration(all_input_images{:}) 
    input('Press any key to continue')
end

spm_jobman('run', jobs, inputs{:});

inputs = cell(1, nrun);
inputs{1, 1} =  {[char(stats_folder) '/SPM.mat']};

jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_estimate.m'};
jobs = repmat(jobfile, 1, nrun);

spm_jobman('run', jobs, inputs{:});

jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_contrast.m'};
jobs = repmat(jobfile, 1, nrun);

spm_jobman('run', jobs, inputs{:});

jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_results.m'};
jobs = repmat(jobfile, 1, nrun);

spm_jobman('run', jobs, inputs{:});

%% Now do single subject stats with TIV file and age covariates in the ANOVA - for later correlation against behavioural measures
nrun = length(group1_mrilist);
jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_factorial_TIV_age_singlesubj.m'};
jobs = repmat(jobfile, 1, 1);
inputs = cell(6, nrun);

split_stem_group2 = regexp(group2_mrilist, '/', 'split');
split_stem_group1 = regexp(group1_mrilist, '/', 'split');

for crun = 1:nrun
    inputs{1, crun} = {[preprocessfolder filesep 'VBM_stats/factorial_single_subject/patient_' num2str(crun)]};
    inputs{2, crun} = cellstr(['/' fullfile(split_stem_group1{crun}{1:end-1}) '/smwc1' split_stem_group1{crun}{end}]);
    inputs{3, crun} = cell(length(group2_mrilist),1);
    for i = 1:length(group2_mrilist)
        inputs{3,crun}(i) = cellstr(['/' fullfile(split_stem_group2{i}{1:end-1}) '/smwc1' split_stem_group2{i}{end}]);
    end
    try
        inputs{4, crun} = [tiv(crun); tiv(nrun+1:end)];
    catch
        tivstem = [preprocessfolder filesep 'volumes_VBM'];
        filename = [tivstem '.csv'];
        delimiter = ',';
        startRow = 2;
        endRow = inf;
        formatSpec = '%s%f%f%f%[^\n\r]';
        fileID = fopen(filename,'r');
        dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
        fclose(fileID);
        tiv= dataArray{2}+dataArray{3}+dataArray{4};
        inputs{4, crun} = [tiv(crun); tiv(nrun+1:end)];
    end
    inputs{5, crun} = [group1_ages(crun); group2_ages];
    inputs{6, crun} =  {'control_majority_smoothed_mask_c1_thr0.1_cons0.8.img'};
end

if ~exist(char(inputs{6, 1}),'file')
    core_imagepaths = [group1_mrilist; group2_mrilist(1:length(group1_mrilist))];
    split_stem_template = regexp(core_imagepaths, '/', 'split');
    path_to_template_6 = cellstr(['/' fullfile(split_stem_template{1}{1:end-1}) '/Template_6.nii']);
    make_VBM_explicit_mask(group2_mrilist, path_to_template_6, 'control')
end

parfor crun = 1:nrun
    spm('defaults', 'PET');
    spm_jobman('initcfg')
    try
        spm_jobman('run', jobs, inputs{:,crun});
    catch
    end
end

inputs = inputs(1,:);
for i = 1:nrun
    inputs{1,i} = {[char(inputs{1,i}) '/SPM.mat']};
end

jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_estimate.m'};
jobs = repmat(jobfile, 1, 1);

parfor crun = 1:nrun
    spm('defaults', 'PET');
    spm_jobman('initcfg')
    try
        spm_jobman('run', jobs, inputs{:,crun});
    catch
    end
end

jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_contrast.m'};
jobs = repmat(jobfile, 1, 1);

parfor crun = 1:nrun
    spm('defaults', 'PET');
    spm_jobman('initcfg')
    try
        spm_jobman('run', jobs, inputs{:,crun});
    catch
    end
end

jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_results.m'};
jobs = repmat(jobfile, 1, 1);

parfor crun = 1:nrun
    spm('defaults', 'PET');
    spm_jobman('initcfg')
    try
        spm_jobman('run', jobs, inputs{:,crun});
    catch
    end
end

%% Make average brain by normalising MPRAGES to template and then averaging them

nrun = length(mrilist);
jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_normalise_unmodulated_unsmoothed.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(3, nrun);

core_imagepaths = [group1_mrilist; group2_mrilist(1:length(group1_mrilist))];
split_stem = regexp(mrilist, '/', 'split');
split_stem_template = regexp(core_imagepaths, '/', 'split');
path_to_template_6 = cellstr(['/' fullfile(split_stem_template{1}{1:end-1}) '/Template_6.nii']);



for crun = 1:nrun
    inputs{1, crun} = path_to_template_6; % Normalise to MNI Space: Dartel Template - cfg_files
    inputs{2, crun} = cellstr(['/' fullfile(split_stem{crun}{1:end-1}) '/u_rc1' split_stem{crun}{end}(1:end-4) '_Template.nii']);
    if ~exist(char(inputs{2,crun}),'file')
        inputs{2, crun} = cellstr(['/' fullfile(split_stem{crun}{1:end-1}) '/u_rc1' split_stem{crun}{end}]); %Deal with the difference in naming convention depending if part of the dartel templating or not
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
jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_imcalc_average.m'};
inputs = cell(2, nrun);

split_stem = regexp(core_imagepaths, '/', 'split');
inputs{1,1} = cell(length(core_imagepaths),1);

for i = 1:length(core_imagepaths)
    inputs{1,1}(i) = cellstr(['/' fullfile(split_stem{i}{1:end-1}) '//w' split_stem{i}{end}]);
end

inputs{2,1} = ['average_matched_control_patient_T1head'];

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
jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_normalise.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(3, nrun);

split_stem = regexp(mrilist, '/', 'split');
split_stem_template = regexp(core_imagepaths, '/', 'split');
path_to_template_6 = cellstr(['/' fullfile(split_stem_template{1}{1:end-1}) '/Template_6.nii']);

for crun = 1:nrun
    inputs{1, crun} = path_to_template_6; % Normalise to MNI Space: Dartel Template - cfg_files
    inputs{2, crun} = cellstr(['/' fullfile(split_stem{crun}{1:end-1}) '/u_rc1' split_stem{crun}{end}(1:end-4) '_Template.nii']);
    if ~exist(char(inputs{2,crun}),'file')
        inputs{2, crun} = cellstr(['/' fullfile(split_stem{crun}{1:end-1}) '/u_rc1' split_stem{crun}{end}]); %Deal with the difference in naming convention depending if part of the dartel templating or not
    end
    inputs{3, crun} = cellstr(['/' fullfile(split_stem{crun}{1:end-1}) '/c2' split_stem{crun}{end}]);
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


%% Now read in TIV file

if exist('tiv_filename','var')
    filename =tiv_filename;
else
    filename = [preprocessedpathstem filesep 'volumes_VBM.csv'];
end
delimiter = ',';
startRow = 2;
endRow = inf;
formatSpec = '%s%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
fclose(fileID);
tiv= dataArray{2}+dataArray{3}+dataArray{4};

%% Now do white matter group stats with TIV and age file as covariates in the ANOVA
nrun = 1;
jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_factorial_TIV_age.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(6, nrun);

stats_folder = {[preprocessfolder filesep 'VBM_stats/WM_factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask']};
split_stem_group2 = regexp(group2_mrilist, '/', 'split');
split_stem_group1 = regexp(group1_mrilist, '/', 'split');

inputs{1, 1} = stats_folder;

for crun = 1:nrun
    inputs{2, 1} = cell(length(group1_mrilist),1);
    for i = 1:length(group1_mrilist)
        inputs{2,crun}(i) = cellstr(['/' fullfile(split_stem_group1{i}{1:end-1}) '/smwc2' split_stem_group1{i}{end}]);
    end
    inputs{3, 1} = cell(length(group2_mrilist),1);
    for i = 1:length(group2_mrilist)
        inputs{3,crun}(i) = cellstr(['/' fullfile(split_stem_group2{i}{1:end-1}) '/smwc2' split_stem_group2{i}{end}]);
    end
end

try
    inputs{4, 1} = tiv;
catch
    tivstem = [preprocessfolder filesep 'volumes_VBM.csv'];
    filename = [tivstem '.csv'];
    delimiter = ',';
    startRow = 2;
    endRow = inf;
    formatSpec = '%s%f%f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
    fclose(fileID);
    tiv= dataArray{2}+dataArray{3}+dataArray{4};
    inputs{4, 1} = tiv;
end
load('VBM_ages.mat')
inputs{5, 1} = [group1_ages; controlages];
inputs{6, 1} = {'control_majority_unsmoothed_mask_c2_thr0.05_cons0.8.img'};

if ~exist(char(inputs{6, 1}),'file')
    core_imagepaths = [group1_mrilist; group2_mrilist(1:length(group1_mrilist))];
    split_stem_template = regexp(core_imagepaths, '/', 'split');
    path_to_template_6 = cellstr(['/' fullfile(split_stem_template{1}{1:end-1}) '/Template_6.nii']);
    make_WM_VBM_explicit_mask(group2_mrilist, path_to_template_6, 'control')
end

spm_jobman('run', jobs, inputs{:});

inputs = cell(1, nrun);
inputs{1, 1} =  {[char(stats_folder) '/SPM.mat']};

jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_estimate.m'};
jobs = repmat(jobfile, 1, nrun);

spm_jobman('run', jobs, inputs{:});

jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_contrast.m'};
jobs = repmat(jobfile, 1, nrun);

spm_jobman('run', jobs, inputs{:});

jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_results.m'};
jobs = repmat(jobfile, 1, nrun);

spm_jobman('run', jobs, inputs{:});

%% END
disp('VBM done!')