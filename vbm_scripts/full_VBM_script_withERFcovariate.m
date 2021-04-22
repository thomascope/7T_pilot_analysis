% % List of open inputs
% % Segment: Volumes - cfg_files
cd('/imaging/tc02/vespa/scans/PNFA_VBM/tom')
load('image_paths.mat')

% %% First segment all the images
% 
% nrun = length(all_imagepaths);
% jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_segment.m'};
% jobs = repmat(jobfile, 1, nrun);
% inputs = cell(1, nrun);
% 
% for crun = 1:nrun
%     inputs{1, crun} = cellstr(all_imagepaths{crun}); % for dartel templating
% end
% 
% segmentworkedcorrectly = zeros(1,nrun);
% jobs = repmat(jobfile, 1, 1);
% 
% cbupool(nrun)
% 
% parfor crun = 1:nrun
%     spm('defaults', 'PET');
%     spm_jobman('initcfg')
%     try
%         spm_jobman('run', jobs, inputs{:,crun});
%         segmentworkedcorrectly(crun) = 1;
%     catch
%         segmentworkedcorrectly(crun) = 0;
%     end
% end
% 
% % Calculate TIVs
% nrun = 1;
% jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_TIV.m'};
% jobs = repmat(jobfile, 1, nrun);
% inputs = cell(1, nrun);
% 
% for crun = 1:nrun
%     inputs{1, crun} = cell(length(all_imagepaths),1);
%     for i = 1:length(all_imagepaths)
%         %inputs{1, crun} = cellstr([all_imagepaths{crun}(1:end-4) '_seg8.mat'); % for dartel templating
%     end
% end
% 
% TIVworkedcorrectly = zeros(1,nrun);
% jobs = repmat(jobfile, 1, 1);
% 
% parfor crun = 1:nrun %Even though single shot, submit to parfor to avoid overloading login node
%     spm('defaults', 'PET');
%     spm_jobman('initcfg')
%     try
%         spm_jobman('run', jobs, inputs{:,crun});
%         segmentworkedcorrectly(crun) = 1;
%     catch
%         segmentworkedcorrectly(crun) = 0;
%     end
% end
% 
% % Now calculate the TIV and then rename all of the mc files as will be
% % overwritten by DARTEL (if not smoothed)
% nrun = 1;
% jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_TIV.m'};
% jobs = repmat(jobfile, 1, nrun);
% inputs = cell(2, nrun);
% 
% for crun = 1:nrun
%     inputs{1, crun} = cell(length(all_imagepaths),1);
%     for i = 1:length(all_imagepaths)
%         inputs{1, crun}(i) = cellstr([all_imagepaths{i}(1:end-4) '_seg8.mat']); % for dartel templating
%     end
% end
% 
% inputs{2,1} = ['/imaging/tc02/vespa/scans/PNFA_VBM/tom/volumes_PNFA_VBM'];
% tiv_filename = [inputs{2,1} '.csv'];
% 
% TIVworkedcorrectly = zeros(1,nrun);
% jobs = repmat(jobfile, 1, 1);
% 
% parfor crun = 1:nrun %Even though single shot, submit to parfor to avoid overloading login node
%     spm('defaults', 'PET');
%     spm_jobman('initcfg')
%     try
%         spm_jobman('run', jobs, inputs{:,crun});
%         TIVworkedcorrectly(crun) = 1;
%     catch
%         TIVworkedcorrectly(crun) = 0;
%     end
% end
% 
% split_stem = regexp(all_imagepaths, '/mprage_125/', 'split');
% old_imagepaths = cell(length(all_imagepaths),3);
% new_imagepaths = cell(length(all_imagepaths),3);
% for i = 1:length(all_imagepaths)
%     for j = 1:3
%         old_imagepaths(i,j) = cellstr([split_stem{i}{1} '/mprage_125/mwc' num2str(j) split_stem{i}{2}]); 
%         new_imagepaths(i,j) = cellstr([split_stem{i}{1} '/mprage_125/mwc' num2str(j) '_ns_' split_stem{i}{2}]); 
%         movefile(char(old_imagepaths(i,j)),char(new_imagepaths(i,j)))
%     end
% end
% 
% 
% 
% %% Then make a DARTEL template based on only the core images (9 patients and 9 matched controls) and then
% core_imagepaths = [PNFA_full_imagepath; Matched_control_full_imagepath];
% 
% nrun = 1; % enter the number of runs here
% jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_dartel.m'};
% jobs = repmat(jobfile, 1, nrun);
% inputs = cell(2,1);
% inputs{1,1} = cell(length(core_imagepaths),1);
% inputs{2,1} = cell(length(core_imagepaths),1);
% split_stem = regexp(core_imagepaths, '/mprage_125/', 'split');
% 
% for i = 1:length(core_imagepaths)
%    inputs{1,1}(i) = cellstr([split_stem{i}{1} '/mprage_125/rc1' split_stem{i}{2}]); 
%    inputs{2,1}(i) = cellstr([split_stem{i}{1} '/mprage_125/rc2' split_stem{i}{2}]); 
% end
% 
% dartelworkedcorrectly = zeros(1,nrun);
% jobs = repmat(jobfile, 1, 1);
% 
% parfor crun = 1:nrun %Still submit as a parfor to avoid overloading a login node
%     spm('defaults', 'PET');
%     spm_jobman('initcfg')
%     try
%         spm_jobman('run', jobs, inputs{:,crun});
%         dartelworkedcorrectly(crun) = 1;
%     catch
%         dartelworkedcorrectly(crun) = 0;
%     end
% end
% 
% %% Then apply the DARTEL template to the remaining controls
% 
% %Find the controls not processed above
% 
% [C,ia,ib] = intersect(Larger_control_full_imagepath,Matched_control_full_imagepath,'rows');
% total_control_numbers = 1:length(Larger_control_full_imagepath);
% indexes_to_process = setdiff(total_control_numbers,ia);
% remaining_control_imagepaths = Larger_control_full_imagepath(indexes_to_process);
% 
% nrun = 1; % enter the number of runs here
% jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_templated_dartel.m'};
% jobs = repmat(jobfile, 1, nrun);
% inputs = cell(8, nrun);
% inputs{1,1} = cell(length(remaining_control_imagepaths),1);
% inputs{2,1} = cell(length(remaining_control_imagepaths),1);
% 
% split_stem = regexp(remaining_control_imagepaths, '/mprage_125/', 'split');
% for i = 1:length(remaining_control_imagepaths)
%     inputs{1,1}(i) = cellstr([split_stem{i}{1} '/mprage_125/rc1' split_stem{i}{2}]); 
%     inputs{2,1}(i) = cellstr([split_stem{i}{1} '/mprage_125/rc2' split_stem{i}{2}]);
% end
% 
% split_stem = regexp(core_imagepaths, '/mprage_125/', 'split');
% inputs{3,1} = cellstr([split_stem{1}{1} '/mprage_125/Template_1.nii']);
% inputs{4,1} = cellstr([split_stem{1}{1} '/mprage_125/Template_2.nii']);
% inputs{5,1} = cellstr([split_stem{1}{1} '/mprage_125/Template_3.nii']);
% inputs{6,1} = cellstr([split_stem{1}{1} '/mprage_125/Template_4.nii']);
% inputs{7,1} = cellstr([split_stem{1}{1} '/mprage_125/Template_5.nii']);
% inputs{8,1} = cellstr([split_stem{1}{1} '/mprage_125/Template_6.nii']);
% 
% split_stem_template = regexp(core_imagepaths, '/mprage_125/', 'split');
% path_to_template_6 = cellstr([split_stem_template{1}{1} '/mprage_125/Template_6.nii']);
% 
% templateddartelworkedcorrectly = zeros(1,nrun);
% jobs = repmat(jobfile, 1, 1);
% 
% parfor crun = 1:nrun %Still submit as a parfor to avoid overloading a login node
%     spm('defaults', 'PET');
%     spm_jobman('initcfg')
%     try
%         spm_jobman('run', jobs, inputs{:,crun});
%         templateddartelworkedcorrectly(crun) = 1;
%     catch
%         templateddartelworkedcorrectly(crun) = 0;
%     end
% end
% 
% %% Now normalise all scans 
% 
% nrun = length(all_imagepaths);
% jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_normalise.m'};
% jobs = repmat(jobfile, 1, nrun);
% inputs = cell(3, nrun);
% 
% split_stem = regexp(all_imagepaths, '/mprage_125/', 'split');
% split_stem_template = regexp(core_imagepaths, '/mprage_125/', 'split');
% path_to_template_6 = cellstr([split_stem_template{1}{1} '/mprage_125/Template_6.nii']);
% 
% for crun = 1:nrun
%     inputs{1, crun} = path_to_template_6; % Normalise to MNI Space: Dartel Template - cfg_files
%     inputs{2, crun} = cellstr([split_stem{crun}{1} '/mprage_125/u_rc1' split_stem{crun}{2}(1:end-4) '_Template.nii']); 
%     if ~exist(char(inputs{2,crun}),'file')
%         inputs{2, crun} = cellstr([split_stem{crun}{1} '/mprage_125/u_rc1' split_stem{crun}{2}]); %Deal with the difference in naming convention depending if part of the dartel templating or not
%     end
%     inputs{3, crun} = cellstr([split_stem{crun}{1} '/mprage_125/c1' split_stem{crun}{2}]); 
% end
% 
% normaliseworkedcorrectly = zeros(1,nrun);
% jobs = repmat(jobfile, 1, 1);
% 
% parfor crun = 1:nrun
%     spm('defaults', 'PET');
%     spm_jobman('initcfg')
%     try
%         spm_jobman('run', jobs, inputs{:,crun});
%         normaliseworkedcorrectly(crun) = 1;
%     catch
%         normaliseworkedcorrectly(crun) = 0;
%     end
% end
% 
% % %% Now do stats based on grey matter volume normalisation NB: SUBOPTIMAL due to significant loss of grey matter in patient group - better to normalise on TIVs as below.
% % %split_stem_controls = regexp(Larger_control_full_imagepath, '/mprage_125/', 'split');
% % % Problem with normalisation for larger group: Temporary fix:
% % nrun = 1;
% % jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_factorial.m'};
% % jobs = repmat(jobfile, 1, nrun);
% % inputs = cell(3, nrun);
% % 
% % stats_folder = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/Stats/factorial_full_group_vbm_greynormalised'};
% % split_stem_controls = regexp(Larger_control_full_imagepath, '/mprage_125/', 'split');
% % split_stem_PNFA = regexp(PNFA_full_imagepath, '/mprage_125/', 'split');
% % 
% % inputs{1, 1} = stats_folder;
% % 
% % for crun = 1:nrun
% %     inputs{2, 1} = cell(length(PNFA_full_imagepath),1);
% %     for i = 1:length(PNFA_full_imagepath)
% %         inputs{2,crun}(i) = cellstr([split_stem_PNFA{i}{1} '/mprage_125/smwc1' split_stem_PNFA{i}{2}]);
% %     end
% %     inputs{3, 1} = cell(length(Larger_control_full_imagepath),1);
% %     for i = 1:length(Larger_control_full_imagepath)
% %         inputs{3,crun}(i) = cellstr([split_stem_controls{i}{1} '/mprage_125/smwc1' split_stem_controls{i}{2}]);
% %     end
% % end
% 
% % spm_jobman('run', jobs, inputs{:});
% % 
% % inputs = cell(1, nrun);
% % inputs{1, 1} =  {[char(stats_folder) '/SPM.mat']};
% % 
% % jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_estimate.m'};
% % jobs = repmat(jobfile, 1, nrun);
% % 
% % spm_jobman('run', jobs, inputs{:});
% % 
% % jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_contrast.m'};
% % jobs = repmat(jobfile, 1, nrun);
% % 
% % spm_jobman('run', jobs, inputs{:});
% % 
% % jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_results.m'};
% % jobs = repmat(jobfile, 1, nrun);
% % 
% % spm_jobman('run', jobs, inputs{:});
% 
% %% Now read in TIV file
% if exist(tiv_filename,'var')
%     filename =tiv_filename;
% else
%     filename = ['/imaging/tc02/vespa/scans/PNFA_VBM/tom/volumes_PNFA_VBM.csv'];
% end
% delimiter = ',';
% startRow = 2;
% endRow = inf;
% formatSpec = '%s%f%f%f%[^\n\r]';
% fileID = fopen(filename,'r');
% dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
% fclose(fileID);
% tiv= dataArray{2}+dataArray{3}+dataArray{4};
% 
%% Now do group stats with TIV and age file as covariates in the ANOVA
nrun = 1;
jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_factorial_TIV_age_ERF_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(5, nrun);

stats_folder = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/Stats/ERF_factorial_full_group_vbm_TIVnormalised_agecovaried_smoothedmask'};
split_stem_controls = regexp(Larger_control_full_imagepath, '/mprage_125/', 'split');
split_stem_PNFA = regexp(PNFA_full_imagepath, '/mprage_125/', 'split');

inputs{1, 1} = stats_folder;

for crun = 1:nrun
    inputs{2, 1} = cell(length(PNFA_full_imagepath),1);
    for i = 1:length(PNFA_full_imagepath)
        inputs{2,crun}(i) = cellstr([split_stem_PNFA{i}{1} '/mprage_125/smwc1' split_stem_PNFA{i}{2}]);
    end
end

try
    inputs{3, 1} = tiv(1:9);
catch
    tivstem = ['/imaging/tc02/vespa/scans/PNFA_VBM/tom/volumes_PNFA_VBM'];
    filename = [tivstem '.csv'];
    delimiter = ',';
    startRow = 2;
    endRow = inf;
    formatSpec = '%s%f%f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
    fclose(fileID);
    tiv= dataArray{2}+dataArray{3}+dataArray{4};
    inputs{3, 1} = tiv;
end
load('VBM_ages.mat')
inputs{4, 1} = [pnfaages];

PNFA_ERFs = [1048	624	688	584	412	416	516	608 732]';
Control_ERF_mean = 348;

inputs{5, 1} = [PNFA_ERFs];

inputs{6, 1} = {'control_majority_smoothed_mask_c1_thr0.05_cons0.8.img'};

if ~exist(char(inputs{6, 1}),'file')
    core_imagepaths = [PNFA_full_imagepath; Matched_control_full_imagepath];
    split_stem_template = regexp(core_imagepaths, '/mprage_125/', 'split');
    path_to_template_6 = cellstr([split_stem_template{1}{1} '/mprage_125/Template_6.nii']);
    make_VBM_explicit_mask(Larger_control_full_imagepath, path_to_template_6, 'control')
end

spm_jobman('run', jobs, inputs{:});

inputs = cell(1, nrun);
inputs{1, 1} =  {[char(stats_folder) '/SPM.mat']};

jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_estimate.m'};
jobs = repmat(jobfile, 1, nrun);

spm_jobman('run', jobs, inputs{:});

jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_contrast_ERF.m'};
jobs = repmat(jobfile, 1, nrun);

spm_jobman('run', jobs, inputs{:});

jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_results.m'};
jobs = repmat(jobfile, 1, nrun);

spm_jobman('run', jobs, inputs{:});

% %% Now do group stats with TIV file as covariates in the ANOVA
% nrun = 1;
% jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_factorial_TIV_ERF_job.m'};
% jobs = repmat(jobfile, 1, nrun);
% inputs = cell(5, nrun);
% 
% stats_folder = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/Stats/ERF_factorial_full_group_vbm_TIVnormalised_unsmoothedmask'};
% split_stem_controls = regexp(Larger_control_full_imagepath, '/mprage_125/', 'split');
% split_stem_PNFA = regexp(PNFA_full_imagepath, '/mprage_125/', 'split');
% 
% inputs{1, 1} = stats_folder;
% 
% for crun = 1:nrun
%     inputs{2, 1} = cell(length(PNFA_full_imagepath),1);
%     for i = 1:length(PNFA_full_imagepath)
%         inputs{2,crun}(i) = cellstr([split_stem_PNFA{i}{1} '/mprage_125/smwc1' split_stem_PNFA{i}{2}]);
%     end
% end
% 
% try
%     inputs{3, 1} = tiv(1:9);
% catch
%     tivstem = ['/imaging/tc02/vespa/scans/PNFA_VBM/tom/volumes_PNFA_VBM'];
%     filename = [tivstem '.csv'];
%     delimiter = ',';
%     startRow = 2;
%     endRow = inf;
%     formatSpec = '%s%f%f%f%[^\n\r]';
%     fileID = fopen(filename,'r');
%     dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
%     fclose(fileID);
%     tiv= dataArray{2}+dataArray{3}+dataArray{4};
%     inputs{3, 1} = tiv;
% end
% 
% PNFA_ERFs = [1048	624	688	584	412	416	516	608 732]';
% Control_ERF_mean = 348;
% 
% inputs{4, 1} = [PNFA_ERFs];
% 
% inputs{5, 1} = {'control_majority_unsmoothed_mask_c1_thr0.05_cons0.8.img'};
% 
% if ~exist(char(inputs{5, 1}),'file')
%     core_imagepaths = [PNFA_full_imagepath; Matched_control_full_imagepath];
%     split_stem_template = regexp(core_imagepaths, '/mprage_125/', 'split');
%     path_to_template_6 = cellstr([split_stem_template{1}{1} '/mprage_125/Template_6.nii']);
%     make_VBM_explicit_mask(Larger_control_full_imagepath, path_to_template_6, 'control')
% end
% 
% spm_jobman('run', jobs, inputs{:});
% 
% inputs = cell(1, nrun);
% inputs{1, 1} =  {[char(stats_folder) '/SPM.mat']};
% 
% jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_estimate.m'};
% jobs = repmat(jobfile, 1, nrun);
% 
% spm_jobman('run', jobs, inputs{:});
% 
% jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_contrast_ERF_noage.m'};
% jobs = repmat(jobfile, 1, nrun);
% 
% spm_jobman('run', jobs, inputs{:});
% 
% jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_results.m'};
% jobs = repmat(jobfile, 1, nrun);
% 
% spm_jobman('run', jobs, inputs{:});
% 



%% Now repeat for white matter

nrun = 1;
jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_factorial_TIV_age_ERF_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(5, nrun);

stats_folder = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/Stats/ERF_WM_factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask'};
split_stem_controls = regexp(Larger_control_full_imagepath, '/mprage_125/', 'split');
split_stem_PNFA = regexp(PNFA_full_imagepath, '/mprage_125/', 'split');

inputs{1, 1} = stats_folder;

for crun = 1:nrun
    inputs{2, 1} = cell(length(PNFA_full_imagepath),1);
    for i = 1:length(PNFA_full_imagepath)
        inputs{2,crun}(i) = cellstr([split_stem_PNFA{i}{1} '/mprage_125/smwc2' split_stem_PNFA{i}{2}]);
    end
end

try
    inputs{3, 1} = tiv(1:9);
catch
    tivstem = ['/imaging/tc02/vespa/scans/PNFA_VBM/tom/volumes_PNFA_VBM'];
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
inputs{4, 1} = [pnfaages];

PNFA_ERFs = [1048	624	688	584	412	416	516	608 732]';
Control_ERF_mean = 348;

inputs{5, 1} = [PNFA_ERFs];

inputs{6, 1} = {'control_majority_unsmoothed_mask_c2_thr0.05_cons0.8.img'};

if ~exist(char(inputs{6, 1}),'file')
    core_imagepaths = [PNFA_full_imagepath; Matched_control_full_imagepath];
    split_stem_template = regexp(core_imagepaths, '/mprage_125/', 'split');
    path_to_template_6 = cellstr([split_stem_template{1}{1} '/mprage_125/Template_6.nii']);
    make_VBM_explicit_mask(Larger_control_full_imagepath, path_to_template_6, 'control')
end

spm_jobman('run', jobs, inputs{:});

inputs = cell(1, nrun);
inputs{1, 1} =  {[char(stats_folder) '/SPM.mat']};

jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_estimate.m'};
jobs = repmat(jobfile, 1, nrun);

spm_jobman('run', jobs, inputs{:});

jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_contrast_ERF.m'};
jobs = repmat(jobfile, 1, nrun);

spm_jobman('run', jobs, inputs{:});

jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_results.m'};
jobs = repmat(jobfile, 1, nrun);

spm_jobman('run', jobs, inputs{:});
% 
% %% Now do group stats against just the subjects in the dartel with TIV file and age as covariates in the ANOVA
% nrun = 1;
% jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_factorial_TIV_age.m'};
% jobs = repmat(jobfile, 1, nrun);
% inputs = cell(6, nrun);
% 
% stats_folder = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/Stats/factorial_matched_group_vbm_TIVnormalised_agecovaried_unsmoothedmask'};
% split_stem_controls = regexp(Matched_control_full_imagepath, '/mprage_125/', 'split');
% split_stem_PNFA = regexp(PNFA_full_imagepath, '/mprage_125/', 'split');
% 
% [C,ia,ib] = intersect(Larger_control_full_imagepath,Matched_control_full_imagepath,'rows');
% 
% inputs{1, 1} = stats_folder;
% 
% for crun = 1:nrun
%     inputs{2, 1} = cell(length(PNFA_full_imagepath),1);
%     for i = 1:length(PNFA_full_imagepath)
%         inputs{2,crun}(i) = cellstr([split_stem_PNFA{i}{1} '/mprage_125/smwc1' split_stem_PNFA{i}{2}]);
%     end
%     inputs{3, 1} = cell(length(Matched_control_full_imagepath),1);
%     for i = 1:length(Matched_control_full_imagepath)
%         inputs{3,crun}(i) = cellstr([split_stem_controls{i}{1} '/mprage_125/smwc1' split_stem_controls{i}{2}]);
%     end
% end
% 
% try
%     inputs{4, 1} = [tiv(1:length(PNFA_full_imagepath)); tiv(ia)];
% catch
%     tivstem = ['/imaging/tc02/vespa/scans/PNFA_VBM/tom/volumes_PNFA_VBM'];
%     filename = [tivstem '.csv'];
%     delimiter = ',';
%     startRow = 2;
%     endRow = inf;
%     formatSpec = '%s%f%f%f%[^\n\r]';
%     fileID = fopen(filename,'r');
%     dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
%     fclose(fileID);
%     tiv= dataArray{2}+dataArray{3}+dataArray{4};
%     inputs{4, 1} = [tiv(1:length(PNFA_full_imagepath)); tiv(ia)];
% end
% load('VBM_ages.mat')
% inputs{5, 1} = [pnfaages; controlages(ia)];
% inputs{6, 1} = {'matched_control_majority_unsmoothed_mask_c1_thr0.05_cons0.8.img'};
% 
% if ~exist(char(inputs{6, 1}),'file')
%     core_imagepaths = [PNFA_full_imagepath; Matched_control_full_imagepath];
%     split_stem_template = regexp(core_imagepaths, '/mprage_125/', 'split');
%     path_to_template_6 = cellstr([split_stem_template{1}{1} '/mprage_125/Template_6.nii']);
%     make_VBM_explicit_mask(Matched_control_full_imagepath, path_to_template_6, 'matched_control')
% end
% 
% spm_jobman('run', jobs, inputs{:});
% 
% inputs = cell(1, nrun);
% inputs{1, 1} =  {[char(stats_folder) '/SPM.mat']};
% 
% jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_estimate.m'};
% jobs = repmat(jobfile, 1, nrun);
% 
% spm_jobman('run', jobs, inputs{:});
% 
% jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_contrast.m'};
% jobs = repmat(jobfile, 1, nrun);
% 
% spm_jobman('run', jobs, inputs{:});
% 
% jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_results.m'};
% jobs = repmat(jobfile, 1, nrun);
% 
% spm_jobman('run', jobs, inputs{:});
% 
% 
% %% Now do single subject stats with TIV file and age covariates in the ANOVA
% nrun = length(PNFA_full_imagepath);
% jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_factorial_TIV_age_singlesubj.m'};
% jobs = repmat(jobfile, 1, 1);
% inputs = cell(6, nrun);
% 
% split_stem_controls = regexp(Larger_control_full_imagepath, '/mprage_125/', 'split');
% split_stem_PNFA = regexp(PNFA_full_imagepath, '/mprage_125/', 'split');
% load('VBM_ages.mat')
% for crun = 1:nrun
%     inputs{1, crun} = {['/imaging/tc02/vespa/scans/PNFA_VBM/tom/Stats/factorial_single_subject/' split_stem_PNFA{crun}{2}(1:end-4)]};
%     inputs{2, crun} = cellstr([split_stem_PNFA{crun}{1} '/mprage_125/smwc1' split_stem_PNFA{crun}{2}]);
%     inputs{3, crun} = cell(length(Larger_control_full_imagepath),1);
%     for i = 1:length(Larger_control_full_imagepath)
%         inputs{3,crun}(i) = cellstr([split_stem_controls{i}{1} '/mprage_125/smwc1' split_stem_controls{i}{2}]);
%     end
%     try
%         inputs{4, crun} = [tiv(crun); tiv(nrun+1:end)];
%     catch
%         tivstem = ['/imaging/tc02/vespa/scans/PNFA_VBM/tom/volumes_PNFA_VBM'];
%         filename = [tivstem '.csv'];
%         delimiter = ',';
%         startRow = 2;
%         endRow = inf;
%         formatSpec = '%s%f%f%f%[^\n\r]';
%         fileID = fopen(filename,'r');
%         dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
%         fclose(fileID);
%         tiv= dataArray{2}+dataArray{3}+dataArray{4};
%         inputs{4, crun} = [tiv(crun); tiv(nrun+1:end)];
%     end
%     inputs{5, crun} = [pnfaages(crun); controlages];
%     inputs{6, crun} =  {'control_majority_smoothed_mask_c1_thr0.1_cons0.8.img'};
% end
% 
% if ~exist(char(inputs{6, 1}),'file')
%     core_imagepaths = [PNFA_full_imagepath; Matched_control_full_imagepath];
%     split_stem_template = regexp(core_imagepaths, '/mprage_125/', 'split');
%     path_to_template_6 = cellstr([split_stem_template{1}{1} '/mprage_125/Template_6.nii']);
%     make_VBM_explicit_mask(Larger_control_full_imagepath, path_to_template_6, 'control')
% end
% 
% parfor crun = 1:nrun
%     spm('defaults', 'PET');
%     spm_jobman('initcfg')
%     try
%         spm_jobman('run', jobs, inputs{:,crun});
%     catch
%     end
% end
% 
% inputs = inputs(1,:);
% for i = 1:nrun
%     inputs{1,i} = {[char(inputs{1,i}) '/SPM.mat']};
% end
% 
% jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_estimate.m'};
% jobs = repmat(jobfile, 1, 1);
% 
% parfor crun = 1:nrun
%     spm('defaults', 'PET');
%     spm_jobman('initcfg')
%     try
%         spm_jobman('run', jobs, inputs{:,crun});
%     catch
%     end
% end
% 
% jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_contrast.m'};
% jobs = repmat(jobfile, 1, 1);
% 
% parfor crun = 1:nrun
%     spm('defaults', 'PET');
%     spm_jobman('initcfg')
%     try
%         spm_jobman('run', jobs, inputs{:,crun});
%     catch
%     end
% end
% 
% jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_results.m'};
% jobs = repmat(jobfile, 1, 1);
% 
% parfor crun = 1:nrun
%     spm('defaults', 'PET');
%     spm_jobman('initcfg')
%     try
%         spm_jobman('run', jobs, inputs{:,crun});
%     catch
%     end
% end
% 
% %% Repeat single subject stats against just the controls in the dartel with TIV file and age covariates in the ANOVA
% nrun = length(PNFA_full_imagepath);
% jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_factorial_TIV_age_singlesubj.m'};
% jobs = repmat(jobfile, 1, 1);
% inputs = cell(6, nrun);
% 
% [C,ia,ib] = intersect(Larger_control_full_imagepath,Matched_control_full_imagepath,'rows');
% 
% split_stem_controls = regexp(Matched_control_full_imagepath, '/mprage_125/', 'split');
% split_stem_PNFA = regexp(PNFA_full_imagepath, '/mprage_125/', 'split');
% load('VBM_ages.mat')
% for crun = 1:nrun
%     inputs{1, crun} = {['/imaging/tc02/vespa/scans/PNFA_VBM/tom/Stats/factorial_single_subject_against_matched_controls/' split_stem_PNFA{crun}{2}(1:end-4)]};
%     inputs{2, crun} = cellstr([split_stem_PNFA{crun}{1} '/mprage_125/smwc1' split_stem_PNFA{crun}{2}]);
%     inputs{3, crun} = cell(length(Matched_control_full_imagepath),1);
%     for i = 1:length(Matched_control_full_imagepath)
%         inputs{3,crun}(i) = cellstr([split_stem_controls{i}{1} '/mprage_125/smwc1' split_stem_controls{i}{2}]);
%     end
%     try
%         inputs{4, crun} = [tiv(crun); tiv(ia)];
%     catch
%         tivstem = ['/imaging/tc02/vespa/scans/PNFA_VBM/tom/volumes_PNFA_VBM'];
%         filename = [tivstem '.csv'];
%         delimiter = ',';
%         startRow = 2;
%         endRow = inf;
%         formatSpec = '%s%f%f%f%[^\n\r]';
%         fileID = fopen(filename,'r');
%         dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
%         fclose(fileID);
%         tiv= dataArray{2}+dataArray{3}+dataArray{4};
%         inputs{4, crun} = [tiv(crun); tiv(nrun+1:end)];
%     end
%     inputs{5, crun} = [pnfaages(crun); controlages(ia)];
%     inputs{6, crun} =  {'matched_control_majority_smoothed_mask_c1_thr0.1_cons0.8.img'};
% end
% 
% if ~exist(char(inputs{6, 1}),'file')
%     core_imagepaths = [PNFA_full_imagepath; Matched_control_full_imagepath];
%     split_stem_template = regexp(core_imagepaths, '/mprage_125/', 'split');
%     path_to_template_6 = cellstr([split_stem_template{1}{1} '/mprage_125/Template_6.nii']);
%     make_VBM_explicit_mask(Matched_control_full_imagepath, path_to_template_6, 'matched_control')
% end
% 
% parfor crun = 1:nrun
%     spm('defaults', 'PET');
%     spm_jobman('initcfg')
%     try
%         spm_jobman('run', jobs, inputs{:,crun});
%     catch
%     end
% end
% 
% inputs = inputs(1,:);
% for i = 1:nrun
%     inputs{1,i} = {[char(inputs{1,i}) '/SPM.mat']};
% end
% 
% jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_estimate.m'};
% jobs = repmat(jobfile, 1, 1);
% 
% parfor crun = 1:nrun
%     spm('defaults', 'PET');
%     spm_jobman('initcfg')
%     try
%         spm_jobman('run', jobs, inputs{:,crun});
%     catch
%     end
% end
% 
% jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_contrast.m'};
% jobs = repmat(jobfile, 1, 1);
% 
% parfor crun = 1:nrun
%     spm('defaults', 'PET');
%     spm_jobman('initcfg')
%     try
%         spm_jobman('run', jobs, inputs{:,crun});
%     catch
%     end
% end
% 
% jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_results.m'};
% jobs = repmat(jobfile, 1, 1);
% 
% parfor crun = 1:nrun
%     spm('defaults', 'PET');
%     spm_jobman('initcfg')
%     try
%         spm_jobman('run', jobs, inputs{:,crun});
%     catch
%     end
% end
% 
% 
% %% View single subjects results cell - does not work as a standalone
% for crun = 1:nrun
%         spm_jobman('run', jobs, inputs{:,crun});
%         
%         input('press enter for next spm')
%         
% end
% 
% %% Make average brain by normalising MPRAGES to template and then averaging them
% 
% nrun = length(all_imagepaths);
% jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_normalise_unmodulated_unsmoothed.m'};
% jobs = repmat(jobfile, 1, nrun);
% inputs = cell(3, nrun);
% 
% split_stem = regexp(all_imagepaths, '/mprage_125/', 'split');
% split_stem_template = regexp(core_imagepaths, '/mprage_125/', 'split');
% path_to_template_6 = cellstr([split_stem_template{1}{1} '/mprage_125/Template_6.nii']);
% 
% for crun = 1:nrun
%     inputs{1, crun} = path_to_template_6; % Normalise to MNI Space: Dartel Template - cfg_files
%     inputs{2, crun} = cellstr([split_stem{crun}{1} '/mprage_125/u_rc1' split_stem{crun}{2}(1:end-4) '_Template.nii']); 
%     if ~exist(char(inputs{2,crun}),'file')
%         inputs{2, crun} = cellstr([split_stem{crun}{1} '/mprage_125/u_rc1' split_stem{crun}{2}]); %Deal with the difference in naming convention depending if part of the dartel templating or not
%     end
%     inputs{3, crun} = cellstr([split_stem{crun}{1} '/mprage_125/' split_stem{crun}{2}]); 
% end
% 
% normaliseforaverageworkedcorrectly = zeros(1,nrun);
% jobs = repmat(jobfile, 1, 1);
% 
% parfor crun = 1:nrun
%     spm('defaults', 'PET');
%     spm_jobman('initcfg')
%     try
%         spm_jobman('run', jobs, inputs{:,crun});
%         normaliseforaverageworkedcorrectly(crun) = 1;
%     catch
%         normaliseforaverageworkedcorrectly(crun) = 0;
%     end
% end
% 
% nrun = 1;
% jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_imcalc_average.m'};
% inputs = cell(2, nrun);
% 
% split_stem = regexp(core_imagepaths, '/mprage_125/', 'split');
% inputs{1,1} = cell(length(core_imagepaths),1);
% 
% for i = 1:length(core_imagepaths)
%    inputs{1,1}(i) = cellstr([split_stem{i}{1} '/mprage_125/w' split_stem{i}{2}]); 
% end
% 
% inputs{2,1} = ['average_9control_9PNFA_T1head'];
% 
% imcalcworkedcorrectly = zeros(1,nrun);
% jobs = repmat(jobfile, 1, 1);
% 
% parfor crun = 1:nrun %Still submit as a parfor to avoid overloading a login node
%     spm('defaults', 'PET');
%     spm_jobman('initcfg')
%     try
%         spm_jobman('run', jobs, inputs{:,crun});
%         imcalcworkedcorrectly(crun) = 1;
%     catch
%         imcalcworkedcorrectly(crun) = 0;
%     end
% end
% 
% %% Now repeat for white matter
% 
% %% Now normalise all white matter scans 
% 
% nrun = length(all_imagepaths);
% jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_normalise.m'};
% jobs = repmat(jobfile, 1, nrun);
% inputs = cell(3, nrun);
% 
% split_stem = regexp(all_imagepaths, '/mprage_125/', 'split');
% split_stem_template = regexp(core_imagepaths, '/mprage_125/', 'split');
% path_to_template_6 = cellstr([split_stem_template{1}{1} '/mprage_125/Template_6.nii']);
% 
% for crun = 1:nrun
%     inputs{1, crun} = path_to_template_6; % Normalise to MNI Space: Dartel Template - cfg_files
%     inputs{2, crun} = cellstr([split_stem{crun}{1} '/mprage_125/u_rc1' split_stem{crun}{2}(1:end-4) '_Template.nii']); 
%     if ~exist(char(inputs{2,crun}),'file')
%         inputs{2, crun} = cellstr([split_stem{crun}{1} '/mprage_125/u_rc1' split_stem{crun}{2}]); %Deal with the difference in naming convention depending if part of the dartel templating or not
%     end
%     inputs{3, crun} = cellstr([split_stem{crun}{1} '/mprage_125/c2' split_stem{crun}{2}]); 
% end
% 
% normaliseworkedcorrectly = zeros(1,nrun);
% jobs = repmat(jobfile, 1, 1);
% 
% parfor crun = 1:nrun
%     spm('defaults', 'PET');
%     spm_jobman('initcfg')
%     try
%         spm_jobman('run', jobs, inputs{:,crun});
%         normaliseworkedcorrectly(crun) = 1;
%     catch
%         normaliseworkedcorrectly(crun) = 0;
%     end
% end
% 
% 
% %% Now read in TIV file
% 
% if exist('tiv_filename','var')
%     filename =tiv_filename;
% else
%     filename = ['/imaging/tc02/vespa/scans/PNFA_VBM/tom/volumes_PNFA_VBM.csv'];
% end
% delimiter = ',';
% startRow = 2;
% endRow = inf;
% formatSpec = '%s%f%f%f%[^\n\r]';
% fileID = fopen(filename,'r');
% dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
% fclose(fileID);
% tiv= dataArray{2}+dataArray{3}+dataArray{4};
% 
% %% Now do white matter group stats with TIV and age file as covariates in the ANOVA
% nrun = 1;
% jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_factorial_TIV_age.m'};
% jobs = repmat(jobfile, 1, nrun);
% inputs = cell(6, nrun);
% 
% stats_folder = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/Stats/WM_factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask'};
% split_stem_controls = regexp(Larger_control_full_imagepath, '/mprage_125/', 'split');
% split_stem_PNFA = regexp(PNFA_full_imagepath, '/mprage_125/', 'split');
% 
% inputs{1, 1} = stats_folder;
% 
% for crun = 1:nrun
%     inputs{2, 1} = cell(length(PNFA_full_imagepath),1);
%     for i = 1:length(PNFA_full_imagepath)
%         inputs{2,crun}(i) = cellstr([split_stem_PNFA{i}{1} '/mprage_125/smwc2' split_stem_PNFA{i}{2}]);
%     end
%     inputs{3, 1} = cell(length(Larger_control_full_imagepath),1);
%     for i = 1:length(Larger_control_full_imagepath)
%         inputs{3,crun}(i) = cellstr([split_stem_controls{i}{1} '/mprage_125/smwc2' split_stem_controls{i}{2}]);
%     end
% end
% 
% try
%     inputs{4, 1} = tiv;
% catch
%     tivstem = ['/imaging/tc02/vespa/scans/PNFA_VBM/tom/volumes_PNFA_VBM'];
%     filename = [tivstem '.csv'];
%     delimiter = ',';
%     startRow = 2;
%     endRow = inf;
%     formatSpec = '%s%f%f%f%[^\n\r]';
%     fileID = fopen(filename,'r');
%     dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
%     fclose(fileID);
%     tiv= dataArray{2}+dataArray{3}+dataArray{4};
%     inputs{4, 1} = tiv;
% end
% load('VBM_ages.mat')
% inputs{5, 1} = [pnfaages; controlages];
% inputs{6, 1} = {'control_majority_unsmoothed_mask_c2_thr0.05_cons0.8.img'};
% 
% if ~exist(char(inputs{6, 1}),'file')
%     core_imagepaths = [PNFA_full_imagepath; Matched_control_full_imagepath];
%     split_stem_template = regexp(core_imagepaths, '/mprage_125/', 'split');
%     path_to_template_6 = cellstr([split_stem_template{1}{1} '/mprage_125/Template_6.nii']);
%     make_WM_VBM_explicit_mask(Larger_control_full_imagepath, path_to_template_6, 'control')
% end
% 
% spm_jobman('run', jobs, inputs{:});
% 
% inputs = cell(1, nrun);
% inputs{1, 1} =  {[char(stats_folder) '/SPM.mat']};
% 
% jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_estimate.m'};
% jobs = repmat(jobfile, 1, nrun);
% 
% spm_jobman('run', jobs, inputs{:});
% 
% jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_contrast.m'};
% jobs = repmat(jobfile, 1, nrun);
% 
% spm_jobman('run', jobs, inputs{:});
% 
% jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_results.m'};
% jobs = repmat(jobfile, 1, nrun);
% 
% spm_jobman('run', jobs, inputs{:});
% 
% %% Now do GLM within PNFA group with metrics as covariates. NB: MUST ACCOUNT FOR DIFFERENT ORDERING IN THE PNFAMETRICS FILE AND SCAN LIST
% % 
% PNFA_ERFs = [1048	624	688	584	412	416	516	608 732]';
% 
% split_stem_PNFA = regexp(PNFA_full_imagepath, '/mprage_125/', 'split');
% 
% LOCA = zeros(10,1);
% LOCB = zeros(9,1);
% all_numbers = zeros(9,1);
% for i = 1:9
% [~,LOCB(i)] = ismember(str2num(split_stem_PNFA{i}{1}(end-4:end)),(PNFAmetrics.PNFA_WBIC_NUMBERS));
% all_numbers(i) = str2num(split_stem_PNFA{i}{1}(end-4:end));
% end
% for i = 1:10
% [~,LOCA(i)] = ismember((PNFAmetrics.PNFA_WBIC_NUMBERS(i)),all_numbers);
% end
% 
% %PNFAmetrics.PNFA_WBIC_NUMBERS(LOCB)
% covariate_types = {
%     'ERF latency'
%     };
% 
% tissue_types = {
%     'GM'
%     'WM'
%     };
% 
% mask_types = {
%     'WholeBrain'
%     'AtrophicROIs'
%     };
% 
% start_directory = pwd;
% 
% for covariate_number = 1:4
%     for tissue_number = 1:2
%         for mask_number = 1:2
%             cd(start_directory)
%             nrun = 1;
%             jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_TIV_regression.m'};
%             jobs = repmat(jobfile, 1, nrun);
%             
%             inputs = cell(6, nrun);
%             stats_folder = {['/imaging/tc02/vespa/scans/PNFA_VBM/tom/Stats/regression/' mask_types{mask_number} '_' tissue_types{tissue_number} '_' covariate_types{covariate_number}]};
%             split_stem_controls = regexp(Larger_control_full_imagepath, '/mprage_125/', 'split');
%             
%             inputs{1, 1} = stats_folder;
%             
%             if tissue_number == 1
%                 for crun = 1:nrun
%                     inputs{2, 1} = cell(length(PNFA_full_imagepath),1);
%                     for i = 1:length(PNFA_full_imagepath)
%                         inputs{2,crun}(i) = cellstr([split_stem_PNFA{i}{1} '/mprage_125/smwc1' split_stem_PNFA{i}{2}]);
%                     end
%                 end
%             elseif tissue_number == 2
%                 for crun = 1:nrun
%                     inputs{2, 1} = cell(length(PNFA_full_imagepath),1);
%                     for i = 1:length(PNFA_full_imagepath)
%                         inputs{2,crun}(i) = cellstr([split_stem_PNFA{i}{1} '/mprage_125/smwc2' split_stem_PNFA{i}{2}]);
%                     end
%                 end
%             end
%             
%             unorderedcovariate = eval(['PNFAmetrics.' covariate_types{covariate_number}]);
%             inputs{3,crun} = unorderedcovariate(LOCB);
%             
%             inputs{4,crun} = [covariate_types{covariate_number}];
%             
%             try
%                 inputs{5, 1} = tiv(1:length(PNFA_full_imagepath));
%             catch
%                 tivstem = ['/imaging/tc02/vespa/scans/PNFA_VBM/tom/volumes_PNFA_VBM'];
%                 filename = [tivstem '.csv'];
%                 delimiter = ',';
%                 startRow = 2;
%                 endRow = inf;
%                 formatSpec = '%s%f%f%f%[^\n\r]';
%                 fileID = fopen(filename,'r');
%                 dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
%                 fclose(fileID);
%                 tiv= dataArray{2}+dataArray{3}+dataArray{4};
%                 inputs{5, 1} = tiv(1:length(PNFA_full_imagepath));
%             end
%             load('VBM_ages.mat')
%             inputs{6, 1} = pnfaages;
%             if mask_number == 1
%                 if tissue_number == 1
%                     inputs{7, 1} = {'control_majority_unsmoothed_mask_c1_thr0.05_cons0.8.img'};
%                 elseif tissue_number == 2
%                     inputs{7, 1} = {'control_majority_unsmoothed_mask_c2_thr0.05_cons0.8.img'};
%                 end
%             else
%                 if tissue_number == 1
%                     inputs{7, 1} = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/Stats/VBM_illustration_files/GM_clusterwise_PNFA.nii'};
%                 elseif tissue_number == 2
%                     inputs{7, 1} = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/Stats/VBM_illustration_files/WM_clusterwise_PNFA.nii'};
%                 end
%             end
%             
%             if ~exist(char(inputs{7, 1}),'file')
%                 core_imagepaths = [PNFA_full_imagepath; Matched_control_full_imagepath];
%                 split_stem_template = regexp(core_imagepaths, '/mprage_125/', 'split');
%                 path_to_template_6 = cellstr([split_stem_template{1}{1} '/mprage_125/Template_6.nii']);
%                 make_WM_VBM_explicit_mask(Larger_control_full_imagepath, path_to_template_6, 'control')
%             end
%             
%             spm_jobman('run', jobs, inputs{:});
%             
%             inputs = cell(1, nrun);
%             inputs{1, 1} =  {[char(stats_folder) '/SPM.mat']};
%             
%             jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_estimate.m'};
%             jobs = repmat(jobfile, 1, nrun);
%             
%             spm_jobman('run', jobs, inputs{:});
%             
%             jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_regression_contrast.m'};
%             jobs = repmat(jobfile, 1, nrun);
%             
%             spm_jobman('run', jobs, inputs{:});
%             
%             jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_results.m'};
%             jobs = repmat(jobfile, 1, nrun);
%             
%             spm_jobman('run', jobs, inputs{:});
%         end
%     end
% end
% %% Now extract ROI values for GM - WILL NOT WORK WITHOUT PRESTEPS, SEE COMMENTS
% 
% % Need to normalise atlas to template - change files in /home/tc02/Regionofinterestscripts/dartel_mni_icbm/warp_ICBM_atlas_to_template.m
% 
% % Then normalise to MNI - /home/tc02/Regionofinterestscripts/dartel_mni_icbm/dartel_roi_to_mni.m
% 
% % aal regional lookup table: /home/tc02/PickAtlas/WFU_PickAtlas_3.0.5b/wfu_pickatlas/MNI_atlas_templates
% 
% split_stem_PNFA = regexp(PNFA_full_imagepath, '/mprage_125/', 'split');
% filepath = cell(length(PNFA_full_imagepath),1);
% for crun = 1:length(PNFA_full_imagepath)
%     filepath{crun} = [split_stem_PNFA{crun}{1} '/mprage_125/smwc1' split_stem_PNFA{crun}{2}];
% end
% 
% mask =  {'matched_control_majority_smoothed_mask_c1_thr0.1_cons0.8.img'};
% 
% %spm_summarise(char(filepath), '/home/tc02/Regionofinterestscripts/rwaal_nfvPPA_2mni.nii')
% addpath('/home/tc02/Regionofinterestscripts/masaroi/')
% masaroi(0,Inf, ['a'], char(filepath), ['/home/tc02/Regionofinterestscripts/rwaal_nfvPPA_2mni.nii'],char(mask),[],1,'aal') ;
% 
% filename = [pwd '/nfvPPA_roi_density.txt'];
% movefile([pwd '/all_roi_results.txt'],filename)
% 
% roidata = importdata('nfvPPA_roi_density.txt');
% 
% try
%     all_tivs = tiv;
% catch
%     tivstem = ['/imaging/tc02/vespa/scans/PNFA_VBM/tom/volumes_PNFA_VBM'];
%     filename = [tivstem '.csv'];
%     delimiter = ',';
%     startRow = 2;
%     endRow = inf;
%     formatSpec = '%s%f%f%f%[^\n\r]';
%     fileID = fopen(filename,'r');
%     dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
%     fclose(fileID);
%     tiv= dataArray{2}+dataArray{3}+dataArray{4};
%     all_tivs = tiv;
% end
% 
% PNFAmetrics = load('PNFAmetrics_NBORDERINGNOTSAMEASSCANS');
% LOCA = zeros(10,1);
% LOCB = zeros(9,1);
% all_numbers = zeros(9,1);
% for i = 1:9
% [~,LOCB(i)] = ismember(str2num(split_stem_PNFA{i}{1}(end-4:end)),(PNFAmetrics.PNFA_WBIC_NUMBERS));
% all_numbers(i) = str2num(split_stem_PNFA{i}{1}(end-4:end));
% end
% for i = 1:10
% [~,LOCA(i)] = ismember((PNFAmetrics.PNFA_WBIC_NUMBERS(i)),all_numbers);
% end
% 
% all_ROInames = importdata('/home/tc02/PickAtlas/WFU_PickAtlas_3.0.5b/wfu_pickatlas/MNI_atlas_templates/aal_MNI_V4.txt');
% 
% 
% ROInames = {'Frontal_Inf_Tri_L', 'Frontal_Inf_Oper_L', 'Rolandic_Oper_L', 'Putamen_L', 'Caudate_L', 'Frontal_Inf_Tri_R', 'Frontal_Inf_Oper_R', 'Rolandic_Oper_R', 'Putamen_R', 'Caudate_R'};
% ROInumbers = zeros(size(ROInames));
% for i = 1:length(ROInames)
%     thisROI_place = strfind(all_ROInames,ROInames{i});
%     thisROI_index = find(not(cellfun('isempty', thisROI_place)));
%     thisROI_number = all_ROInames{thisROI_index}(1:(thisROI_place{thisROI_index}-1));
%     ROInumbers(i) = str2num(thisROI_number);
% end
% 
% %ROInumber is in column 1, mean is column 7, median is column 8, nonzero mean is column 14, masked mean is column 28 and median is column 29
% extracted_roi_values = zeros(10,length(ROInames));
% extracted_masked_roi_values = zeros(10,length(ROInames));
% extracted_nonzero_roi_values = zeros(10,length(ROInames));
% correctedgrey = zeros(10,1);
% for i = 1:10
%     if LOCA(i) ~= 0 
%         numbertoworkon = all_numbers(LOCA(i));
%         thistiv = all_tivs(LOCA(i)); %Remember to correct for TIV!
%         thistotalgrey = dataArray{2}(LOCA(i));
%     else
%         numbertoworkon = NaN;
%         thistiv = NaN;
%         thistotalgrey = NaN;
%     end
%     
%     IndexC = strfind(roidata.textdata(:,1),num2str(numbertoworkon));
%     Index_thispatient = find(not(cellfun('isempty', IndexC)));
%     Index_thispatient = Index_thispatient-1; %account for header row
%     
%     for j = 1:length(ROInames)
%         IndexC = roidata.data(:,1) == ROInumbers(j); 
%         current_row = Index_thispatient(IndexC(Index_thispatient));
%         if isempty(current_row)
%             extracted_roi_values(i,j) = NaN;
%             extracted_nonzero_roi_values(i,j) = NaN;
%             extracted_masked_roi_values(i,j) = NaN;
%         else
%             extracted_roi_values(i,j) = roidata.data(current_row,7)/thistiv;
%             extracted_nonzero_roi_values(i,j) = roidata.data(current_row,14)/thistiv;
%             extracted_masked_roi_values(i,j) = roidata.data(current_row,28)/thistiv;
%         end
%     end
%     
%     correctedgrey(i) = thistotalgrey/thistiv;
%   
% end
%     
% 
% %%
% matlabpool close