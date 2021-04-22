% Makes an explicit mask for input data, smoothed and unsmoothed but unmodulated
function make_VBM_explicit_mask(scan_paths, path_to_template_6, prefix)

nrun = length(scan_paths);
jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_normalise_unmodulated_unsmoothed.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(3, nrun);

split_stem = regexp(scan_paths, '/mprage_125/', 'split');

for crun = 1:nrun
    inputs{1, crun} = path_to_template_6; % Normalise to MNI Space: Dartel Template - cfg_files
    inputs{2, crun} = cellstr([split_stem{crun}{1} '/mprage_125/u_rc1' split_stem{crun}{2}(1:end-4) '_Template.nii']); 
    if ~exist(char(inputs{2,crun}),'file')
        inputs{2, crun} = cellstr([split_stem{crun}{1} '/mprage_125/u_rc1' split_stem{crun}{2}]); %Deal with the difference in naming convention depending if part of the dartel templating or not
    end
    inputs{3, crun} = cellstr([split_stem{crun}{1} '/mprage_125/c2' split_stem{crun}{2}]); 
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

oldfilenames = cell(length(scan_paths),1);
newfilenames = cell(length(scan_paths),1);
for i = 1:length(scan_paths)
    oldfilenames{i}= cellstr([split_stem{i}{1} '/mprage_125/wc2' split_stem{i}{2}]);
    newfilenames{i}= cellstr([split_stem{i}{1} '/mprage_125/s_0_wc2' split_stem{i}{2}]);
    try
        movefile(char(oldfilenames{i}),char(newfilenames{i}));
    catch
    end
end

addpath([pwd '/Masking'])
%  make_majority_mask(thresholds, consensus, outputfilname, files)
for i = 1:length(scan_paths)
    newfilenames{i}= char(newfilenames{i});
end
make_majority_mask([0.2 0.1 0.05 0.001], 0.8, [prefix '_majority_unsmoothed_mask_c2'], char(newfilenames))


nrun = length(scan_paths);
jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/VBM_batch_normalise_unmodulated_smoothed.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(3, nrun);

split_stem = regexp(scan_paths, '/mprage_125/', 'split');

for crun = 1:nrun
    inputs{1, crun} = path_to_template_6; % Normalise to MNI Space: Dartel Template - cfg_files
    inputs{2, crun} = cellstr([split_stem{crun}{1} '/mprage_125/u_rc1' split_stem{crun}{2}(1:end-4) '_Template.nii']); 
    if ~exist(char(inputs{2,crun}),'file')
        inputs{2, crun} = cellstr([split_stem{crun}{1} '/mprage_125/u_rc1' split_stem{crun}{2}]); %Deal with the difference in naming convention depending if part of the dartel templating or not
    end
    inputs{3, crun} = cellstr([split_stem{crun}{1} '/mprage_125/c2' split_stem{crun}{2}]); 
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

oldfilenames = cell(length(scan_paths),1);
newfilenames = cell(length(scan_paths),1);
for i = 1:length(scan_paths)
    oldfilenames{i}= cellstr([split_stem{i}{1} '/mprage_125/swc2' split_stem{i}{2}]);
    newfilenames{i}= cellstr([split_stem{i}{1} '/mprage_125/s_8_wc2' split_stem{i}{2}]);
    try
        movefile(char(oldfilenames{i}),char(newfilenames{i}));
    catch
    end
end

addpath([pwd '/Masking'])
%  make_majority_mask(thresholds, consensus, outputfilname, files)
for i = 1:length(scan_paths)
    newfilenames{i}= char(newfilenames{i});
end
make_majority_mask([0.2 0.1 0.05 0.001], 0.8, [prefix '_majority_smoothed_mask_c2'], char(newfilenames))

