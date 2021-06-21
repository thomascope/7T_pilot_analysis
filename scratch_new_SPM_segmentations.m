
%% First copy targets to a new folder, coregister and crop necks
parfor crun = 1:nrun
    mkdir([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '_newSPMsegment'])
    mkdir([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'INV2'))} '_newSPMsegment'])
    copyfile([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '/' blocksin{crun}{find(strcmp(blocksout{crun},'structural'))}],[rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '_newSPMsegment/' blocksin{crun}{find(strcmp(blocksout{crun},'structural'))}]);
    copyfile([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'INV2'))} '/' blocksin{crun}{find(strcmp(blocksout{crun},'INV2'))}],[rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'INV2'))} '_newSPMsegment/' blocksin{crun}{find(strcmp(blocksout{crun},'INV2'))}]);
%     inputs{1, crun} = cellstr([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '_newSPMsegment/' blocksin{crun}{find(strcmp(blocksout{crun},'structural'))}]);
%     inputs{2, crun} = cellstr([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'INV2'))} '_newSPMsegment/' blocksin{crun}{find(strcmp(blocksout{crun},'INV2'))}]);
%     crop_images([inputs{1, crun};inputs{2, crun}], 2) % Had to type out longhand for parfor
    crop_images([cellstr([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '_newSPMsegment/' blocksin{crun}{find(strcmp(blocksout{crun},'structural'))}]);cellstr([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'INV2'))} '_newSPMsegment/' blocksin{crun}{find(strcmp(blocksout{crun},'INV2'))}])], 2)
        %     %ensure orientations the same
%     vol1 = spm_vol(char(inputs{1, crun}));
%             vol2 = spm_vol(char(inputs{2, crun}));
%             new_vol = vol2;
%             new_vol.mat = vol1.mat;
%             spm_write_vol(new_vol,spm_read_vols(vol2));  
    
end




%% Then skullstrip structural using SPM but with reduced bias cutoff to 30 and bias regularisation to 0.0001
nrun = size(subjects,2); % enter the number of runs here
jobfile = {[scriptdir 'module_skullstrip_INV2_job.m']};
inputs = cell(2, nrun);

for crun = 1:nrun
    %     inputs{1, crun} = cellstr([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '_newSPMsegment/p' blocksin{crun}{find(strcmp(blocksout{crun},'structural'))}]); % Cropped
    %     inputs{2, crun} = cellstr([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'INV2'))} '_newSPMsegment/p' blocksin{crun}{find(strcmp(blocksout{crun},'INV2'))}]); % Cropped
    %Have to use uncropped or dimensions non idential.
    inputs{1, crun} = cellstr([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '_newSPMsegment/' blocksin{crun}{find(strcmp(blocksout{crun},'structural'))}]);
    inputs{2, crun} = cellstr([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'INV2'))} '_newSPMsegment/' blocksin{crun}{find(strcmp(blocksout{crun},'INV2'))}]);
            %     %ensure orientations the same
            vol1 = spm_vol(char(inputs{1, crun}));
            vol2 = spm_vol(char(inputs{2, crun}));
            new_vol = vol2;
            new_vol.mat = vol1.mat;
            spm_write_vol(new_vol,spm_read_vols(vol2));
    
    inputs{3, crun} = cellstr([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '_newSPMsegment/' blocksin{crun}{find(strcmp(blocksout{crun},'structural'))}]);
    inputs{4, crun} = 'structural_newSPMsegment';
    inputs{5, crun} = cellstr([preprocessedpathstem subjects{crun} '/']);
    if ~exist(inputs{5, crun}{1})
        mkdir(inputs{5, crun}{1});
    end
    inputs{6, crun} = 'structural_csf_newSPMsegment';
    inputs{7, crun} = cellstr([preprocessedpathstem subjects{crun} '/']);
    inputs{8, crun} = 'c1structural_newSPMsegment';
    inputs{9, crun} = cellstr([preprocessedpathstem subjects{crun} '/']);
    inputs{10, crun} = 'c2structural_newSPMsegment';
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
    this_scan(crun) = cellstr([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '_newSPMsegment/' blocksin{crun}{find(strcmp(blocksout{crun},'structural'))}]);
    %this_scan(crun) = cellstr([preprocessedpathstem subjects{crun} '/structural.nii']);
    if segmented
        this_segmented(crun) = cellstr([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '_newSPMsegment/c1' blocksin{crun}{find(strcmp(blocksout{crun},'structural'))}]);
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

module_vbm_job(group1_mrilist, group1_ages, group2_mrilist, group2_ages, [preprocessedpathstem '/newSPMsegment'], segmented)