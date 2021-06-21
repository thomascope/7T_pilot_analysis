visualise_data = 0;
addpath ./crop
nrun = size(subjects,2); % enter the number of runs here

%% First copy scans to a new folder and regularise the UNI image
Regularised_structurals = {};
parfor crun = 1:nrun
    Regularised_structurals{crun} = [rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '_reg_SPM00001/reg' blocksin{crun}{find(strcmp(blocksout{crun},'structural'))}];
    if crun>1 % Subject 1 was a pilot and annoyingly did not have INV1 saved
        mkdir([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '_reg_SPM00001'])
        mkdir([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'INV2'))} '_reg_SPM00001'])
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
        mkdir([inv1folder '_reg_SPM00001'])
        copyfile([inv1folder '/' inv1rawfile],[inv1folder '_reg_SPM00001/' inv1rawfile]);
        copyfile([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '/' blocksin{crun}{find(strcmp(blocksout{crun},'structural'))}],[rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '_reg_SPM00001/' blocksin{crun}{find(strcmp(blocksout{crun},'structural'))}]);
        copyfile([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'INV2'))} '/' blocksin{crun}{find(strcmp(blocksout{crun},'INV2'))}],[rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'INV2'))} '_reg_SPM00001/' blocksin{crun}{find(strcmp(blocksout{crun},'INV2'))}]);
        INV1_file = [inv1folder '_reg_SPM00001/' inv1rawfile];
        UNI_file = [rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '_reg_SPM00001/' blocksin{crun}{find(strcmp(blocksout{crun},'structural'))}];
        INV2_file = [rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'INV2'))} '_reg_SPM00001/' blocksin{crun}{find(strcmp(blocksout{crun},'INV2'))}];
        module_regularise_MP2RAGE(UNI_file,INV1_file,INV2_file,Regularised_structurals{crun})
        crop_images([cellstr([Regularised_structurals{crun}]);cellstr([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'INV2'))} '_reg_SPM00001/' blocksin{crun}{find(strcmp(blocksout{crun},'INV2'))}])],2);
    else
        %Can't regularise subject 1 - will need to manually check every segmentation very carefully.
        mkdir([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '_reg_SPM00001'])
        mkdir([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'INV2'))} '_reg_SPM00001'])
        copyfile([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '/' blocksin{crun}{find(strcmp(blocksout{crun},'structural'))}],Regularised_structurals{crun});
        copyfile([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'INV2'))} '/' blocksin{crun}{find(strcmp(blocksout{crun},'INV2'))}],[rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'INV2'))} '_reg_SPM00001/' blocksin{crun}{find(strcmp(blocksout{crun},'INV2'))}]);
        crop_images([cellstr([Regularised_structurals{crun}]);cellstr([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'INV2'))} '_reg_SPM00001/' blocksin{crun}{find(strcmp(blocksout{crun},'INV2'))}])],2);
    end
end

if visualise_data
    spm_check_registration(char(['/group/language/data/thomascope/spm12_fil_r7771/tpm/TPM.nii';Regularised_structurals']))
end


%% Then skullstrip structural using SPM but with reduced bias cutoff to 30 and bias regularisation to 0.0001
nrun = size(subjects,2); % enter the number of runs here
jobfile = {[scriptdir 'module_skullstrip_INV2_job.m']};
inputs = cell(2, nrun);

for crun = 1:nrun
    %     inputs{1, crun} = cellstr([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '_reg_SPM00001/p' blocksin{crun}{find(strcmp(blocksout{crun},'structural'))}]); % Cropped
    %     inputs{2, crun} = cellstr([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'INV2'))} '_reg_SPM00001/p' blocksin{crun}{find(strcmp(blocksout{crun},'INV2'))}]); % Cropped
    %Have to use uncropped or dimensions non idential.
    inputs{1, crun} = cellstr([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '_reg_SPM00001/reg' blocksin{crun}{find(strcmp(blocksout{crun},'structural'))}]);
    inputs{2, crun} = cellstr([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'INV2'))} '_reg_SPM00001/' blocksin{crun}{find(strcmp(blocksout{crun},'INV2'))}]);
    %     %ensure orientations the same
    vol1 = spm_vol(char(inputs{1, crun}));
    vol2 = spm_vol(char(inputs{2, crun}));
    new_vol = vol2;
    new_vol.mat = vol1.mat;
    spm_write_vol(new_vol,spm_read_vols(vol2));
    
    inputs{3, crun} = cellstr([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '_reg_SPM00001/reg' blocksin{crun}{find(strcmp(blocksout{crun},'structural'))}]);
    inputs{4, crun} = 'structural_reg_SPM00001';
    inputs{5, crun} = cellstr([preprocessedpathstem subjects{crun} '/']);
    if ~exist(inputs{5, crun}{1})
        mkdir(inputs{5, crun}{1});
    end
    inputs{6, crun} = 'structural_csf_reg_SPM00001';
    inputs{7, crun} = cellstr([preprocessedpathstem subjects{crun} '/']);
    inputs{8, crun} = 'c1structural_reg_SPM00001';
    inputs{9, crun} = cellstr([preprocessedpathstem subjects{crun} '/']);
    inputs{10, crun} = 'c2structural_reg_SPM00001';
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

if visualise_data
    for crun = 1:nrun
        this_scan(crun) = cellstr([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '_reg_SPM00001/reg' blocksin{crun}{find(strcmp(blocksout{crun},'structural'))}]);
        this_segmented(crun) = cellstr([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '_reg_SPM00001/c1reg' blocksin{crun}{find(strcmp(blocksout{crun},'structural'))}]);
        cmd = ['freeview -v ' this_scan{crun} ' ' this_segmented{crun} ':heatscale=0.1,1:colormap=Heat:opacity=0.2 ']
        system(cmd)
    end
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
    this_scan(crun) = cellstr([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '_reg_SPM00001/reg' blocksin{crun}{find(strcmp(blocksout{crun},'structural'))}]);
    %this_scan(crun) = cellstr([preprocessedpathstem subjects{crun} '/structural.nii']);
    if segmented
        this_segmented(crun) = cellstr([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '_reg_SPM00001/c1reg' blocksin{crun}{find(strcmp(blocksout{crun},'structural'))}]);
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

module_vbm_job(group1_mrilist, group1_ages, group2_mrilist, group2_ages, [preprocessedpathstem '/regSPM_00001'], segmented)