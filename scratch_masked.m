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
    GLMDir = [preprocessedpathstem subjects{crun} '/stats5_multi_3_nowritten2_masked'];
    try
        TDTCrossnobisAnalysis_1Subj(GLMDir,downsamp_ratio)
        mahalanobisworkedcorrectly(crun) = 1;
    catch
        mahalanobisworkedcorrectly(crun) = 0;
    end
end
%% Do an RSA analysis separately if you want (already integrated into previous step for vowels, but now can compare new models etc without repeating the time consuming cross-nobis)
nrun = size(subjects,2); % enter the number of runs here
RSAnobisworkedcorrectly = zeros(1,nrun);
downsamp_ratio = 2; %Downsampling in each dimension, must be an integer, 2 is 8 times faster than 1 (2 cubed). 
parfor crun = 1:nrun
    addpath(genpath('./RSA_scripts'))
    GLMDir = [preprocessedpathstem subjects{crun} '/stats5_multi_3_nowritten2_masked'];
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
    GLMDir = [preprocessedpathstem subjects{crun} '/stats5_multi_3_nowritten2_masked'];
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
    GLMDir = [preprocessedpathstem subjects{crun} '/stats5_multi_3_nowritten2_masked'];
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

GLMDir = [preprocessedpathstem subjects{crun} '/stats5_multi_3_nowritten2_masked']; %Template, first subject
outpath = [preprocessedpathstem '/stats5_multi_3_nowritten2_masked/searchlight/downsamp_' num2str(downsamp_ratio) filesep 'second_level']; %Results directory

searchlightsecondlevel = [];
searchlightsecondlevel = module_searchlight_secondlevel(GLMDir,subjects,group,age_lookup,outpath,downsamp_ratio);
