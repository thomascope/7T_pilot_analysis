%% Do an RSA analysis separately if you want (already integrated into previous step for vowels, but now can compare new models etc without repeating the time consuming cross-nobis)
nrun = size(subjects,2); % enter the number of runs here
RSAnobisworkedcorrectly = zeros(1,nrun);
downsamp_ratio = 1; %Downsampling in each dimension, much be an integer, 2 is 8 times faster than 1 (2 cubed). 
parfor crun = 1:nrun
    addpath(genpath('./RSA_scripts'))
    GLMDir = [preprocessedpathstem subjects{crun} '/stats4_multi_3_nowritten2'];
    try
        module_make_essential_effect_maps(GLMDir,downsamp_ratio)
        RSAnobisworkedcorrectly(crun) = 1;
    catch
        RSAnobisworkedcorrectly(crun) = 0;
    end
end

%% Do a partial-correlation based RSA analysis to tease apart the written and spoken word representations
nrun = size(subjects,2); % enter the number of runs here
partialRSAnobisworkedcorrectly = zeros(1,nrun);
downsamp_ratio = 1; %Downsampling in each dimension, much be an integer, 2 is 8 times faster than 1 (2 cubed). 
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
downsamp_ratio = 1; %Downsampling in each dimension, much be an integer, 2 is 8 times faster than 1 (2 cubed). 
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
downsamp_ratio = 1; %Downsampling in each dimension, much be an integer, 2 is 8 times faster than 1 (2 cubed). 
rmpath('/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/RSA_scripts/es_scripts_fMRI') %Stops SPM getting defaults for second level if on path

GLMDir = [preprocessedpathstem subjects{crun} '/stats4_multi_3_nowritten2']; %Template, first subject
outpath = [preprocessedpathstem '/stats4_multi_3_nowritten2/searchlight/downsamp_' num2str(downsamp_ratio) filesep 'second_level']; %Results directory

searchlightsecondlevel = [];
searchlightsecondlevel = module_searchlight_secondlevel(GLMDir,subjects,group,age_lookup,outpath,downsamp_ratio);