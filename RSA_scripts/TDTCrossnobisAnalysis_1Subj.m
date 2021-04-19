function []=TDTCrossnobisAnalysis_1Subj(GLMDir)

addpath('/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/RSA_scripts/es_scripts_fMRI')
addpath('/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/RSA_scripts/decoding_toolbox_v3.999')
addpath(genpath('/group/language/data/ediz.sohoglu/matlab/rsatoolbox'));

%GLMDir = fullfile(GLMAnalPD,SubjID);
%GLMDir = '/imaging/mlr/users/tc02/PINFA_preprocessed_2021/P7C05/stats4_multi_3_nowritten2'; %For testing

%StrDir = fullfile(ProcDataDir,SubjID,'Structural');

%% Compute searchlight crossnobis distance

% This script is a template that can be used for an encoding analysis on 
% brain image data using cross-validated Mahalanobis distance. It is for 
% people who have betas available from an SPM.mat (for AFNI, see
% decoding_tutorial) and want to automatically extract the relevant images
% used for calculation of the cross-validated Mahalanobis distance, as well
% as corresponding labels and decoding chunk numbers (e.g. run numbers). If
% you don't have this available, then inspect the differences between
% decoding_template and decoding_template_nobetas and adapt this template
% to use it without betas.

% Set defaults
cfg = decoding_defaults;

% Set the analysis that should be performed (default is 'searchlight')
cfg.analysis = 'searchlight';

% Set the output directory where data will be saved, e.g. 'c:\exp\results\buttonpress'
cfg.results.dir = fullfile(GLMDir,'TDTcrossnobis');

% Set the filepath where your SPM.mat and all related betas are, e.g. 'c:\exp\glm\model_button'
beta_loc = GLMDir;

% Set the filename of your brain mask (or your ROI masks as cell matrix)
% for searchlight or wholebrain e.g. 'c:\exp\glm\model_button\mask.img' OR
% for ROI e.g. {'c:\exp\roi\roimaskleft.img', 'c:\exp\roi\roimaskright.img'}
% You can also use a mask file with multiple masks inside that are
% separated by different integer values (a "multi-mask")
cfg.files.mask = fullfile(GLMDir,'mask.nii');

% Set the label names to the regressor names which you want to use for
% your similarity analysis, e.g.
%labelnames = {'Strong+M_Set1_Item1','Strong+M_Set1_Item2','Strong+M_Set1_Item3','Strong+M_Set1_Item4','Strong+M_Set2_Item1','Strong+M_Set2_Item2','Strong+M_Set2_Item3','Strong+M_Set2_Item4','Strong+M_Set3_Item1','Strong+M_Set3_Item2','Strong+M_Set3_Item3','Strong+M_Set3_Item4','Strong+M_Set4_Item1','Strong+M_Set4_Item2','Strong+M_Set4_Item3','Strong+M_Set4_Item4','Strong+M_Set5_Item1','Strong+M_Set5_Item2','Strong+M_Set5_Item3','Strong+M_Set5_Item4','Strong+M_Set6_Item1','Strong+M_Set6_Item2','Strong+M_Set6_Item3','Strong+M_Set6_Item4','Strong+M_Set7_Item1','Strong+M_Set7_Item2','Strong+M_Set7_Item3','Strong+M_Set7_Item4','Strong+M_Set8_Item1','Strong+M_Set8_Item2','Strong+M_Set8_Item3','Strong+M_Set8_Item4','Weak+M_Set1_Item1','Weak+M_Set1_Item2','Weak+M_Set1_Item3','Weak+M_Set1_Item4','Weak+M_Set2_Item1','Weak+M_Set2_Item2','Weak+M_Set2_Item3','Weak+M_Set2_Item4','Weak+M_Set3_Item1','Weak+M_Set3_Item2','Weak+M_Set3_Item3','Weak+M_Set3_Item4','Weak+M_Set4_Item1','Weak+M_Set4_Item2','Weak+M_Set4_Item3','Weak+M_Set4_Item4','Weak+M_Set5_Item1','Weak+M_Set5_Item2','Weak+M_Set5_Item3','Weak+M_Set5_Item4','Weak+M_Set6_Item1','Weak+M_Set6_Item2','Weak+M_Set6_Item3','Weak+M_Set6_Item4','Weak+M_Set7_Item1','Weak+M_Set7_Item2','Weak+M_Set7_Item3','Weak+M_Set7_Item4','Weak+M_Set8_Item1','Weak+M_Set8_Item2','Weak+M_Set8_Item3','Weak+M_Set8_Item4','Strong+MM_Set1_Item1','Strong+MM_Set1_Item2','Strong+MM_Set1_Item3','Strong+MM_Set1_Item4','Strong+MM_Set2_Item1','Strong+MM_Set2_Item2','Strong+MM_Set2_Item3','Strong+MM_Set2_Item4','Strong+MM_Set3_Item1','Strong+MM_Set3_Item2','Strong+MM_Set3_Item3','Strong+MM_Set3_Item4','Strong+MM_Set4_Item1','Strong+MM_Set4_Item2','Strong+MM_Set4_Item3','Strong+MM_Set4_Item4','Strong+MM_Set5_Item1','Strong+MM_Set5_Item2','Strong+MM_Set5_Item3','Strong+MM_Set5_Item4','Strong+MM_Set6_Item1','Strong+MM_Set6_Item2','Strong+MM_Set6_Item3','Strong+MM_Set6_Item4','Strong+MM_Set7_Item1','Strong+MM_Set7_Item2','Strong+MM_Set7_Item3','Strong+MM_Set7_Item4','Strong+MM_Set8_Item1','Strong+MM_Set8_Item2','Strong+MM_Set8_Item3','Strong+MM_Set8_Item4','Weak+MM_Set1_Item1','Weak+MM_Set1_Item2','Weak+MM_Set1_Item3','Weak+MM_Set1_Item4','Weak+MM_Set2_Item1','Weak+MM_Set2_Item2','Weak+MM_Set2_Item3','Weak+MM_Set2_Item4','Weak+MM_Set3_Item1','Weak+MM_Set3_Item2','Weak+MM_Set3_Item3','Weak+MM_Set3_Item4','Weak+MM_Set4_Item1','Weak+MM_Set4_Item2','Weak+MM_Set4_Item3','Weak+MM_Set4_Item4','Weak+MM_Set5_Item1','Weak+MM_Set5_Item2','Weak+MM_Set5_Item3','Weak+MM_Set5_Item4','Weak+MM_Set6_Item1','Weak+MM_Set6_Item2','Weak+MM_Set6_Item3','Weak+MM_Set6_Item4','Weak+MM_Set7_Item1','Weak+MM_Set7_Item2','Weak+MM_Set7_Item3','Weak+MM_Set7_Item4','Weak+MM_Set8_Item1','Weak+MM_Set8_Item2','Weak+MM_Set8_Item3','Weak+MM_Set8_Item4','Strong+Noise_Set1_Item1','Strong+Noise_Set1_Item2','Strong+Noise_Set1_Item3','Strong+Noise_Set1_Item4','Strong+Noise_Set2_Item1','Strong+Noise_Set2_Item2','Strong+Noise_Set2_Item3','Strong+Noise_Set2_Item4','Strong+Noise_Set3_Item1','Strong+Noise_Set3_Item2','Strong+Noise_Set3_Item3','Strong+Noise_Set3_Item4','Strong+Noise_Set4_Item1','Strong+Noise_Set4_Item2','Strong+Noise_Set4_Item3','Strong+Noise_Set4_Item4','Strong+Noise_Set5_Item1','Strong+Noise_Set5_Item2','Strong+Noise_Set5_Item3','Strong+Noise_Set5_Item4','Strong+Noise_Set6_Item1','Strong+Noise_Set6_Item2','Strong+Noise_Set6_Item3','Strong+Noise_Set6_Item4','Strong+Noise_Set7_Item1','Strong+Noise_Set7_Item2','Strong+Noise_Set7_Item3','Strong+Noise_Set7_Item4','Strong+Noise_Set8_Item1','Strong+Noise_Set8_Item2','Strong+Noise_Set8_Item3','Strong+Noise_Set8_Item4','Weak+Noise_Set1_Item1','Weak+Noise_Set1_Item2','Weak+Noise_Set1_Item3','Weak+Noise_Set1_Item4','Weak+Noise_Set2_Item1','Weak+Noise_Set2_Item2','Weak+Noise_Set2_Item3','Weak+Noise_Set2_Item4','Weak+Noise_Set3_Item1','Weak+Noise_Set3_Item2','Weak+Noise_Set3_Item3','Weak+Noise_Set3_Item4','Weak+Noise_Set4_Item1','Weak+Noise_Set4_Item2','Weak+Noise_Set4_Item3','Weak+Noise_Set4_Item4','Weak+Noise_Set5_Item1','Weak+Noise_Set5_Item2','Weak+Noise_Set5_Item3','Weak+Noise_Set5_Item4','Weak+Noise_Set6_Item1','Weak+Noise_Set6_Item2','Weak+Noise_Set6_Item3','Weak+Noise_Set6_Item4','Weak+Noise_Set7_Item1','Weak+Noise_Set7_Item2','Weak+Noise_Set7_Item3','Weak+Noise_Set7_Item4','Weak+Noise_Set8_Item1','Weak+Noise_Set8_Item2','Weak+Noise_Set8_Item3','Weak+Noise_Set8_Item4','Noise+Speech_Set1_Item1','Noise+Speech_Set1_Item2','Noise+Speech_Set1_Item3','Noise+Speech_Set1_Item4','Noise+Speech_Set2_Item1','Noise+Speech_Set2_Item2','Noise+Speech_Set2_Item3','Noise+Speech_Set2_Item4','Noise+Speech_Set3_Item1','Noise+Speech_Set3_Item2','Noise+Speech_Set3_Item3','Noise+Speech_Set3_Item4','Noise+Speech_Set4_Item1','Noise+Speech_Set4_Item2','Noise+Speech_Set4_Item3','Noise+Speech_Set4_Item4','Noise+Speech_Set5_Item1','Noise+Speech_Set5_Item2','Noise+Speech_Set5_Item3','Noise+Speech_Set5_Item4','Noise+Speech_Set6_Item1','Noise+Speech_Set6_Item2','Noise+Speech_Set6_Item3','Noise+Speech_Set6_Item4','Noise+Speech_Set7_Item1','Noise+Speech_Set7_Item2','Noise+Speech_Set7_Item3','Noise+Speech_Set7_Item4','Noise+Speech_Set8_Item1','Noise+Speech_Set8_Item2','Noise+Speech_Set8_Item3','Noise+Speech_Set8_Item4'};
temp = load([GLMDir filesep 'SPM.mat']);
labelnames = {};
for i = 1:length(temp.SPM.Sess(1).U)
    if ~strncmp(temp.SPM.Sess(1).U(i).name,{'Match','Mismatch','Written'},5)
        continue
    else
        labelnames(end+1) = temp.SPM.Sess(1).U(i).name;
    end
end        
labels = 1:length(labelnames);

% set everything to calculate (dis)similarity estimates
cfg.decoding.software = 'distance';
cfg.decoding.method = 'classification';
cfg.decoding.train.classification.model_parameters = 'cveuclidean';

% This option averages across (dis)similarity matrices of each
% cross-validation iteration and across all cells of the lower diagonal
% (i.e. all distance comparisons). If you want the entire matrix, consider
% using 'other_average' which only averages across cross-validation
% iterations. Alternatively, you could use the output 'RSA_beta' which is
% more general purpose, but a little more complex.
cfg.results.output = 'other_average';

% These parameters carry out the multivariate noise normalization using the
% residuals
cfg.scale.method = 'cov'; % we scale by noise covariance
cfg.scale.estimation = 'separate'; % we scale all data for each run separately while iterating across searchlight spheres
cfg.scale.shrinkage = 'lw2'; % Ledoit-Wolf shrinkage retaining variances

% The crossnobis distance is identical to the cross-validated Euclidean
% distance after prewhitening (multivariate noise normalization). It has
% been shown that a good estimate for the multivariate noise is provided
% by the residuals of the first-level model, in addition with Ledoit-Wolf
% regularization. Here we calculate those residuals. If you have them
% available already, you can load them into misc.residuals using only the
% voxels from cfg.files.mask
[misc.residuals,cfg.files.residuals.chunk] = residuals_from_spm(fullfile(beta_loc,'SPM.mat'),cfg.files.mask); % this only needs to be run once and can be saved and loaded

% Set additional parameters manually if you want (see decoding.m or
% decoding_defaults.m). Below some example parameters that you might want
% to use:

cfg.searchlight.unit = 'mm';
cfg.searchlight.radius = 8; % this will yield a searchlight radius of 12mm.
cfg.searchlight.spherical = 1;
cfg.verbose = 2; % you want all information to be printed on screen
downsamp_ratio = 1;
cfg.searchlight.subset = 1:downsamp_ratio:1000000;

% Decide whether you want to see the searchlight/ROI/... during decoding
cfg.plot_selected_voxels = 100; % 0: no plotting, 1: every step, 2: every second step, 100: every hundredth step...

%% Nothing needs to be changed below for standard dissimilarity estimates using all data

% The following function extracts all beta names and corresponding run
% numbers from the SPM.mat
regressor_names = design_from_spm(beta_loc);

% Extract all information for the cfg.files structure (labels will be [1 -1] )
cfg = decoding_describe_data(cfg,labelnames,labels,regressor_names,beta_loc);

% This creates a design in which cross-validation is done between the distance estimates
cfg.design = make_design_similarity_cv(cfg);

% Run decoding
cfg.results.overwrite = 1;
try
    results = decoding(cfg,[],misc);
catch
    assert(~~exist([cfg.results.dir filesep 'res_other_average.mat'],'file'),'Something went wrong with the decoding - the results do not exist')
end

%% Make effect-maps (by correlating neural RDMs to model RDMs)

version = 'spearman'; % how to assess accuracy of model RDMs (pearson, spearman, weighted average)

outputDir = fullfile(GLMDir,'TDTcrossnobis',version);
if exist(outputDir,'dir'); rmdir(outputDir,'s'); mkdir(outputDir); else; mkdir(outputDir); end

clear models

basemodels.vowels = zeros(16,16);
basemodels.vowels(1:17:end) = 1;
basemodels.vowels(2:68:end) = 1/3;
basemodels.vowels(3:68:end) = 1/3;
basemodels.vowels(4:68:end) = 1/3;
basemodels.vowels(17:68:end) = 1/3;
basemodels.vowels(19:68:end) = 1/3;
basemodels.vowels(20:68:end) = 1/3;
basemodels.vowels(33:68:end) = 1/3;
basemodels.vowels(34:68:end) = 1/3;
basemodels.vowels(36:68:end) = 1/3;
basemodels.vowels(49:68:end) = 1/3;
basemodels.vowels(50:68:end) = 1/3;
basemodels.vowels(51:68:end) = 1/3;
basemodels.vowels = 1-basemodels.vowels;
basemodelNames = {'vowels'};

modeltemplate = NaN(size(results.other_average.output{1}));

labelnames_denumbered = {};
for i = 1:length(labelnames)
    labelnames_denumbered{i} = labelnames{i}(isletter(labelnames{i})|isspace(labelnames{i}));
end
modelNames = unique(labelnames_denumbered,'stable');

for i = 1:length(modelnames)
    models{i} = modeltemplate;
    models{i}(strcmp(modelnames{i},labelnames_denumbered),strcmp(modelnames{i},labelnames_denumbered))=basemodels.vowels
end

% % models = modelRDMs; close all
% %modelNames = fieldnames(models);
% 
% %modelNames = {'ProbM' 'ProbMM' 'EntropyM' 'EntropyMM'};

load(fullfile(GLMDir,'TDTcrossnobis','res_other_average.mat'));
data = results.other_average.output;

V = spm_vol(fullfile(GLMDir,'mask.nii'));
mask = spm_read_vols(V);
mask_index = results.mask_index;

clear results % to free memory

for m=1:length(modelNames)
    fprintf('\nComputing effect-map for model %s\n',modelNames{m});
    
    modelRDM = vectorizeRDMs(models{m})';
    effectMap = NaN(size(mask));
    for vx=1:numel(data)
        neuralRDM = vectorizeRDMs(data{vx})';
        if isempty(neuralRDM)
            effectMap(mask_index(vx)) = effectMap(mask_index(notempty)); %Duplicate values for downsampled searchlight.
            continue
        end
        notempty = vx;
        if ~isempty(strfind(version,'pearson'))
            effectMap(mask_index(vx)) = fisherTransform(corr(modelRDM,neuralRDM,'type','Pearson','Rows','pairwise'));
        elseif ~isempty(strfind(version,'spearman'))
            effectMap(mask_index(vx)) = fisherTransform(corr(modelRDM,neuralRDM,'type','Spearman','Rows','pairwise'));
        elseif ~isempty(strfind(version,'average'))
            %effectMap(mask_index(vx)) = mean(neuralRDM(find(~isnan(modelRDM)),:),1);
            effectMap(mask_index(vx)) = mean(neuralRDM(find(modelRDM==1),:),1);
        end
        if ~mod(vx,100)
            disp(['Processing voxel ' num2str(vx) ' of ' num2str(numel(data))])
        end
    end
    saveMRImage(effectMap,fullfile(outputDir,['effect-map_' modelNames{m} '.nii']),V.mat);
end

%% Normalise effect-maps to MNI template
% 
% clear version
% 
% cnt = 0;
% 
% cnt = cnt + 1;
% version{cnt} = 'spearman';
% 
% spm('FMRI');
% 
% for v=1:length(version)
%     
%     versionCurrent = version{v};
%     
%     clear matlabbatch
%     
%     load([TempPD '/Normalize.mat']);
%     
%     % Gather images for current subject
%     images = cellstr(spm_select('FPList', [GLMDir '/TDTcrossnobis/' versionCurrent '/'], '^effect-map_.*.nii'));
%     
%     % Write out masks
%     images_mask = {};
%     for i=1:length(images)
%         V = spm_vol(images{i});
%         Y = spm_read_vols(V);
%         Y(~isnan(Y)) = 1;
%         Y(isnan(Y)) = 0;
%         images_mask{i,1} = strrep(images{i},'effect-map','nativeSpaceMask');
%         saveMRImage(Y,images_mask{i,1},V.mat);
%     end
%     
%     % Normalize
%     DefFile = cellstr(spm_select('FPList', StrDir, '^y.*.nii'));
%     matlabbatch{1}.spm.spatial.normalise.write.subj.resample = [images; images_mask];
%     matlabbatch{1}.spm.spatial.normalise.write.subj.def = DefFile;
%     
%     save(fullfile(GLMDir,'TDTcrossnobis',versionCurrent,'NormalizeTDTcrossnobis.mat'), 'matlabbatch');
%     spm_jobman('initcfg')
%     spm_jobman('run', matlabbatch);
%     
%     clear images images_mask
%     
%     % Smooth, mask and write out normalised images
%     
%     % Gather images for current subject
%     images = cellstr(spm_select('FPList', [GLMDir '/TDTcrossnobis/' versionCurrent '/'], '^weffect-map_.*.nii'));
%     images_mask = cellstr(spm_select('FPList', [GLMDir '/TDTcrossnobis/' versionCurrent '/'], '^wnativeSpaceMask_.*.nii'));
%     
%     % Mask and smooth normalised data
%     mask_threshold = .05;
%     for i=1:length(images)
%         % Fix normalised mask
%         V = spm_vol(images_mask{i});
%         fname_mask = V.fname;
%         Y_mask = spm_read_vols(V);
%         Y_mask(Y_mask<mask_threshold) = 0;
%         Y_mask(isnan(Y_mask)) = 0;
%         saveMRImage(Y_mask,fname_mask,V.mat);
%         
%         % Fix normalised r-map
%         V = spm_vol(images{i});
%         fname_image = V.fname;
%         Y = spm_read_vols(V);
%         Y(Y_mask<mask_threshold) = 0;
%         Y_mask(isnan(Y_mask)) = 0;
%         saveMRImage(Y,fname_image,V.mat);
%         
%         % Smooth and (re)mask normalised r-map
%         fname_smoothed = strrep(images{i},'weffect-map_','sweffect-map_');
%         spm_smooth(images{i},fname_smoothed,[6 6 6]);
%         V = spm_vol(fname_smoothed);
%         Y = spm_read_vols(V);
%         Y(Y_mask<mask_threshold) = NaN;
%         saveMRImage(Y,fname_smoothed,V.mat);
%     end
%     
% end
%     
    
