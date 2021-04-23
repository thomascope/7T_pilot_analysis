function module_make_effect_maps(GLMDir,downsamp_ratio)
%For taking already calculated crossnobis distances and doing RSA

if ~exist('downsamp_ratio','var')
    downsamp_ratio = 1;
end

addpath('/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/RSA_scripts/es_scripts_fMRI')
addpath('/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/RSA_scripts/decoding_toolbox_v3.999')
addpath(genpath('/group/language/data/ediz.sohoglu/matlab/rsatoolbox'));

%Define input data location
if downsamp_ratio == 1
    cfg.results.dir = fullfile(GLMDir,'TDTcrossnobis');
else
    cfg.results.dir = fullfile(GLMDir,['TDTcrossnobis_downsamp_' num2str(downsamp_ratio)]);
end

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

%% Make effect-maps (by correlating neural RDMs to model RDMs)

version = 'spearman'; % how to assess accuracy of model RDMs (pearson, spearman, weighted average)

if downsamp_ratio == 1
    outputDir = fullfile(GLMDir,'TDTcrossnobis',version);
else
    outputDir = fullfile(GLMDir,['TDTcrossnobis_downsamp_' num2str(downsamp_ratio)],version);
end
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

load(fullfile(cfg.results.dir,'res_other_average.mat'));
data = results.other_average.output;
notempty_data = find(~cellfun(@isempty,results.other_average.output));
modeltemplate = NaN(size(results.other_average.output{notempty_data(1)}));

labelnames_denumbered = {};
for i = 1:length(labelnames)
    labelnames_denumbered{i} = labelnames{i}(isletter(labelnames{i})|isspace(labelnames{i}));
end
modelNames = unique(labelnames_denumbered,'stable');

for i = 1:length(modelNames)
    models{i} = modeltemplate;
    models{i}(strcmp(modelNames{i},labelnames_denumbered),strcmp(modelNames{i},labelnames_denumbered))=basemodels.vowels;
end

% % models = modelRDMs; close all
% %modelNames = fieldnames(models);
% 
% %modelNames = {'ProbM' 'ProbMM' 'EntropyM' 'EntropyMM'};


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
    dims = size(effectMap);
    downsamped_effectMap = effectMap(1:downsamp_ratio:dims(1),1:downsamp_ratio:dims(2),1:downsamp_ratio:dims(3));
    downsamped_V.mat = V.mat;
    downsamped_V.mat(1:3,1:3)=downsamped_V.mat(1:3,1:3)*downsamp_ratio;
    
    saveMRImage(downsamped_effectMap,fullfile(outputDir,['effect-map_' modelNames{m} '.nii']),downsamped_V.mat);
end