function []=UnivariateCorrelations_1Subj(SubjID,GLMAnalPD,ProcDataDir,TempPD)

addpath /group/language/data/ediz.sohoglu/projects/fMRI_2017/scripts_fMRI
addpath /imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7219
addpath(genpath('/group/language/data/ediz.sohoglu/matlab/rsatoolbox'));

GLMDir = fullfile(GLMAnalPD,SubjID);

% %% Make effect-maps (correlations between stim and fMRI)
% 
% predictorDir = '/group/language/data/ediz.sohoglu/projects/fMRI_2017/scripts_fMRI';
% 
% outputDir = fullfile(GLMDir,'univariateCorrelations');
% if exist(outputDir,'dir'); rmdir(outputDir,'s'); mkdir(outputDir); else; mkdir(outputDir); end
% 
% load(fullfile(predictorDir,'predictor_prob_univariate.mat'),'regressor');
% reg.prob = regressor;
% load(fullfile(predictorDir,'predictor_entropy_univariate.mat'),'regressor');
% reg.entropy = regressor;
% 
% % get data
% image_names = {}; for i=1:192; image_names{i} = sprintf('con_%04d',i); end
% data = [];
% for i=1:length(image_names)
%     V = spm_vol(fullfile(GLMDir,[image_names{i} '.nii']));
%     data(i,:,:,:) = spm_read_vols(V);
% end
% 
% % prob
% 
% % correlations for Speech+M
% r = fisherTransform(corr(reg.prob(1:64),data(1:64,:),'type','Spearman','Rows','pairwise'));
% effectMap = reshape(r,V.dim);
% saveMRImage(effectMap,fullfile(outputDir,'effect-map_Prob_M.nii'),V.mat);
% 
% % correlations for Speech+MM
% r = fisherTransform(corr(reg.prob(65:128),data(65:128,:),'type','Spearman','Rows','pairwise'));
% effectMap = reshape(r,V.dim);
% saveMRImage(effectMap,fullfile(outputDir,'effect-map_Prob_MM.nii'),V.mat);
% 
% % correlations for Speech+Noise
% r = fisherTransform(corr(reg.prob(1:64),data(129:192,:),'type','Spearman','Rows','pairwise'));
% effectMap = reshape(r,V.dim);
% saveMRImage(effectMap,fullfile(outputDir,'effect-map_Prob_Noise.nii'),V.mat);
% 
% % entropy
% 
% % correlations for Speech+M
% r = fisherTransform(corr(reg.entropy(1:64),data(1:64,:),'type','Spearman','Rows','pairwise'));
% effectMap = reshape(r,V.dim);
% saveMRImage(effectMap,fullfile(outputDir,'effect-map_Entropy_M.nii'),V.mat);
% 
% % correlations for Speech+MM
% r = fisherTransform(corr(reg.entropy(65:128),data(65:128,:),'type','Spearman','Rows','pairwise'));
% effectMap = reshape(r,V.dim);
% saveMRImage(effectMap,fullfile(outputDir,'effect-map_Entropy_MM.nii'),V.mat);
% 
% % correlations for Speech+Noise
% r = fisherTransform(corr(reg.entropy(1:64),data(129:192,:),'type','Spearman','Rows','pairwise'));
% effectMap = reshape(r,V.dim);
% saveMRImage(effectMap,fullfile(outputDir,'effect-map_Entropy_Noise.nii'),V.mat);

% %% Make effect-maps (partial correlations between stim and fMRI)
% 
% predictorDir = '/group/language/data/ediz.sohoglu/projects/fMRI_2017/scripts_fMRI';
% 
% outputDir = fullfile(GLMDir,'univariatePartialCorrelations');
% if exist(outputDir,'dir'); rmdir(outputDir,'s'); mkdir(outputDir); else; mkdir(outputDir); end
% 
% clear reg
% load(fullfile(predictorDir,'predictor_prob_univariate.mat'),'regressor');
% reg.prob = regressor;
% load(fullfile(predictorDir,'predictor_entropy_univariate.mat'),'regressor');
% reg.entropy = regressor;
% 
% % get data
% image_names = {}; for i=1:192; image_names{i} = sprintf('con_%04d',i); end
% data = [];
% for i=1:length(image_names)
%     V = spm_vol(fullfile(GLMDir,[image_names{i} '.nii']));
%     data(i,:,:,:) = spm_read_vols(V);
% end
% 
% indGood = find(~isnan(squeeze(data(1,:,:,:))));
% 
% % prob
% 
% clear r
% r{1} = zeros(1,prod(V.dim));
% r{2} = zeros(1,prod(V.dim));
% r{3} = zeros(1,prod(V.dim));
% r{4} = zeros(1,prod(V.dim));
% r{5} = zeros(1,prod(V.dim));
% r{6} = zeros(1,prod(V.dim));
% for vx=1:1000:numel(indGood)
%     indCurrent = vx:vx+1000-1;
%     if indCurrent(end) > numel(indGood)
%         indCurrent = vx:numel(indGood);
%     end
%     
%     % prob
%     
%     % partial correlations for Speech+M
%     r{1}(indGood(indCurrent)) = fisherTransform(partialcorr(data(1:64,indGood(indCurrent)),reg.prob(1:64),reg.entropy(1:64),'type','Spearman','Rows','pairwise'));
%     % partial correlations for Speech+MM
%     r{2}(indGood(indCurrent)) = fisherTransform(partialcorr(data(65:128,indGood(indCurrent)),reg.prob(65:128),reg.entropy(65:128),'type','Spearman','Rows','pairwise'));  
%     % partial correlations for Speech+Noise
%     r{3}(indGood(indCurrent)) = fisherTransform(partialcorr(data(129:192,indGood(indCurrent)),reg.prob(1:64),reg.entropy(1:64),'type','Spearman','Rows','pairwise'));
%     
%     % entropy
%     
%     % partial correlations for Speech+M
%     r{4}(indGood(indCurrent)) = fisherTransform(partialcorr(data(1:64,indGood(indCurrent)),reg.entropy(1:64),reg.prob(1:64),'type','Spearman','Rows','pairwise'));
%     % partial correlations for Speech+MM
%     r{5}(indGood(indCurrent)) = fisherTransform(partialcorr(data(65:128,indGood(indCurrent)),reg.entropy(65:128),reg.prob(65:128),'type','Spearman','Rows','pairwise'));
%     % partial correlations for Speech+Noise
%     r{6}(indGood(indCurrent)) = fisherTransform(partialcorr(data(129:192,indGood(indCurrent)),reg.entropy(1:64),reg.prob(1:64),'type','Spearman','Rows','pairwise'));
% end
% 
% effectMap = reshape(r{1},V.dim);
% saveMRImage(effectMap,fullfile(outputDir,'effect-map_Prob_M.nii'),V.mat);
% effectMap = reshape(r{2},V.dim);
% saveMRImage(effectMap,fullfile(outputDir,'effect-map_Prob_MM.nii'),V.mat);
% effectMap = reshape(r{3},V.dim);
% saveMRImage(effectMap,fullfile(outputDir,'effect-map_Prob_Noise.nii'),V.mat);
% effectMap = reshape(r{4},V.dim);
% saveMRImage(effectMap,fullfile(outputDir,'effect-map_Entropy_M.nii'),V.mat);
% effectMap = reshape(r{5},V.dim);
% saveMRImage(effectMap,fullfile(outputDir,'effect-map_Entropy_MM.nii'),V.mat);
% effectMap = reshape(r{6},V.dim);
% saveMRImage(effectMap,fullfile(outputDir,'effect-map_Entropy_Noise.nii'),V.mat);

%% Make effect-maps (correlations between MEG and fMRI)

predictorDir = '/group/language/data/ediz.sohoglu/projects/fMRI_2017/scripts_fMRI';

outputDir = fullfile(GLMDir,'univariateCorrelationsMEGfMRI');
if exist(outputDir,'dir'); rmdir(outputDir,'s'); mkdir(outputDir); else; mkdir(outputDir); end

load(fullfile(predictorDir,'predictor_MEG_univariate.mat'),'regressor');
reg.MEG = regressor;

% get data
image_names = {}; for i=1:192; image_names{i} = sprintf('con_%04d',i); end
data = [];
for i=1:length(image_names)
    V = spm_vol(fullfile(GLMDir,[image_names{i} '.nii']));
    data(i,:,:,:) = spm_read_vols(V);
end

% correlations 
for hem = 1:2
    for t=1:size(reg.MEG,2)
        r = fisherTransform(corr(reg.MEG(1:64,t,hem),data(1:64,:),'type','Spearman','Rows','pairwise'));
        effectMap = reshape(r,V.dim);
        saveMRImage(effectMap,fullfile(outputDir,sprintf('effect-map_M_tWin%d_hem%d.nii',t,hem)),V.mat);
        r = fisherTransform(corr(reg.MEG(65:128,t,hem),data(65:128,:),'type','Spearman','Rows','pairwise'));
        effectMap = reshape(r,V.dim);
        saveMRImage(effectMap,fullfile(outputDir,sprintf('effect-map_MM_tWin%d_hem%d.nii',t,hem)),V.mat);
        r = fisherTransform(corr(reg.MEG(1:128,t,hem),data(1:128,:),'type','Spearman','Rows','pairwise'));
        effectMap = reshape(r,V.dim);
        saveMRImage(effectMap,fullfile(outputDir,sprintf('effect-map_All_tWin%d_hem%d.nii',t,hem)),V.mat);
    end
end

%% Make effect-maps (commonality between stim, MEG and fMRI)
% 
% predictorDir = '/group/language/data/ediz.sohoglu/projects/fMRI_2017/scripts_fMRI';
% 
% outputDir = fullfile(GLMDir,'univariateCommonality');
% if exist(outputDir,'dir'); rmdir(outputDir,'s'); mkdir(outputDir); else; mkdir(outputDir); end
% 
% clear reg
% load(fullfile(predictorDir,'predictor_prob_univariate.mat'),'regressor');
% reg.stim = regressor;
% load(fullfile(predictorDir,'predictor_MEG_univariate.mat'),'regressor');
% reg.meg = regressor;
% 
% % get data
% image_names = {}; for i=1:192; image_names{i} = sprintf('con_%04d',i); end
% data = [];
% for i=1:length(image_names)
%     V = spm_vol(fullfile(GLMDir,[image_names{i} '.nii']));
%     data(i,:,:,:) = spm_read_vols(V);
% end
% 
% % commonality analysis for Speech+M
% R2c = es_varDecomp(data(1:64,:),[reg.stim(1:64) reg.meg(1:64)]); 
% effectMap = reshape(R2c,V.dim);
% saveMRImage(effectMap,fullfile(outputDir,'effect-map_M.nii'),V.mat);
% 
% % commonality analysis for Speech+MM
% R2c = es_varDecomp(data(1:64,:),[reg.stim(65:128) reg.meg(65:128)]); 
% effectMap = reshape(R2c,V.dim);
% saveMRImage(effectMap,fullfile(outputDir,'effect-map_MM.nii'),V.mat);
% 
% % commonality analysis for Speech+Noise
% R2c = es_varDecomp(data(1:64,:),[reg.stim(1:64) reg.meg(129:192)]);
% effectMap = reshape(R2c,V.dim);
% saveMRImage(effectMap,fullfile(outputDir,'effect-map_Noise.nii'),V.mat);
