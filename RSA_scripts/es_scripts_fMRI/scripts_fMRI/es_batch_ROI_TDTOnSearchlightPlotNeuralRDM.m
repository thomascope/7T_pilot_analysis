clearvars

PreProcPD = '/imaging/es03/fMRI_2017/PreprocessAnalysis';
GLMAnalPD = '/imaging/es03/fMRI_2017/GLMAnalysisNative';

Subj = {'subj1' 'subj2' 'subj3' 'subj4' 'subj5' 'subj6' 'subj7' 'subj8' 'subj9' 'subj10' 'subj11' 'subj12' 'subj13' 'subj14' 'subj15' 'subj16' 'subj17' 'subj18' 'subj19' 'subj20' 'subj21' 'subj22' 'subj23' 'subj24' 'subj25' 'subj26' 'subj27'};
SubjToAnalyze = [3,4,5,6,8,10,11,12,13,15,16,17,18,19,20,21,22,23,24,26,27];

roi_filenames_left = {'Left_HG_Syl2Crossnobis' 'Left_STG_SylEvans' 'Left_STG_PosteriorPKBlank' 'Left_STG_AnteriorPKBlank'};
roi_filenames_right = {'Right_HG_Syl2Crossnobis' 'Right_STG_SylEvans' 'Right_STG_PosteriorPKBlank' 'Right_STG_AnteriorPKBlank'};
roi_labels = {'HG' 'STG Central' 'STG Post' 'STG Ant'};

clear neuralRDM
for s=1:length(SubjToAnalyze)
    
    for r=1:length(roi_labels)
        
        SubjCurrent = Subj{SubjToAnalyze(s)};
        SubjDir = fullfile(PreProcPD,SubjCurrent);
        ROIDir = fullfile(SubjDir, 'ROIs');
        
        GLMDir = fullfile(GLMAnalPD,SubjCurrent);
        load(fullfile(GLMDir,'TDTcrossnobis','res_other_average.mat'));
        data = results.other_average.output;
        mask_index = results.mask_index;
        
        V = spm_vol(fullfile(ROIDir,['w' roi_filenames_left{r} '.nii']));
        roi = spm_read_vols(V);
        [~,ind] = intersect(mask_index,find(roi));
        
        neuralRDM(:,:,r,1,s) = nanmean(cat(3,data{ind}),3);
        
        V = spm_vol(fullfile(ROIDir,['w' roi_filenames_right{r} '.nii']));
        roi = spm_read_vols(V);
        [~,ind] = intersect(mask_index,find(roi));
        
        neuralRDM(:,:,r,2,s) = nanmean(cat(3,data{ind}),3);
        
    end
    
end

%%

% neuralRDM_old = neuralRDM;
% neuralRDM(:,:,:,:,2) = [];
% SubjToAnalyze(2) = [];

%%

roi2plot = 1;
hem2plot = 1;

figure;
for m=1:length(modelNames)
    [~,~,~,stats] = ttest(permute(squeeze(neuralRDM([1:32]+(32*(m-1)),[1:32]+(32*(m-1)),roi2plot,hem2plot,:)),[3 1 2]));
    subplot(1,length(modelNames),m); imagesc(squeeze(stats.tstat),[-5 5]);
    title(modelNames{m});
end

%%

% modelNames = {'Syl2StrongM' 'Syl2WeakM' 'Syl2StrongMM' 'Syl2WeakMM' 'Syl2StrongNoise' 'Syl2WeakNoise' 'Syl2NoiseSpeech'};
% modelNamesForGraph = {'+M' '-M' '+MM' '-MM' '+N' '-N' 'NSp'};
modelNames = {'Syl1Strong' 'Syl1Weak'};
modelNamesForGraph = {'+' '-'};

clear models
models = modelRDMs; close all

effect = [];
for s=1:length(SubjToAnalyze)
    
    for r=1:length(roi_labels)
        
        for m=1:length(modelNames)
            
            subj2keep = setdiff(1:length(SubjToAnalyze),s);
            %subj2keep = 1:length(SubjToAnalyze);
            
            effect(s,r,1,m) = corr(vectorizeRDMs(neuralRDM(:,:,r,1,s))',vectorizeRDMs(mean(neuralRDM(:,:,r,1,subj2keep),5).*~isnan(models.(modelNames{m})))','type','spearman','rows','pairwise');
            effect(s,r,2,m) = corr(vectorizeRDMs(neuralRDM(:,:,r,2,s))',vectorizeRDMs(mean(neuralRDM(:,:,r,2,subj2keep),5).*~isnan(models.(modelNames{m})))','type','spearman','rows','pairwise');
            %effect(s,r,1,m) = corr(vectorizeRDMs(neuralRDM(:,:,r,1,s))',vectorizeRDMs(mean(neuralRDM(:,:,r,1,subj2keep),5))','type','spearman','rows','pairwise');
            %effect(s,r,2,m) = corr(vectorizeRDMs(neuralRDM(:,:,r,2,s))',vectorizeRDMs(mean(neuralRDM(:,:,r,2,subj2keep),5))','type','spearman','rows','pairwise');
            
        end
        
    end
    
end

effectAdj = [];
for r=1:length(roi_labels)
    effectAdj(:,r,1,:) = es_removeBetween(squeeze(effect(:,r,1,:)));
    effectAdj(:,r,2,:) = es_removeBetween(squeeze(effect(:,r,2,:)));
end

means = squeeze(mean(effectAdj,1));
sems = squeeze(std(effectAdj,0,1) ./ sqrt(length(SubjToAnalyze)));

figure; errorbar(means',sems','LineWidth',1.5);
set(gca,'xtick',1:length(modelNamesForGraph));
set(gca,'xticklabel',modelNamesForGraph);
set(gca,'FontSize',15);

for r=1:length(roi_labels)
    fprintf('\n\nROI %d\n\n',r);
    repanova([squeeze(effect(:,r,1,1:4)) squeeze(effect(:,r,2,1:4))],[2 2 2],{'Hem' 'Cong' 'Strength'});
    repanova(squeeze(effect(:,r,1,1:4)),[2 2],{'Cong' 'Strength'});
    repanova(squeeze(effect(:,r,2,1:4)),[2 2],{'Cong' 'Strength'});
end

