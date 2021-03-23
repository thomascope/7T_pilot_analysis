clearvars

addpath('/group/language/data/ediz.sohoglu/matlab/utilities');

%% Extract data

project_path = '/imaging/es03/fMRI_2017/';          
%data_path = [project_path 'GLMAnalysisMNISmooth6mm/'];
%data_path = [project_path 'GLMAnalysisNative/'];
data_path = [project_path 'GLMAnalysisMNISmooth6mmByItem/'];

%image_names = {'con_0001' 'con_0002' 'con_0003' 'con_0004'};
%roi_filenames_left = {'Left_STG_UnivariatePriorStrength'};
%roi_filenames_right = {'Right_STG_UnivariatePriorStrength'};
%roi_labels = {'STG'};

% image_names = {}; for i=1:128; image_names{i} = sprintf('con_%04d',i); end
% roi_filenames_left = {'Left_HG_Syl2Crossnobis' 'Left_STG_SylEvans' 'Left_STG_PosteriorPKBlank' 'Left_STG_AnteriorPKBlank'};
% roi_filenames_right = {'Right_HG_Syl2Crossnobis' 'Right_STG_SylEvans' 'Right_STG_PosteriorPKBlank' 'Right_STG_AnteriorPKBlank'};
% roi_labels = {'HG' 'STG Central' 'STG Post' 'STG Ant'};

image_names = {}; for i=1:192; image_names{i} = sprintf('con_%04d',i); end
roi_filenames_left = {'Left_HG_Syl2GradedNormCrossnobis'};
roi_filenames_right = {'Right_HG_Syl2GradedNormCrossnobis'};
roi_labels = {'HG'};

Subj = {'subj1' 'subj2' 'subj3' 'subj4' 'subj5' 'subj6' 'subj7' 'subj8' 'subj9' 'subj10' 'subj11' 'subj12' 'subj13' 'subj14' 'subj15' 'subj16' 'subj17' 'subj18' 'subj19' 'subj20' 'subj21' 'subj22' 'subj23' 'subj24' 'subj25' 'subj26' 'subj27'};
SubjToAnalyze = [3,4,5,6,8,10,11,12,13,15,16,17,18,19,20,21,22,23,24,26,27];

clear effect nVoxels
for s=1:length(SubjToAnalyze)
    
    for i=1:length(image_names)
        
        for r=1:length(roi_labels)
            
            Vi = spm_vol([data_path Subj{SubjToAnalyze(s)} '/' image_names{i} '.nii']);
            image = spm_read_vols(Vi);
            
            %Vi = spm_vol([project_path 'PreprocessAnalysis/' Subj{SubjToAnalyze(s)} '/ROIs/w' roi_filenames_left{r} '.nii']); % Native
            Vi = spm_vol([project_path 'ROIs/' roi_filenames_left{r} '.nii']); % MNI
            roi = spm_read_vols(Vi);
            roi(roi<.01) = 0;
            nVoxels(s,r,1) = numel(find(roi));
            
            effect(s,r,1,i) = nanmean(image(find(roi)));
            
            %Vi = spm_vol([project_path 'PreprocessAnalysis/' Subj{SubjToAnalyze(s)} '/ROIs/w' roi_filenames_right{r} '.nii']); % Native
            Vi = spm_vol([project_path 'ROIs/' roi_filenames_right{r} '.nii']); % MNI
            roi = spm_read_vols(Vi);
            roi(roi<.01) = 0;
            nVoxels(s,r,2) = numel(find(roi));
            
            effect(s,r,2,i) = nanmean(image(find(roi)));
            
        end
        
    end
    
end

%% Plotting and statistics

figure;

models = {'+M' '-M' '+MM' '-MM'};
for r=1:length(roi_labels)
    effectAdj = [];
    for hem=1:2
        effectAdj(:,r,hem,:) = es_removeBetween(squeeze(effect(:,r,hem,:)));
    end
    means = squeeze(mean(effectAdj(:,r,:,:),1));
    sems = squeeze(std(effectAdj(:,r,:,:),0,1)/sqrt(length(SubjToAnalyze)));

    subplot(1,length(roi_labels),r);
    %errorbar(means',sems','LineWidth',1.5); hold all
    barwitherr(sems',means','LineWidth',1.5); hold all
    set(gca,'xticklabel',models);
    set(gca,'xtick',[1:length(models)]);
    title(roi_labels{r},'FontSize',15,'FontWeight','Bold');
    set(gca,'FontSize',15);
end
legend({'Left' 'Right'});

for r=1:length(roi_labels)
    fprintf('\nRegion %s\n',roi_labels{r});
    
    % 2 X 2 X 2 anova (hemisphere X strength X congruency)
    repanova([squeeze(effect(:,r,1,1:4)) squeeze(effect(:,r,2,1:4))],[2 2 2],{'Hem' 'Cong' 'Strength'});
    % 2 X 2 anova (LH: strength X congruency)
    repanova(squeeze(effect(:,r,1,1:4)),[2 2],{'Cong' 'Strength'});
    % 2 X 2 anova (RH: strength X congruency)
    repanova(squeeze(effect(:,r,2,1:4)),[2 2],{'Cong' 'Strength'});

end


%% Multiple regression between stimulus predictor variables and fMRI

cnt = 0;
clear X
%cnt = cnt + 1;
%X = ones(128,cnt);
cnt = cnt + 1;
load(['./PredictorsFromSimSegVersion/predictor_sim_univariate_peAbs_syl1.mat'],'predictorForfMRI');
X(:,cnt) = predictorForfMRI;
%X(X(:,cnt)==0) = 1e-3;
%X(:,cnt) = -log(predictorForfMRI);
%X(:,cnt) = zscore(X(:,cnt));
cnt = cnt + 1;
load(['./PredictorsFromSimSegVersion/predictor_sim_univariate_peAbs_syl2seg2.mat'],'predictorForfMRI');
X(:,cnt) = predictorForfMRI;
%X(X(:,cnt)==0) = 1e-3;
%X(:,cnt) = -log(predictorForfMRI);
%X(:,cnt) = zscore(X(:,cnt));

% effectRegress = [];
% for s=1:length(SubjToAnalyze)
%     for r=1:length(roi_labels)
%         B = regress(squeeze(effect(s,r,1,:)),X);
%         effectRegress(s,r,1,:) = B(2:end);
%         B = regress(squeeze(effect(s,r,2,:)),X);
%         effectRegress(s,r,2,:) = B(2:end);
% 
%     end
% end

effectRegress = [];
for s=1:length(SubjToAnalyze)
    for r=1:length(roi_labels)
        for p=1:size(X,2)
            effectRegress(s,r,1,p) = corr(squeeze(effect(s,r,1,:)),X(:,p),'type','Spearman','Rows','pairwise');
            effectRegress(s,r,2,p) = corr(squeeze(effect(s,r,2,:)),X(:,p),'type','Spearman','Rows','pairwise');
        end
    end
end

figure;

predictorNames = {'peAbs'};
for r=1:length(roi_labels)
    means = squeeze(mean(effectRegress(:,r,:,:),1));
    sems = squeeze(std(effectRegress(:,r,:,:),0,1)/sqrt(length(SubjToAnalyze)));

    subplot(1,length(roi_labels),r);
    %errorbar(means',sems','LineWidth',1.5); hold all
    barwitherr(sems',means','LineWidth',1.5); hold all
    set(gca,'xticklabel',predictorNames );
    set(gca,'xtick',[1:length(predictorNames )]);
    xlim([0 length(predictorNames )+1]);
    title(roi_labels{r},'FontSize',15,'FontWeight','Bold');
    set(gca,'FontSize',15);
end
legend({'Left' 'Right'});

for r=1:length(roi_labels)
    fprintf('\nRegion %s\n',roi_labels{r});
    
    [~,p_left] = ttest(squeeze(effectRegress(:,r,1,:)),0,'tail','right')
    [~,p_right] = ttest(squeeze(effectRegress(:,r,2,:)),0,'tail','right')
    
    repanova([squeeze(effectRegress(:,r,1,:)) squeeze(effectRegress(:,r,2,:))],[2 length(predictorNames)],{'Hem' 'Predictor'});

end
