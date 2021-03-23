clearvars

addpath(genpath('/group/language/data/ediz.sohoglu/matlab/rsatoolbox'));
addpath '/group/language/data/ediz.sohoglu/matlab/utilities';

%% Extract data

project_path = '/imaging/es03/fMRI_2017/';        
data_path = [project_path 'GLMAnalysisNative/RSA_parametric_NaNforCrossSplice/Maps/'];
project_name = 'fMRI_2017';
image_names = {'Syl2StrongM' 'Syl2WeakM' 'Syl2StrongMM' 'Syl2WeakMM'};
%image_names = {'Syl2GradedStrongM' 'Syl2GradedWeakM' 'Syl2GradedStrongMM' 'Syl2GradedWeakMM' 'Syl2GradedNoiseSpeech'};
%image_names = {'Syl2' 'Syl2Graded'};
for i=1:length(image_names)
    image_filenames{i} = sprintf('sw%s_rMap_mask_%s',project_name,image_names{i});
end
roi_filenames_left = {'Left_HG_Oxford_30%' 'Left_PP_Oxford_30%' 'Left_PT_Oxford_30%' 'Left_STG_SylEvans' 'Left_STG_PosteriorPKBlank' 'Left_STG_AnteriorPKBlank'};
roi_filenames_right = {'Right_HG_Oxford_30%' 'Right_PP_Oxford_30%'  'Right_PT_Oxford_30%' 'Right_STG_SylEvans' 'Right_STG_PosteriorPKBlank' 'Right_STG_AnteriorPKBlank'};
roi_labels = {'HG' 'PP' 'PT' 'STG Central' 'STG Post' 'STG Ant'};

Subj = {'subj1' 'subj2' 'subj3' 'subj4' 'subj5' 'subj6' 'subj7' 'subj8' 'subj9' 'subj10' 'subj11' 'subj12' 'subj13' 'subj14' 'subj15' 'subj16' 'subj17' 'subj18' 'subj19' 'subj20' 'subj21' 'subj22' 'subj23' 'subj24' 'subj25' 'subj26' 'subj27'};
SubjToAnalyze = [3,4,5,6,8,10,11,12,13,15,16,17,18,19,20,21,22,23,24,26,27];

clear effect nVoxels
for s=1:length(SubjToAnalyze)
    
    for i=1:length(image_filenames)
        
        for r=1:length(roi_labels)
            
            Vi = spm_vol([data_path image_filenames{i} '_' Subj{SubjToAnalyze(s)} '.img']);
            image = spm_read_vols(Vi);
            
            %Vi = spm_vol([project_path 'PreprocessAnalysis/' Subj{SubjToAnalyze(s)} '/ROIs/w' roi_filenames_left{r} '.nii']); % Native
            Vi = spm_vol([project_path '/ROIs/' roi_filenames_left{r} '.nii']); % MNI
            roi = spm_read_vols(Vi);
            roi(roi<.01) = 0;
            nVoxels(s,r,1) = numel(find(roi));
            
            if sum(size(roi)==size(image))
                effect(s,r,1,i) = nanmean(nanmean(nanmean(image(find(roi)))));
            else
                error('Data and ROI image dimensions dont match');
            end
            
            %Vi = spm_vol([project_path 'PreprocessAnalysis/' Subj{SubjToAnalyze(s)} '/ROIs/w' roi_filenames_right{r} '.nii']); % Native
            Vi = spm_vol([project_path '/ROIs/' roi_filenames_right{r} '.nii']); % MNI
            roi = spm_read_vols(Vi);
            roi(roi<.01) = 0;
            nVoxels(s,r,2) = numel(find(roi));
            
            if sum(size(roi)==size(image))
                effect(s,r,2,i) = nanmean(nanmean(nanmean(image(find(roi)))));
            else
                error('Data and ROI image dimensions dont match');
            end            
        end
        
    end
        
end

%% Plotting and statistics

figure;

means = squeeze(mean(effect,1));
sems = squeeze(std(effect,0,1)/sqrt(length(SubjToAnalyze)));

models = {'+M' '-M' '+MM' '-MM'};
%models = {'Syl2' 'Syl2Graded'};
for hem=1:2
    subplot(1,2,hem);
    errorbar(squeeze(means(:,hem,:))',squeeze(sems(:,hem,:))','LineWidth',1.5); hold all
    ylim([min(means(:)) max(means(:))]);
    set(gca,'xticklabel',models);
    set(gca,'xtick',[1:length(models)]);
    %title(roi_labels{r},'FontSize',15,'FontWeight','Bold');
    set(gca,'FontSize',15);
end
legend(roi_labels);

for r=1:length(roi_labels)
    fprintf('\nRegion %s\n',roi_labels{r});
    
%     [~,p_left] = ttest(squeeze(effect(:,r,1,:)),0,'tail','right')
%     [~,p_right] = ttest(squeeze(effect(:,r,2,:)),0,'tail','right')
    
%     [~,p_left] = ttest(squeeze(effect(:,r,1,5)),0,'tail','right');
%     [~,p_right] = ttest(squeeze(effect(:,r,2,5)),0,'tail','right');
%     fprintf('\nNoise+Speech: Left p = %0.3f, Right p = %0.3f\n',p_left,p_right);
%     
    [~,p_left] = ttest(squeeze(mean(effect(:,r,1,:),4)),0,'tail','right');
    [~,p_right] = ttest(squeeze(mean(effect(:,r,2,:),4)),0,'tail','right');
    fprintf('\nAverage: Left p = %0.3f, Right p = %0.3f\n',p_left,p_right);
    
    repanova([squeeze(effect(:,r,1,1:4)) squeeze(effect(:,r,2,1:4))],[2 2 2],{'Hem' 'Cong' 'Strength'});
    %repanova(squeeze(effect(:,r,1,1:4)),[2 2],{'Cong' 'Strength'});
end
