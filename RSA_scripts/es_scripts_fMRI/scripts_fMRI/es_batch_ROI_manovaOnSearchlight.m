clearvars

addpath(genpath('/group/language/effect/ediz.sohoglu/matlab/rsatoolbox'));
addpath '/group/language/data/ediz.sohoglu/matlab/utilities';

%% Extract effect

project_path = '/imaging/es03/fMRI_2017/';        
data_path = [project_path 'GLMAnalysisNativeMasked/'];
% image_indices = [7:10];
% for i=1:length(image_indices)
%     image_filenames{i} = sprintf('swspmDs_C%04d_P0001',image_indices(i));
% end
images2analyse = {'Syl2StrongM' 'Syl2WeakM' 'Syl2StrongMM' 'Syl2WeakMM'};
for i=1:length(images2analyse)
    image_filenames{i} = sprintf('sweffect-map_%s',images2analyse{i});
end
%roi_filenames_left = {'Left_STG_SylEvans' 'Left_STG_PosteriorPKBlank' 'Left_STG_AnteriorPKBlank'};
%roi_filenames_right = {'Right_STG_SylEvans' 'Right_STG_PosteriorPKBlank' 'Right_STG_AnteriorPKBlank'};
%roi_labels = {'STG Central' 'STG Post' 'STG Ant'};
roi_filenames_left = {'Left_HG_ManovaSyl2_NativeMasked' 'Left_STG_ManovaSyl2_NativeMasked'};
roi_filenames_right = {'Left_HG_ManovaSyl2_NativeMasked' 'Left_STG_ManovaSyl2_NativeMasked'};
roi_labels = {'HG' 'STG'};

Subj = {'subj1' 'subj2' 'subj3' 'subj4' 'subj5' 'subj6' 'subj7' 'subj8' 'subj9' 'subj10' 'subj11' 'subj12' 'subj13' 'subj14' 'subj15' 'subj16' 'subj17' 'subj18' 'subj19' 'subj20' 'subj21' 'subj22' 'subj23' 'subj24' 'subj25' 'subj26' 'subj27'};
SubjToAnalyze = [3,4,5,6,8,10,11,12,13,15,16,17,18,19,20,21,22,23,24,26,27];

clear effect nVoxels
for s=1:length(SubjToAnalyze)
    
    for i=1:length(image_filenames)
        
        for r=1:length(roi_labels)
            
            Vi = spm_vol([data_path Subj{SubjToAnalyze(s)} '/' image_filenames{i} '.nii']);
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

models = {'StrongM' 'WeakM' 'StrongMM' 'WeakMM'};
for r=1:length(roi_labels)
    subplot(length(roi_labels),1,r);
    errorbar(squeeze(means(r,:,:))',squeeze(sems(r,:,:))','LineWidth',1.5); hold all
    %set(gca,'xticklabel',models);
    %set(gca,'xtick',[1:length(models)]);
    title(roi_labels{r});
    
    repanova(squeeze(effect(:,r,1,1:4)),[2 2]);
end
legend({'Left' 'Right'});