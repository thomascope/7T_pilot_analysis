clearvars

addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7219');
addpath(genpath('/group/language/data/ediz.sohoglu/matlab/rsatoolbox'));
addpath '/group/language/data/ediz.sohoglu/matlab/utilities';

%% Extract data

project_path = '/imaging/es03/fMRI_2017/';        
data_path = [project_path 'GLMAnalysisNative/'];
project_name = 'fMRI_2017';
version = {'average'};
%image_names = {'Syl2GradedNormStrongM' 'Syl2GradedNormWeakM' 'Syl2GradedNormStrongMM' 'Syl2GradedNormWeakMM'};
%image_names = {'Syl2GradedNormStrongM' 'Syl2GradedNormWeakM' 'Syl2GradedNormStrongMM' 'Syl2GradedNormWeakMM' 'Syl2GradedNormStrongNoise' 'Syl2GradedNormWeakNoise' 'Syl2GradedNormNoiseSpeech'};
image_names = {'Syl2StrongM' 'Syl2WeakM' 'Syl2StrongMM' 'Syl2WeakMM'};
%image_names = {'Syl1GradedNormStrong' 'Syl1GradedNormWeak'};
%image_names = {'ProbM' 'ProbMM'};
% roi_filenames_left = {'Left_HG_Syl2GradedNormCrossnobis' 'Left_STG_SylEvans' 'Left_STG_PosteriorPKBlank' 'Left_STG_AnteriorPKBlank'};
% roi_filenames_right = {'Right_HG_Syl2GradedNormCrossnobis' 'Right_STG_SylEvans' 'Right_STG_PosteriorPKBlank' 'Right_STG_AnteriorPKBlank'};
% roi_labels = {'HG' 'STG Central' 'STG Post' 'STG Ant'};
roi_filenames_left = {'Left_HG_Syl2GradedNormCrossnobis'};
roi_filenames_right = {'Right_HG_Syl2GradedNormCrossnobis'};
roi_labels = {'HG'};
% roi_filenames_left = {'Left_HG_Syl2Crossnobis' 'Left_PP_Syl2Crossnobis'};
% roi_filenames_right = {'Right_HG_Syl2Crossnobis' 'Left_PP_Syl2Crossnobis'};
% roi_labels = {'HG' 'PP'};
% roi_filenames_left = {'Right_STG_M-MM_univariatePartialCorrelation'};
% roi_filenames_right = {'Right_STG_M-MM_univariatePartialCorrelation'};
% roi_labels = {'STG'};

Subj = {'subj1' 'subj2' 'subj3' 'subj4' 'subj5' 'subj6' 'subj7' 'subj8' 'subj9' 'subj10' 'subj11' 'subj12' 'subj13' 'subj14' 'subj15' 'subj16' 'subj17' 'subj18' 'subj19' 'subj20' 'subj21' 'subj22' 'subj23' 'subj24' 'subj25' 'subj26' 'subj27'};
SubjToAnalyze = [3,4,5,6,8,10,11,12,13,15,16,17,18,19,20,21,22,23,24,26,27];

clear effect nVoxels
for s=1:length(SubjToAnalyze)
    
    for r=1:length(roi_labels)
               
        counter = 0;
        for v=1:length(version)
            
            for i=1:length(image_names)
                
                counter = counter + 1;
                
                Vi = spm_vol([data_path '/' Subj{SubjToAnalyze(s)} '/TDTcrossnobis/' version{v} '/sweffect-map_' image_names{i} '.nii']);
                image = spm_read_vols(Vi);
                
                %Vi = spm_vol([project_path 'PreprocessAnalysis/' Subj{SubjToAnalyze(s)} '/ROIs/w' roi_filenames_left{r} '.nii']); % Native
                Vi = spm_vol([project_path '/ROIs/' roi_filenames_left{r} '.nii']); % MNI
                roi = spm_read_vols(Vi);
                roi(roi<.01) = 0;
                nVoxels(s,r,1) = numel(find(roi));
                
                if sum(size(roi)==size(image))
                    effect(s,r,1,counter) = nanmean(nanmean(nanmean(image(find(roi)))));
                else
                    error('Data and ROI image dimensions dont match');
                end
                
                %Vi = spm_vol([project_path 'PreprocessAnalysis/' Subj{SubjToAnalyze(s)} '/ROIs/w' roi_filenames_right{r} '.nii']); % Native
                Vi = spm_vol([project_path '/ROIs/' roi_filenames_right{r} '.nii']); % MNI
                roi = spm_read_vols(Vi);
                roi(roi<.01) = 0;
                nVoxels(s,r,2) = numel(find(roi));
                
                if sum(size(roi)==size(image))
                    effect(s,r,2,counter) = nanmean(nanmean(nanmean(image(find(roi)))));
                else
                    error('Data and ROI image dimensions dont match');
                end
            end
            
        end
        
    end
    
end

%% Plotting and statistics

figure;

models = {'+M' '-M' '+MM' '-MM'};
%models = {'+M' '-M' '+MM' '-MM' '+N' '-N' 'NSp'};
%models = {'+N' '-N'};
%models = {'Sensory' 'Pred' 'PEabs' 'Ps' 'Sensory' 'Pred' 'PEabs' 'Ps' 'Sensory' 'Pred' 'PEabs' 'Ps'};
%models = {'SensorySyl1' 'PredSyl1' 'PEabsSyl1' 'PsSyl1' 'SensorySyl2' 'PredSyl2' 'PEabsSyl2' 'PsSyl2'};
for r=1:length(roi_labels)
    effectAdj = [];
    for hem=1:2
        effectAdj(:,r,hem,:) = es_removeBetween(squeeze(effect(:,r,hem,:)));
    end
    means = squeeze(mean(effectAdj(:,r,:,:),1));
    sems = squeeze(std(effectAdj(:,r,:,:),0,1)/sqrt(length(SubjToAnalyze)));

    subplot(1,length(roi_labels),r);
    errorbar(means',sems','LineWidth',1.5); hold all
    set(gca,'xticklabel',models);
    set(gca,'xtick',[1:length(models)]);
    title(roi_labels{r},'FontSize',15,'FontWeight','Bold');
    set(gca,'FontSize',15);
end
legend({'Left' 'Right'});

for r=1:length(roi_labels)
    fprintf('\nRegion %s\n',roi_labels{r});
    
    % t-test against zero for each condition
    [~,p_left] = ttest(squeeze(effect(:,r,1,:)),0,'tail','right')
    [~,p_right] = ttest(squeeze(effect(:,r,2,:)),0,'tail','right')
    
    % 2 X 2 X 2 anova (hemisphere X strength X congruency)
    repanova([squeeze(effect(:,r,1,1:4)) squeeze(effect(:,r,2,1:4))],[2 2 2],{'Hem' 'Cong' 'Strength'});
    % 2 X 2 anova (LH: strength X congruency)
    repanova(squeeze(effect(:,r,1,1:4)),[2 2],{'Cong' 'Strength'});
    % 2 X 2 anova (RH: strength X congruency)
    repanova(squeeze(effect(:,r,2,1:4)),[2 2],{'Cong' 'Strength'});  
%     % Strong-Weak (M)
%     repanova([squeeze(effect(:,r,1,1:2)) squeeze(effect(:,r,2,1:2))],[2 2],{'Hem' 'Strength'});
%     % Strong-Weak (MM)
%     repanova([squeeze(effect(:,r,1,3:4)) squeeze(effect(:,r,2,3:4))],[2 2],{'Hem' 'Strength'});
    % 2 X 2 anova (hemisphere X Strong+Noise/Weak+noise)
%     repanova([squeeze(effect(:,r,1,5:6)) squeeze(effect(:,r,2,5:6))],[2 2],{'Hem' 'Syl1type+Noise'});
%     % 2 X 3 anova (hemisphere X Strong+Speech/Weak+Speech/Noise+Speech)
%     clear contrast
%     contrast{1} = [.5 0 .5 0 zeros(1,3)];
%     contrast{2} = [0 .5 0 .5 zeros(1,3)];
%     contrast{3} = [zeros(1,6) 1];
%     effectTmp = [];
%     for cont=1:length(contrast)
%         effectTmp(:,cont) = squeeze(effect(:,r,1,:))*contrast{cont}';
%     end
%     for cont=1:length(contrast)
%         effectTmp(:,cont+length(contrast)) = squeeze(effect(:,r,2,:))*contrast{cont}';
%     end
%     repanova(effectTmp,[2 3],{'Hem' 'Syl1type+Clear'});
    
%     % Model comparison of segments, features and syllables
%     repanova([squeeze(effect(:,r,1,:)) squeeze(effect(:,r,2,:))],[2 3 4]);

end
