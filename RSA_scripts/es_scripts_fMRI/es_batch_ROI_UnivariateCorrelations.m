clearvars

addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7219');
addpath(genpath('/group/language/data/ediz.sohoglu/matlab/rsatoolbox'));
addpath '/group/language/data/ediz.sohoglu/matlab/utilities';

%% Extract data

project_path = '/imaging/es03/fMRI_2017/';        
data_path = [project_path 'GLMAnalysisMNISmooth6mmByItem/'];
project_name = 'fMRI_2017';
image_names = {'Prob_M' 'Prob_MM'};
roi_filenames_left = {'Right_STG_M-MM_univariatePartialCorrelation'};
roi_filenames_right = {'Right_STG_M-MM_univariatePartialCorrelation'};
roi_labels = {'STG'};

Subj = {'subj1' 'subj2' 'subj3' 'subj4' 'subj5' 'subj6' 'subj7' 'subj8' 'subj9' 'subj10' 'subj11' 'subj12' 'subj13' 'subj14' 'subj15' 'subj16' 'subj17' 'subj18' 'subj19' 'subj20' 'subj21' 'subj22' 'subj23' 'subj24' 'subj25' 'subj26' 'subj27'};
SubjToAnalyze = [3,4,5,6,8,10,11,12,13,15,16,17,18,19,20,21,22,23,24,26,27];

clear effect nVoxels
for s=1:length(SubjToAnalyze)
    
    for r=1:length(roi_labels)
            
        for i=1:length(image_names)
                        
            Vi = spm_vol([data_path '/' Subj{SubjToAnalyze(s)} '/univariatePartialCorrelations/effect-map_' image_names{i} '.nii']);
            image = spm_read_vols(Vi);
            
            Vi = spm_vol([project_path '/ROIs/' roi_filenames_left{r} '.nii']); % MNI
            roi = spm_read_vols(Vi);
            roi(roi<.01) = 0;
            nVoxels(s,r,1) = numel(find(roi));
            
            if sum(size(roi)==size(image))
                effect(s,r,1,i) = nanmean(nanmean(nanmean(image(find(roi)))));
            else
                error('Data and ROI image dimensions dont match');
            end
            
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

models = {'M' 'MM'};
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
    [~,p_left] = ttest(squeeze(effect(:,r,1,:)))
    [~,p_right] = ttest(squeeze(effect(:,r,2,:)))

end
