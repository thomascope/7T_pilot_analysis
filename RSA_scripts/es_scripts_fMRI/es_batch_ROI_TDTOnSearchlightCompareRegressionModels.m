clearvars

addpath(genpath('/group/language/data/ediz.sohoglu/matlab/rsatoolbox'));
addpath '/group/language/data/ediz.sohoglu/matlab/utilities';

%%

% clear version
% 
% cnt = 0;
% 
% cnt = cnt + 1;
% version{cnt} = 'regressionSyl2GradedGatingWithin';
% cnt = cnt + 1;
% version{cnt} = 'regressionSyl2GradedGatingfMRIWithin';
% cnt = cnt + 1;
% version{cnt} = 'regressionSyl2GradedProbWithin';
% cnt = cnt + 1;
% version{cnt} = 'regressionSyl2GradedProb2Within';
% 
% cnt = cnt + 1;
% version{cnt} = 'regressionSyl2GatingWithin';
% cnt = cnt + 1;
% version{cnt} = 'regressionSyl2GatingfMRIWithin';
% cnt = cnt + 1;
% version{cnt} = 'regressionSyl2ProbWithin';
% cnt = cnt + 1;
% version{cnt} = 'regressionSyl2Prob2Within';

%%


clear version

cnt = 0;

cnt = cnt + 1;
version{cnt} = 'regressionFromSimBySegPEabs';

cnt = cnt + 1;
version{cnt} = 'regressionFromSimBySegSensory';

cnt = cnt + 1;
version{cnt} = 'regressionFromSimBySegPred';

cnt = cnt + 1;
version{cnt} = 'regressionFromSimByFeaPEabs';

cnt = cnt + 1;
version{cnt} = 'regressionFromSimByFeaSensory';

cnt = cnt + 1;
version{cnt} = 'regressionFromSimByFeaPred';

cnt = cnt + 1;
version{cnt} = 'regressionFromSimBySylPEabs';

cnt = cnt + 1;
version{cnt} = 'regressionFromSimBySylSensory';

cnt = cnt + 1;
version{cnt} = 'regressionFromSimBySylPred';

%% Extract data

project_path = '/imaging/es03/fMRI_2017/';        
data_path = [project_path 'GLMAnalysisNative/'];
project_name = 'fMRI_2017';
roi_filenames_left = {'Left_HG_Syl2GradedNormCrossnobis'};
roi_filenames_right = {'Right_HG_Syl2GradedNormCrossnobis'};
roi_labels = {'HG'};

Subj = {'subj1' 'subj2' 'subj3' 'subj4' 'subj5' 'subj6' 'subj7' 'subj8' 'subj9' 'subj10' 'subj11' 'subj12' 'subj13' 'subj14' 'subj15' 'subj16' 'subj17' 'subj18' 'subj19' 'subj20' 'subj21' 'subj22' 'subj23' 'subj24' 'subj25' 'subj26' 'subj27'};
SubjToAnalyze = [3,4,5,6,8,10,11,12,13,15,16,17,18,19,20,21,22,23,24,26,27];

clear effect nVoxels
for s=1:length(SubjToAnalyze)
    
    for v=1:length(version)
        
        for r=1:length(roi_labels)
            
            Vi = spm_vol([data_path '/' Subj{SubjToAnalyze(s)} '/TDTcrossnobis/' version{v} '/sweffect-map_R.nii']);
            image = spm_read_vols(Vi);
            
            Vi = spm_vol([project_path '/ROIs/' roi_filenames_left{r} '.nii']); % MNI
            roi = spm_read_vols(Vi);
            roi(roi<.01) = 0;
            nVoxels(s,r,1) = numel(find(roi));
            
            if sum(size(roi)==size(image))
                effect(s,r,1,v) = nanmean(nanmean(nanmean(image(find(roi)))));
            else
                error('Data and ROI image dimensions dont match');
            end
            
            Vi = spm_vol([project_path '/ROIs/' roi_filenames_right{r} '.nii']); % MNI
            roi = spm_read_vols(Vi);
            roi(roi<.01) = 0;
            nVoxels(s,r,2) = numel(find(roi));
            
            if sum(size(roi)==size(image))
                effect(s,r,2,v) = nanmean(nanmean(nanmean(image(find(roi)))));
            else
                error('Data and ROI image dimensions dont match');
            end
                       
        end
        
    end
        
end

%% Plotting and statistics


figure;

%models = {'G' 'Gf' 'C' 'S'};
models = {'SE' 'SS' 'SP' 'FE' 'FS' 'FP' 'SylE' 'SylS' 'SylP'};
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
    title(strrep(roi_labels{r},'_',' '),'FontSize',15,'FontWeight','Bold');
    set(gca,'FontSize',15);
    %ylim([0.015 .032]);
end
legend({'Left' 'Right'});
