clearvars

addpath(genpath('/group/language/data/ediz.sohoglu/matlab/rsatoolbox'));
addpath '/group/language/data/ediz.sohoglu/matlab/utilities';

%% Extract data

project_path = '/imaging/es03/fMRI_2017/';        
data_path = [project_path 'GLMAnalysisNative/'];
image_indices = [1:224];
%image_indices = [1:80 82:145];
for i=1:length(image_indices) % Don't include Hits/False Alarms and Misses
    image_filenames{i} = sprintf('spmT_%04d',image_indices(i));
    %image_filenames{i} = sprintf('con_%04d',image_indices(i));
end
roi_filenames_left = {'Left_HG_Oxford_30%' 'Left_PP_Oxford_30%' 'Left_PT_Oxford_30%' 'Left_STG_SylEvans' 'Left_STG_PosteriorPKBlank' 'Left_STG_AnteriorPKBlank'};
roi_filenames_right = {'Right_HG_Oxford_30%' 'Right_PP_Oxford_30%'  'Right_PT_Oxford_30%' 'Right_STG_SylEvans' 'Right_STG_PosteriorPKBlank' 'Right_STG_AnteriorPKBlank'};
roi_labels = {'HG' 'PP' 'PT' 'STG Central' 'STG Post' 'STG Ant'};

Subj = {'subj1' 'subj2' 'subj3' 'subj4' 'subj5' 'subj6' 'subj7' 'subj8' 'subj9' 'subj10' 'subj11' 'subj12' 'subj13' 'subj14' 'subj15' 'subj16' 'subj17' 'subj18' 'subj19' 'subj20' 'subj21' 'subj22' 'subj23' 'subj24' 'subj25' 'subj26' 'subj27'};
SubjToAnalyze = [3,4,5,6,8,10,11,12,13,15,16,17,18,19,20,21,22,23,24,26,27];

clear data nVoxels
for s=1:length(SubjToAnalyze)
    
    for i=1:length(image_filenames)
        
        for r=1:length(roi_filenames_left)
            
            Vi = spm_vol([data_path Subj{SubjToAnalyze(s)} '/' image_filenames{i} '.nii']);
            image = spm_read_vols(Vi);
            
            Vi = spm_vol([project_path 'PreprocessAnalysis/' Subj{SubjToAnalyze(s)} '/ROIs/w' roi_filenames_left{r} '.nii']); % Native
            roi = spm_read_vols(Vi);
            roi(isnan(image)) = 0;
            roi(roi<.01) = 0;
            nVoxels(s,r,1) = numel(find(roi));
            
            if sum(size(roi)==size(image))
                data{s}{r,1}(:,i) = image(find(roi));
            else
                error('Data and ROI image dimensions dont match');
            end
            
            Vi = spm_vol([project_path 'PreprocessAnalysis/' Subj{SubjToAnalyze(s)} '/ROIs/w' roi_filenames_right{r} '.nii']); % Native
            roi = spm_read_vols(Vi);
            roi(isnan(image)) = 0;
            roi(roi<.01) = 0;
            nVoxels(s,r,2) = numel(find(roi));
            
            if sum(size(roi)==size(image))
                data{s}{r,2}(:,i) = image(find(roi));
            else
                error('Data and ROI image dimensions dont match');
            end
            
        end
        
    end
    
end

%% Compute (dis)similarity matrix

dist_neural = {};
for s=1:length(SubjToAnalyze)
    
    for r=1:length(roi_filenames_left)
        
        dist_neural{s}{r,1} = 1-corr(data{s}{r,1},'type','Pearson','Rows','Pairwise');
        dist_neural{s}{r,2} = 1-corr(data{s}{r,2},'type','Pearson','Rows','Pairwise');
        
    end
    
end

%% RSA models

clear Models
tmp = modelRDMs_NaNforCrossSplice; close all;
Models(1).RDM = tmp.Syl2StrongM;
Models(2).RDM = tmp.Syl2WeakM;
Models(3).RDM = tmp.Syl2StrongMM;
Models(4).RDM = tmp.Syl2WeakMM;

% Convert to 'parametric' format?
for m=1:length(Models)
   Models(m).RDM(find(Models(m).RDM==1)) = 1/length(find(Models(m).RDM==1)); 
   Models(m).RDM(find(Models(m).RDM==0)) = -1/length(find(Models(m).RDM==0)); 
   Models(m).RDM(find(isnan(Models(m).RDM))) = 0; 
end

% clear Models
% tmp = modelRDMs_fromSimulations; close all;
% Models(1).RDM = tmp.PE;
% Models(2).RDM = tmp.PEabs;

%% Compute contrasts using RSA models

clear effect
for s=1:length(SubjToAnalyze)
    
    for r=1:length(roi_filenames_left)
        
        for m=1:length(Models)
            
            effect(s,r,1,m) = vectorizeRDMs(dist_neural{s}{r,1})*vectorizeRDMs(Models(m).RDM)';
            effect(s,r,2,m) = vectorizeRDMs(dist_neural{s}{r,2})*vectorizeRDMs(Models(m).RDM)';
            %effect(s,r,1,m) = fisherTransform(corr(vectorizeRDMs(dist_neural{s}{r,1})',vectorizeRDMs(Models(m).RDM)','type','spearman','rows','pairwise'));
            %effect(s,r,2,m) = fisherTransform(corr(vectorizeRDMs(dist_neural{s}{r,2})',vectorizeRDMs(Models(m).RDM)','type','spearman','rows','pairwise'));
            %effect(s,r,1,m) = fisherTransform(corr(vectorizeRDMs(dist_neural{s}{r,1}(1:128,1:128))',vectorizeRDMs(Models(m).RDM)','type','pearson','rows','pairwise'));
            %effect(s,r,2,m) = fisherTransform(corr(vectorizeRDMs(dist_neural{s}{r,2}(1:128,1:128))',vectorizeRDMs(Models(m).RDM)','type','pearson','rows','pairwise'));
                        
        end
        
    end
    
end

%% Plotting and statistics

figure;

means = squeeze(mean(effect,1));
sems = squeeze(std(effect,0,1)/sqrt(length(SubjToAnalyze)));

models = {'+M' '-M' '+MM' '-MM'};
%models = {'Syl2' 'Syl1'};
%models = {'Cat' 'Cat(N)' 'Graded' 'Graded(N)'};
%models = {'PE' 'PEabs'};
for hem=1:2
    subplot(1,2,hem);
    errorbar(squeeze(means(:,hem,:))',squeeze(sems(:,hem,:))','LineWidth',1.5); hold all
    %ylim([min(means(:)) max(means(:))]);
    set(gca,'xticklabel',models);
    set(gca,'xtick',[1:length(models)]);
    %title(roi_labels{r},'FontSize',15,'FontWeight','Bold');
    set(gca,'FontSize',15);
end
legend(roi_labels);

for r=1:length(roi_labels)
    fprintf('\nRegion %s\n',roi_labels{r});
    
    %[~,p_left] = ttest(squeeze(effect(:,r,1,5)),0,'tail','right');
    %[~,p_right] = ttest(squeeze(effect(:,r,2,5)),0,'tail','right');
    
    %fprintf('\nNoise+Speech: Left p = %0.3f, Right p = %0.3f\n',p_left,p_right);
    
    [~,p_left] = ttest(squeeze(mean(effect(:,r,1,:),4)),0,'tail','right');
    [~,p_right] = ttest(squeeze(mean(effect(:,r,2,:),4)),0,'tail','right');
    fprintf('\nAverage: Left p = %0.3f, Right p = %0.3f\n',p_left,p_right);
    
    repanova([squeeze(effect(:,r,1,1:4)) squeeze(effect(:,r,2,1:4))],[2 2 2],{'Hem' 'Cong' 'Strength'});
    %repanova(squeeze(effect(:,r,1,1:4)),[2 2],{'Cong' 'Strength'}); 
end