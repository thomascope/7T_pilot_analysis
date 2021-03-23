clearvars

addpath /group/language/data/ediz.sohoglu/matlab/cvMANOVA_v3
addpath(genpath('/group/language/data/ediz.sohoglu/matlab/rsatoolbox'));
addpath '/group/language/data/ediz.sohoglu/matlab/utilities';

%% Extract data

project_path = '/imaging/es03/fMRI_2017/';        
data_path = [project_path 'GLMAnalysisNative/'];
% roi_filenames_left = {'Left_HG_Oxford_30%' 'Left_PP_Oxford_30%' 'Left_PT_Oxford_30%' 'Left_STG_SylEvans' 'Left_STG_PosteriorPKBlank' 'Left_STG_AnteriorPKBlank'};
% roi_filenames_right = {'Right_HG_Oxford_30%' 'Right_PP_Oxford_30%'  'Right_PT_Oxford_30%' 'Right_STG_SylEvans' 'Right_STG_PosteriorPKBlank' 'Right_STG_AnteriorPKBlank'};
% roi_labels = {'HG' 'PP' 'PT' 'STG Central' 'STG Post' 'STG Ant'};
roi_filenames_left = {'Left_HG_Syl2GradedCrossnobis'};
roi_filenames_right = {'Right_HG_Syl2GradedCrossnobis'};
roi_labels = {'HG'};

Subj = {'subj1' 'subj2' 'subj3' 'subj4' 'subj5' 'subj6' 'subj7' 'subj8' 'subj9' 'subj10' 'subj11' 'subj12' 'subj13' 'subj14' 'subj15' 'subj16' 'subj17' 'subj18' 'subj19' 'subj20' 'subj21' 'subj22' 'subj23' 'subj24' 'subj25' 'subj26' 'subj27'};
SubjToAnalyze = [3,4,5,6,8,10,11,12,13,15,16,17,18,19,20,21,22,23,24,26,27];

[X,a] = indicatorMatrix('allpairs',1:128);
for c=1:size(X,1); Cs{c} = X(c,:)'; end

clear dist_neural nVoxels
for s=1:length(SubjToAnalyze)
    
    clear RoiDir
    for r=1:length(roi_filenames_left)
        RoiDir{r} = [project_path 'PreprocessAnalysis/' Subj{SubjToAnalyze(s)} '/ROIs/w' roi_filenames_left{r} '.nii'];
    end
    
    [data, nVox] = cvManovaRegion([data_path Subj{SubjToAnalyze(s)}], RoiDir, Cs, 0);
    for r=1:size(data,3)
        dist_neural{s}{r,1} = squeeze(data(:,:,r))';
        nVoxels(s,r,1) = nVox(r);
        dist_neural{s}{r,1} = dist_neural{s}{r,1} ./ sqrt(nVoxels(s,r,1));
    end
    
    clear RoiDir
    for r=1:length(roi_filenames_left)
        RoiDir{r} = [project_path 'PreprocessAnalysis/' Subj{SubjToAnalyze(s)} '/ROIs/w' roi_filenames_right{r} '.nii'];
    end
    
    [data, nVox] = cvManovaRegion([data_path Subj{SubjToAnalyze(s)}], RoiDir, Cs, 0);
    for r=1:size(data,3)
        dist_neural{s}{r,2} = squeeze(data(:,:,r))';
        nVoxels(s,r,2) = nVox(r);
        dist_neural{s}{r,2} = dist_neural{s}{r,2} ./ sqrt(nVoxels(s,r,2));
    end
    
end

%% Plot neural RDM

for s=1:length(SubjToAnalyze)
    dist2plot(:,:,s) = squareform(dist_neural{s}{1,1});
    figure(100); subplot(4,6,s); imagesc(dist2plot(:,:,s)); set(gca,'clim',[-1e-03 1e-03]);
end

figure(101); imagesc(mean(dist2plot,3)); %set(gca,'clim',[-1e-03 1e-03]);

%% RSA models

clear Models
tmp = modelRDMs; close all
%tmp = modelRDMsReducedWithSyl1Graded; close all
Models(1).RDM = tmp.Syl2GradedStrongM;
Models(2).RDM = tmp.Syl2GradedWeakM;
Models(3).RDM = tmp.Syl2GradedStrongMM;
Models(4).RDM = tmp.Syl2GradedWeakMM;

% Convert to 'parametric' format?
for m=1:length(Models)
   Models(m).RDM(find(Models(m).RDM==1)) = 1/length(find(Models(m).RDM==1)); 
   Models(m).RDM(find(Models(m).RDM==0)) = -1/length(find(Models(m).RDM==0)); 
   Models(m).RDM(find(isnan(Models(m).RDM))) = 0; 
end

%% Compute contrasts using RSA models

clear effect
for s=1:length(SubjToAnalyze)
    
    for r=1:length(roi_filenames_left)
        
        for m=1:length(Models)
            
            effect(s,r,1,m) = vectorizeRDMs(dist_neural{s}{r,1})*vectorizeRDMs(Models(m).RDM(1:128,1:128))';
            effect(s,r,2,m) = vectorizeRDMs(dist_neural{s}{r,2})*vectorizeRDMs(Models(m).RDM(1:128,1:128))';
            %effect(s,r,1,m) = fisherTransform(corr(vectorizeRDMs(dist_neural{s}{r,1})',vectorizeRDMs(Models(m).RDM)','type','pearson','rows','pairwise'));
            %effect(s,r,2,m) = fisherTransform(corr(vectorizeRDMs(dist_neural{s}{r,2})',vectorizeRDMs(Models(m).RDM)','type','pearson','rows','pairwise'));
                        
        end
        
    end
    
end

%% Exlude subjects?

ind2exclude = [];

effect(ind2exclude,:,:,:) = [];

%% Plotting and statistics

figure;

models = {'+M' '-M' '+MM' '-MM'};
for r=1:length(roi_labels)
    means = squeeze(mean(effect(:,r,:,:),1));
    sems = squeeze(std(effect(:,r,:,:),0,1)/sqrt(length(SubjToAnalyze)));

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
    
    % LH: Strong<Weak (M)
    [~,p_left] = ttest(squeeze(effect(:,r,1,2)-effect(:,r,1,1)),0,'tail','right')
    % RH: Strong<Weak (M)
    [~,p_right] = ttest(squeeze(effect(:,r,2,2)-effect(:,r,2,1)),0,'tail','right')
    % LH: Strong>Weak (MM)
    [~,p_left] = ttest(squeeze(effect(:,r,1,3)-effect(:,r,1,4)),0,'tail','right')
    % RH: Strong>Weak (MM)
    [~,p_right] = ttest(squeeze(effect(:,r,2,3)-effect(:,r,2,4)),0,'tail','right')
    
    % 2 X 2 anova (hemisphere X strength)
    repanova([squeeze(effect(:,r,1,5:6)) squeeze(effect(:,r,2,5:6))],[2 2],{'Hem' 'Strength'});
    % 2 anova (LH: strength)
    repanova(squeeze(effect(:,r,1,5:6)),[2],{'Strength'});
    % 2 anova (RH: strength)
    repanova(squeeze(effect(:,r,2,5:6)),[2],{'Strength'});

end

