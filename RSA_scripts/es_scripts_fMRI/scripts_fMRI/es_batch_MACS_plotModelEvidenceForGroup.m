clearvars

Subj = {'subj1' 'subj2' 'subj3' 'subj4' 'subj5' 'subj6' 'subj7' 'subj8' 'subj9' 'subj10' 'subj11' 'subj12' 'subj13' 'subj14' 'subj15' 'subj16' 'subj17' 'subj18' 'subj19' 'subj20' 'subj21' 'subj22' 'subj23' 'subj24' 'subj25' 'subj26' 'subj27'};
SubjToAnalyze = [3,4,5,6,8,10,11,12,13,15,16,17,18,19,20,21,22,23,24,26,27];

GLMAnalPD{1} = '/imaging/es03/fMRI_2017/GLMAnalysisNative';
GLMAnalPD{2} = '/imaging/es03/fMRI_2017/GLMAnalysisNativeWithSyl1';
GLMAnalPD{3} = '/imaging/es03/fMRI_2017/GLMAnalysisNativeReduced';
GLMAnalPD{4} = '/imaging/es03/fMRI_2017/GLMAnalysisNativeReducedWithSyl1';

GLMNames = {'Syl2' 'Syl2+Syl1' 'Syl2Reduced' 'Syl2Reduced+Syl1'};

for k=1:length(SubjToAnalyze)
    for m=1:length(GLMAnalPD)
        V = spm_vol(fullfile(GLMAnalPD{m},Subj{SubjToAnalyze(k)},'MA_cvLME.nii'));
        img = spm_read_vols(V);
        F(k,m) = nanmean(img(:));
        
        V = spm_vol(fullfile(GLMAnalPD{m},Subj{SubjToAnalyze(k)},'MA_cvCom.nii'));
        img = spm_read_vols(V);
        com(k,m) = nanmean(img(:));
        
        V = spm_vol(fullfile(GLMAnalPD{m},Subj{SubjToAnalyze(k)},'MA_cvAcc.nii'));
        img = spm_read_vols(V);
        acc(k,m) = nanmean(img(:));
        
        clear SPM GLM
        load(fullfile(GLMAnalPD{m},Subj{SubjToAnalyze(k)},'SPM.mat'));
        GLM = MA_inspect_GoF('LoadData',SPM);
        R2(k,m) = nanmean(GLM.R2(:));
    end
end

F_group = sum(F,1);
F_group = F_group-min(F_group);

com_group = sum(com,1);
com_group = com_group-min(com_group);

acc_group = sum(acc,1);
acc_group = acc_group-min(acc_group);

figure;
subplot(3,1,1); bar(F_group); set(gca,'xticklabel',GLMNames);
subplot(3,1,2); bar(com_group); set(gca,'xticklabel',GLMNames);
subplot(3,1,3); bar(acc_group); set(gca,'xticklabel',GLMNames);