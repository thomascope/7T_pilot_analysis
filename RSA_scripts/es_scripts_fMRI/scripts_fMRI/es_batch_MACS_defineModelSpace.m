clearvars

spm fmri

Subj = {'subj1' 'subj2' 'subj3' 'subj4' 'subj5' 'subj6' 'subj7' 'subj8' 'subj9' 'subj10' 'subj11' 'subj12' 'subj13' 'subj14' 'subj15' 'subj16' 'subj17' 'subj18' 'subj19' 'subj20' 'subj21' 'subj22' 'subj23' 'subj24' 'subj25' 'subj26' 'subj27'};
SubjToAnalyze = [3,4,5,6,8,10,11,12,13,15,16,17,18,19,20,21,22,23,24,26,27];

dirOutput = '/imaging/es03/fMRI_2017/MACS_ModelSelection/ForRSA';
GLMAnalPD{1} = '/imaging/es03/fMRI_2017/GLMAnalysisNative';
GLMAnalPD{2} = '/imaging/es03/fMRI_2017/GLMAnalysisNativeWithSyl1';
GLMAnalPD{3} = '/imaging/es03/fMRI_2017/GLMAnalysisNativeReduced';
GLMAnalPD{4} = '/imaging/es03/fMRI_2017/GLMAnalysisNativeReducedWithSyl1';

GLMNames = {'Syl2' 'Syl2+Syl1' 'Syl2Reduced' 'Syl2Reduced+Syl1'};

if exist(dirOutput,'dir');
    rmdir(dirOutput,'s');
end
mkdir(dirOutput);

models = {};
for k=1:length(SubjToAnalyze)
    for m=1:length(GLMAnalPD)
        models{k}{m} = {fullfile(GLMAnalPD{m},Subj{SubjToAnalyze(k)},'SPM.mat')};
    end
end

clear matlabbatch
matlabbatch{1}.spm.tools.MACS.MA_model_space.dir = {dirOutput};
matlabbatch{1}.spm.tools.MACS.MA_model_space.models = models;
matlabbatch{1}.spm.tools.MACS.MA_model_space.names = GLMNames;

save(fullfile(dirOutput,'MA_defineModelSpace.mat'), 'matlabbatch');
spm_jobman('initcfg')
spm_jobman('run', matlabbatch);