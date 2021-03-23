clear all
close all
clc

spm fmri

GLMAnalPD = '/imaging/es03/fMRI_2017/GLMAnalysisNative';
%OutputPD = [GLMAnalPD '/TDTcrossnobis_CompareRegressionModelsWithinMatchingOnly2/SPMs/Anova_Masked'];
OutputPD = [GLMAnalPD '/TDTcrossnobis_CompareRegressionModelsFromSim/SPMs/Anova_Masked'];

clear version

cnt = 0;

% cnt = cnt + 1;
% version{cnt} = 'regressionSyl2GradedGatingWithinMatchingOnly';
% cnt = cnt + 1;
% version{cnt} = 'regressionSyl2GradedGatingfMRIWithinMatchingOnly';
% cnt = cnt + 1;
% version{cnt} = 'regressionSyl2GradedProbWithinMatchingOnly';
% cnt = cnt + 1;
% version{cnt} = 'regressionSyl2GradedProb2WithinMatchingOnly';
% 
% % cnt = cnt + 1;
% % version{cnt} = 'regressionSyl2GradedNormGatingWithinMatchingOnly';
% % cnt = cnt + 1;
% % version{cnt} = 'regressionSyl2GradedNormGatingfMRIWithinMatchingOnly';
% % cnt = cnt + 1;
% % version{cnt} = 'regressionSyl2GradedNormProbWithinMatchingOnly';
% % cnt = cnt + 1;
% % version{cnt} = 'regressionSyl2GradedNormProb2WithinMatchingOnly';
% 
% cnt = cnt + 1;
% version{cnt} = 'regressionSyl2GatingWithinMatchingOnly';
% cnt = cnt + 1;
% version{cnt} = 'regressionSyl2GatingfMRIWithinMatchingOnly';
% cnt = cnt + 1;
% version{cnt} = 'regressionSyl2ProbWithinMatchingOnly';
% cnt = cnt + 1;
% version{cnt} = 'regressionSyl2Prob2WithinMatchingOnly';
% 
% % cnt = cnt + 1;
% % version{cnt} = 'regressionSyl2GradedGatingMatchingOnly';
% % cnt = cnt + 1;
% % version{cnt} = 'regressionSyl2GradedGatingfMRIMatchingOnly';
% % cnt = cnt + 1;
% % version{cnt} = 'regressionSyl2GradedProbMatchingOnly';
% % cnt = cnt + 1;
% % version{cnt} = 'regressionSyl2GradedProb2MatchingOnly';
% % 
% % cnt = cnt + 1;
% % version{cnt} = 'regressionSyl2GradedNormGatingMatchingOnly';
% % cnt = cnt + 1;
% % version{cnt} = 'regressionSyl2GradedNormGatingfMRIMatchingOnly';
% % cnt = cnt + 1;
% % version{cnt} = 'regressionSyl2GradedNormProbMatchingOnly';
% % cnt = cnt + 1;
% % version{cnt} = 'regressionSyl2GradedNormProb2MatchingOnly';
% % 
% % cnt = cnt + 1;
% % version{cnt} = 'regressionSyl2GatingMatchingOnly';
% % cnt = cnt + 1;
% % version{cnt} = 'regressionSyl2GatingfMRIMatchingOnly';
% % cnt = cnt + 1;
% % version{cnt} = 'regressionSyl2ProbMatchingOnly';
% % cnt = cnt + 1;
% % version{cnt} = 'regressionSyl2Prob2MatchingOnly';

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

Subj = {'subj1' 'subj2' 'subj3' 'subj4' 'subj5' 'subj6' 'subj7' 'subj8' 'subj9' 'subj10' 'subj11' 'subj12' 'subj13' 'subj14' 'subj15' 'subj16' 'subj17' 'subj18' 'subj19' 'subj20' 'subj21' 'subj22' 'subj23' 'subj24' 'subj25' 'subj26' 'subj27'};
NumSessions = [5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5];

SubjToAnalyze = [3,4,5,6,8,10,11,12,13,15,16,17,18,19,20,21,22,23,24,26,27];

%% Contrasts

% cnt = 0;
% 
% % F contrasts
% cnt = cnt+1;
% contrasts{cnt}.type = 'F';
% contrasts{cnt}.c = kron([1/2 1/2],orth(diff(eye(4))')');
% contrasts{cnt}.name = 'Main effect of prior source';
% cnt = cnt+1;
% contrasts{cnt}.type = 'F';
% contrasts{cnt}.c = kron(orth(diff(eye(2))')',[1/4 1/4 1/4 1/4]);
% contrasts{cnt}.name = 'Main effect of syllable 2 type';
% cnt = cnt+1;
% contrasts{cnt}.type = 'F';
% contrasts{cnt}.c = kron(orth(diff(eye(2))')',orth(diff(eye(4))')');
% contrasts{cnt}.name = 'Interaction';

%% Contrasts

cnt = 0;

% F contrasts
cnt = cnt+1;
contrasts{cnt}.type = 'F';
contrasts{cnt}.c = kron([1/3 1/3 1/3],orth(diff(eye(3))')');
contrasts{cnt}.name = 'Main effect of PE/input/pred';
cnt = cnt+1;
contrasts{cnt}.type = 'F';
contrasts{cnt}.c = kron(orth(diff(eye(3))')',[1/3 1/3 1/3]);
contrasts{cnt}.name = 'Main effect of rep type';
cnt = cnt+1;
contrasts{cnt}.type = 'F';
contrasts{cnt}.c = kron(orth(diff(eye(3))')',orth(diff(eye(3))')');
contrasts{cnt}.name = 'Interaction';

%% Estimate models
    
if ~exist(OutputPD)
    mkdir(OutputPD);
end

files = {};
clear imageVolumes; imageCounter = 0;
for k=1:length(SubjToAnalyze)
    
    SubjCurrent=Subj{SubjToAnalyze(k)};
           
    for v=1:length(version)
        files{1}{k}{v} = sprintf([GLMAnalPD '/' SubjCurrent '/TDTcrossnobis/' version{v} '/' 'sweffect-map_%s.nii'],'R');
        imageCounter = imageCounter + 1;
        imageVolumes(imageCounter) = spm_vol(files{1}{k}{v});
    end
    
end

% set up input structure for batch_spm_anova_vES
S.imgfiles = files;
S.outdir = OutputPD;
mskname = '/imaging/es03/fMRI_2017/Masks/Bilateral_Temporal_Superior_AAL2.nii';
S.uUFp = .005;
if exist('mskname'); S.maskimg = mskname; end;
if exist('contrasts'); S.contrasts = contrasts; end;
if exist('covariates'); S.user_regs = covariates; end;

% estimate model and compute contrasts
batch_spm_anova_es(S);

% compute average over all images
Vo = spm_imcalc(imageVolumes,fullfile(OutputPD,'GrandAverage.nii'),'mean(X,1)',{1,0,0,'float64','GrandAverageR'});

   