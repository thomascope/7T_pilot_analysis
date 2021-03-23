clear all
close all
clc

spm fmri

RsaAnalName = 'fMRI_2017';
RsaAnalPD = '/imaging/es03/fMRI_2017/GLMAnalysisNative/RSA_pearson_fromSimulations';
OutputPD = [RsaAnalPD '/SPMs/Anova_1way_Syl2_Masked'];

Subj = {'subj1' 'subj2' 'subj3' 'subj4' 'subj5' 'subj6' 'subj7' 'subj8' 'subj9' 'subj10' 'subj11' 'subj12' 'subj13' 'subj14' 'subj15' 'subj16' 'subj17' 'subj18' 'subj19' 'subj20' 'subj21' 'subj22' 'subj23' 'subj24' 'subj25' 'subj26' 'subj27'};
NumSessions = [5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5];

SubjToAnalyze = [3,4,5,6,8,10,11,12,13,15,16,17,18,19,20,21,22,23,24,26,27];
models = {'PE' 'PEabs' 'Pred' 'Sensory'};

%% Contrasts

cnt = 0;

% F contrasts
cnt = cnt+1;
contrasts{cnt}.type = 'F';
contrasts{cnt}.c = orth(diff(eye(4))')';
contrasts{cnt}.name = 'Main effect of model';

%% Estimate models
    
if ~exist(OutputPD)
    mkdir(OutputPD);
end

files = {};
for k=1:length(SubjToAnalyze)
    
    SubjCurrent=Subj{SubjToAnalyze(k)};
           
    for m=1:length(models)
        files{1}{k}{m} = [RsaAnalPD '/Maps/sw' RsaAnalName '_rMap_mask_' models{m} '_' SubjCurrent '.img'];
    end         
end

% set up input structure for batch_spm_anova_vES
S.imgfiles = files;
S.outdir = OutputPD;
mskname = '/imaging/es03/fMRI_2017/Masks/Bilateral_Temporal_Superior_AAL2.nii';
S.uUFp = .05;
if exist('mskname'); S.maskimg = mskname; end;
if exist('contrasts'); S.contrasts = contrasts; end;
if exist('covariates'); S.user_regs = covariates; end;

% estimate model and compute contrasts
batch_spm_anova_es(S);
   