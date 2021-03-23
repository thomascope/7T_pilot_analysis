clear all
close all
clc

spm fmri

GLMAnalPD = '/imaging/es03/fMRI_2017/GLMAnalysisMNISmooth6mm';
OutputPD = [GLMAnalPD '/SPMs/Anova_2x3_Masked'];

Subj = {'subj1' 'subj2' 'subj3' 'subj4' 'subj5' 'subj6' 'subj7' 'subj8' 'subj9' 'subj10' 'subj11' 'subj12' 'subj13' 'subj14' 'subj15' 'subj16' 'subj17' 'subj18' 'subj19' 'subj20' 'subj21' 'subj22' 'subj23' 'subj24' 'subj25' 'subj26' 'subj27'};
NumSessions = [5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5];

SubjToAnalyze = [3,4,5,6,8,10,11,12,13,15,16,17,18,19,20,21,22,23,24,26,27];
con2analyse = [9 10 7 5 6 8];
% Order of conditions in SPM: {'Strong+M','Weak+M','Strong+MM','Weak+MM','Strong+Noise','Weak+Noise','Noise+Speech','Noise+Noise','Strong+Clear','Weak+Clear','HitsAndFalseAlarms','Misses'}
% To analyse: {'Strong+Clear','Weak+Clear','Noise+Clear','Strong+Noise' 'Weak+Noise' 'Noise+Noise'}

%% Contrasts

cnt = 0;

% F contrasts
cnt = cnt+1;
contrasts{cnt}.type = 'F';
contrasts{cnt}.c = kron([1/2 1/2],orth(diff(eye(3))')');
contrasts{cnt}.name = 'Main effect of prior strength';
cnt = cnt+1;
contrasts{cnt}.type = 'F';
contrasts{cnt}.c = kron([1 -1],[1/3 1/3 1/3]);
contrasts{cnt}.name = 'Main effect of second syllable type';
cnt = cnt+1;
contrasts{cnt}.type = 'F';
contrasts{cnt}.c = kron([1 -1],orth(diff(eye(3))')');
contrasts{cnt}.name = 'Prior strength X second sylable type';

%% Estimate models
    
if ~exist(OutputPD)
    mkdir(OutputPD);
end

files = {};
for k=1:length(SubjToAnalyze)
    
    SubjCurrent=Subj{SubjToAnalyze(k)};
           
    for c=1:length(con2analyse)
        files{1}{k}{c} = sprintf([GLMAnalPD '/' SubjCurrent '/' 'con_%04d.nii'],con2analyse(c));
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
   