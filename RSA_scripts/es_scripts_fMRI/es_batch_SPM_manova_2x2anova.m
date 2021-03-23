clear all
close all
clc

spm fmri

GLMAnalPD = '/imaging/es03/fMRI_2017/GLMAnalysisNativeMasked';
OutputPD = [GLMAnalPD '/Manova/SPMs/Anova_2x2'];

Subj = {'subj1' 'subj2' 'subj3' 'subj4' 'subj5' 'subj6' 'subj7' 'subj8' 'subj9' 'subj10' 'subj11' 'subj12' 'subj13' 'subj14' 'subj15' 'subj16' 'subj17' 'subj18' 'subj19' 'subj20' 'subj21' 'subj22' 'subj23' 'subj24' 'subj25' 'subj26' 'subj27'};
NumSessions = [5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5];

SubjToAnalyze = [3,4,5,6,8,10,11,12,13,15,16,17,18,19,20,21,22,23,24,26,27];
%con2analyse = [7:10];
con2analyse = {'Syl2StrongM' 'Syl2WeakM' 'Syl2StrongMM' 'Syl2WeakMM'};

%% Contrasts

cnt = 0;

% F contrasts
cnt = cnt+1;
contrasts{cnt}.type = 'F';
contrasts{cnt}.c = kron([1/2 1/2],[1 -1]);
contrasts{cnt}.name = 'Main effect of prior strength';
cnt = cnt+1;
contrasts{cnt}.type = 'F';
contrasts{cnt}.c = kron([1 -1],[1/2 1/2]);
contrasts{cnt}.name = 'Main effect of prior congruency';
cnt = cnt+1;
contrasts{cnt}.type = 'F';
contrasts{cnt}.c = kron([1 -1],[1 -1]);
contrasts{cnt}.name = 'Prior strength X congruency';

%% Estimate models
    
if ~exist(OutputPD)
    mkdir(OutputPD);
end

files = {};
for k=1:length(SubjToAnalyze)
    
    SubjCurrent=Subj{SubjToAnalyze(k)};
           
    for c=1:length(con2analyse)
        %files{1}{k}{c} = sprintf([GLMAnalPD '/' SubjCurrent '/' 'swspmDs_C%04d_P0001.nii'],con2analyse(c));
        files{1}{k}{c} = sprintf([GLMAnalPD '/' SubjCurrent '/' 'sweffect-map_%s.nii'],con2analyse{c});
    end         
end

% set up input structure for batch_spm_anova_vES
S.imgfiles = files;
S.outdir = OutputPD;
%mskname = [];
if exist('mskname'); S.maskimg = mskname; end;
if exist('contrasts'); S.contrasts = contrasts; end;
if exist('covariates'); S.user_regs = covariates; end;
S.uUFp = 1;

% estimate model and compute contrasts
batch_spm_anova_es(S);
   