clearvars
close all
clc

spm fmri

GLMAnalPD = '/imaging/es03/fMRI_2017/GLMAnalysisForRSA2';
PreProcPD = '/imaging/es03/fMRI_2017/PreprocessAnalysis';
TempPD = '/group/language/data/ediz.sohoglu/projects/fMRI_2017/scripts_fMRI/Templates'; % Template PD 

Subj = {'subj1' 'subj2' 'subj3' 'subj4' 'subj5' 'subj6' 'subj7' 'subj8' 'subj9' 'subj10' 'subj11' 'subj12' 'subj13' 'subj14' 'subj15' 'subj16' 'subj17' 'subj18' 'subj19' 'subj20' 'subj21'};
NumSessions=[5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5];

SubjToAnalyze = [10];

for k=1:length(SubjToAnalyze)
    
    SubjCurrent=Subj{SubjToAnalyze(k)};
    NumSess=NumSessions(SubjToAnalyze(k));

    disp(['Preprocessing for Subject:' SubjCurrent]);
    
    GLMAnalDir =[GLMAnalPD '/' SubjCurrent];
    ProcDataDir =[PreProcPD '/' SubjCurrent];    
    
    GLMAnalysisForRSA_1Subj(ProcDataDir,GLMAnalDir,TempPD,SubjCurrent,NumSess)
end