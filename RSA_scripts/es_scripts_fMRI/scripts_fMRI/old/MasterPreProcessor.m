clear all
clc

spm fmri

Subj = {'subj1' 'subj2' 'subj3' 'subj4' 'subj5' 'subj6' 'subj7' 'subj8' 'subj9' 'subj10' 'subj11' 'subj12' 'subj13' 'subj14' 'subj15' 'subj16' 'subj17' 'subj18' 'subj19' 'subj20' 'subj21' 'subj22'};
NumSessions=[5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5];
StartDir=pwd;
DataPD='/imaging/es03/fMRI_2017/PreprocessAnalysis'; % Parent Directory where raw functional data is saved
TempPD='/group/language/data/ediz.sohoglu/projects/fMRI_2017/scripts_fMRI/Templates'; % Template PD
SpmPD = '/imaging/local/software/spm_cbu_svn/releases/spm12_latest';

SubjToAnalyze=[22];

for k=1:length(SubjToAnalyze)
    
    SubjCurrent=Subj{SubjToAnalyze(k)};
    NumSess=NumSessions(SubjToAnalyze(k));
    
    disp(['Preprocessing for Subject:' SubjCurrent]);
    
    SubjDir =fullfile(DataPD,SubjCurrent);
    for sess=1:NumSess
        FuncDir{sess}=fullfile(SubjDir,'Functional',['Sess_' num2str(sess)]);
    end
    
    StrDir =fullfile(SubjDir, 'Structural');
    
    PreProcessData(SubjDir,FuncDir,StrDir,TempPD,SpmPD);
    
    clear SubjDir FuncDir StrDir
end