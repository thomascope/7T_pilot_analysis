clearvars
clc

GLMAnalPD = '/imaging/es03/fMRI_2017/GLMAnalysisMNISmooth6mm';
BehavDataPD ='/group/language/data/ediz.sohoglu/projects/fMRI_2017/analysis_behavioural/pauseDetection';
FuncDataPD ='/imaging/es03/fMRI_2017/PreprocessAnalysis';

Subj = {'subj1' 'subj2' 'subj3' 'subj4' 'subj5' 'subj6' 'subj7' 'subj8' 'subj9' 'subj10' 'subj11' 'subj12' 'subj13' 'subj14' 'subj15' 'subj16' 'subj17' 'subj18' 'subj19' 'subj20' 'subj21' 'subj22' 'subj23' 'subj24' 'subj25' 'subj26' 'subj27'};
NumSessions = [5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5];

SubjToAnalyze = [3,4,5,6,8,10,11,12,13,15,16,17,18,19,20,21,22,23,24,26,27];

for k=1:length(SubjToAnalyze);
    
    SubjCurrent=Subj{SubjToAnalyze(k)};
    NumSess=NumSessions(SubjToAnalyze(k));
    
    GLMAnalDir =[GLMAnalPD '/' SubjCurrent];
    if ~isdir(GLMAnalDir)
        mkdir(GLMAnalDir);
    end
    
    clear BehavDataDir PhysioDataFile MotionDataDir
    for m=1:NumSess
        BehavDataDir{m} = [BehavDataPD '/'];
        %PhysioDataFile{m} =[BehavDataPD '/' SubjCurrent '/' SubjCurrent '_R_session' num2str(m)];
        PhysioDataFile{m} = [];
        MotionDataDir{m} = [FuncDataPD '/' SubjCurrent '/Functional/Sess_' num2str(m)];   
    end
    
    CreateRegressorsForRSA(GLMAnalDir,BehavDataDir,PhysioDataFile,MotionDataDir,SubjCurrent);
    
end
