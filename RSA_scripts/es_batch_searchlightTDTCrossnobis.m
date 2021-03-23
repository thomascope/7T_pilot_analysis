clearvars

Subj = {'subj1' 'subj2' 'subj3' 'subj4' 'subj5' 'subj6' 'subj7' 'subj8' 'subj9' 'subj10' 'subj11' 'subj12' 'subj13' 'subj14' 'subj15' 'subj16' 'subj17' 'subj18' 'subj19' 'subj20' 'subj21' 'subj22' 'subj23' 'subj24' 'subj25' 'subj26' 'subj27'};
SubjToAnalyze = [3,4,5,6,8,10,11,12,13,15,16,17,18,19,20,21,22,23,24,26,27];
%SubjToAnalyze = 3;
PreProcPD = '/imaging/es03/fMRI_2017/PreprocessAnalysis';
GLMAnalPD = '/imaging/es03/fMRI_2017/GLMAnalysisNative';
TempPD = '/group/language/data/ediz.sohoglu/projects/fMRI_2017/scripts_fMRI/Templates'; % Template PD

%% Parallel computing settings

addpath /hpc-software/matlab/cbu/

S = cbu_scheduler();
S.NumWorkers = length(SubjToAnalyze);
S.SubmitArguments = '-l mem=16GB -l walltime=10:00:00';

%% TDT specific set-up

clear J
for k=1:length(SubjToAnalyze)
    
    SubjCurrent = Subj{SubjToAnalyze(k)};
    
    disp(['Searchlight for Subject:' SubjCurrent]);
    
    %TDTCrossnobisAnalysis_1Subj(SubjCurrent,GLMAnalPD,PreProcPD,TempPD);
    J(k).task = @TDTCrossnobisAnalysis_1Subj; % External function name here
    J(k).n_return_values = 0; % important
    J(k).input_args = {SubjCurrent,GLMAnalPD,PreProcPD,TempPD};
    J(k).depends_on = 0;
    
end

cbu_qsub(J, S);
