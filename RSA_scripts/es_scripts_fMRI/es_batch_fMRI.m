clearvars
close all
clc

PreProcPD = '/imaging/es03/fMRI_2017/PreprocessAnalysis';
%GLMAnalPD = '/imaging/es03/fMRI_2017/GLMAnalysisMNISmooth6mm';
GLMAnalPD = '/imaging/es03/fMRI_2017/GLMAnalysisNative';
%GLMAnalPD = '/imaging/es03/fMRI_2017/GLMAnalysisMNISmooth6mmByItem';
BehavDataPD ='/group/language/data/ediz.sohoglu/projects/fMRI_2017/analysis_behavioural/pauseDetection';
SpmPD = '/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7219';
TempPD = '/group/language/data/ediz.sohoglu/projects/fMRI_2017/scripts_fMRI/Templates'; % Template PD
RsaToolboxPD = '/group/language/data/ediz.sohoglu/matlab/rsatoolbox';
ScriptsPD = '/group/language/data/ediz.sohoglu/projects/fMRI_2017/scripts_fMRI';

Subj = {'subj1' 'subj2' 'subj3' 'subj4' 'subj5' 'subj6' 'subj7' 'subj8' 'subj9' 'subj10' 'subj11' 'subj12' 'subj13' 'subj14' 'subj15' 'subj16' 'subj17' 'subj18' 'subj19' 'subj20' 'subj21' 'subj22' 'subj23' 'subj24' 'subj25' 'subj26' 'subj27'};
NumSessions = [5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5];

%SubjToAnalyze = [3,4,5,6,8,10,11,12,13,15,16,17,18,19,20,21,22,23,24,26,27];
SubjToAnalyze = [3];

%% Parallel computing settings

addpath /hpc-software/matlab/cbu/

S = cbu_scheduler();
S.NumWorkers = length(SubjToAnalyze);
S.SubmitArguments = '-l mem=4GB -l walltime=0:30:00';

%% Preprocess

clear J
for k=1:length(SubjToAnalyze)
    
    SubjCurrent=Subj{SubjToAnalyze(k)};
    NumSess=NumSessions(SubjToAnalyze(k));
    
    disp(['Preprocessing for Subject:' SubjCurrent]);
    
    SubjDir =fullfile(PreProcPD,SubjCurrent);
    for sess=1:NumSess
        FuncDir{sess} = fullfile(SubjDir,'Functional',['Sess_' num2str(sess)]);
    end
    
    StrDir = fullfile(SubjDir, 'Structural');
    
    PreProcessData_1Subj(SubjDir,FuncDir,StrDir,TempPD,SpmPD,ScriptsPD);
    J(k).task = @PreProcessData_1Subj; % External function name here
    J(k).n_return_values = 0; % important
    J(k).input_args = {SubjDir,FuncDir,StrDir,TempPD,SpmPD,ScriptsPD};
    J(k).depends_on = 0;
        
end

%cbu_qsub(J, S);

%% Create regressors

for k=1:length(SubjToAnalyze)
    
    SubjCurrent=Subj{SubjToAnalyze(k)};
    NumSess=NumSessions(SubjToAnalyze(k));
    
    disp(['Creating regressors for Subject:' SubjCurrent]);
    
    GLMAnalDir =[GLMAnalPD '/' SubjCurrent];
    if ~isdir(GLMAnalDir)
        mkdir(GLMAnalDir);
    end
    
    clear BehavDataDir PhysioDataFile MotionDataDir
    for m=1:NumSess
        BehavDataDir{m} = [BehavDataPD '/'];
        %PhysioDataFile{m} =[BehavDataPD '/' SubjCurrent '/' SubjCurrent '_R_session' num2str(m)];
        PhysioDataFile{m} = [];
        MotionDataDir{m} = [PreProcPD '/' SubjCurrent '/Functional/Sess_' num2str(m)];   
    end
    
    %CreateRegressors(GLMAnalDir,BehavDataDir,PhysioDataFile,MotionDataDir,SubjCurrent);
    CreateRegressorsForRSA(GLMAnalDir,BehavDataDir,PhysioDataFile,MotionDataDir,SubjCurrent);
    %CreateRegressorsForRSAWithSyl1(GLMAnalDir,BehavDataDir,PhysioDataFile,MotionDataDir,SubjCurrent);
    %CreateRegressorsForRSAReduced(GLMAnalDir,BehavDataDir,PhysioDataFile,MotionDataDir,SubjCurrent);
    %CreateRegressorsForRSAReducedBySyl1(GLMAnalDir,BehavDataDir,PhysioDataFile,MotionDataDir,SubjCurrent);
    %CreateRegressorsForRSAReducedWithSyl1(GLMAnalDir,BehavDataDir,PhysioDataFile,MotionDataDir,SubjCurrent);
    %CreateRegressorsForUnivariateCorrelations3(GLMAnalDir,BehavDataDir,PhysioDataFile,MotionDataDir,SubjCurrent);
    
end

%% Run first-level GLMs in parallel

clear J
for k=1:length(SubjToAnalyze)
    
    SubjCurrent=Subj{SubjToAnalyze(k)};
    NumSess=NumSessions(SubjToAnalyze(k));

    disp(['Running GLM for Subject:' SubjCurrent]);
    
    GLMAnalDir =[GLMAnalPD '/' SubjCurrent];
    ProcDataDir =[PreProcPD '/' SubjCurrent];    
    
    GLMAnalysisForRSA_1Subj(ProcDataDir,GLMAnalDir,TempPD,SubjCurrent,NumSess,SpmPD,ScriptsPD)
    %J(k).task = @GLMAnalysis_1Subj; % External function name here
    %J(k).task = @GLMAnalysisForRSA_1Subj; % External function name here
    %J(k).task = @GLMAnalysisForRSAMasked_1Subj; % External function name here
    J(k).task = @GLMAnalysisByItem_1Subj; % External function name here
    J(k).n_return_values = 0; % important
    J(k).input_args = {ProcDataDir,GLMAnalDir,TempPD,SubjCurrent,NumSess,SpmPD,ScriptsPD};
    J(k).depends_on = 0;
    
end

%cbu_qsub(J, S);

%% Run RSA in parallel

clear J
for k=1:length(SubjToAnalyze)
    
    SubjCurrent=Subj{SubjToAnalyze(k)};

    disp(['Running RSA for Subject:' SubjCurrent]);
    
    
    %RSASearchlightAnalysis_1Subj(PreProcPD,SubjCurrent,RsaToolboxPD,TempPD,SpmPD,ScriptsPD);
    J(k).task = @RSASearchlightAnalysis_1Subj; % External function name here
    J(k).n_return_values = 0; % important
    J(k).input_args = {PreProcPD,SubjCurrent,RsaToolboxPD,TempPD,SpmPD,ScriptsPD};
    J(k).depends_on = 0;
    
end

cbu_qsub(J, S);

%% T-test

spm fmri

RsaAnalName = 'fMRI_2017';
RsaAnalPD = '/imaging/es03/fMRI_2017/GLMAnalysisNative/RSA_pearson_fromSimulations';
%explicitMask = [];
explicitMask = '/imaging/es03/fMRI_2017/Masks/Bilateral_Temporal_Superior_AAL2.nii';

clear images models
for k=1:length(SubjToAnalyze)
    
    SubjCurrent = Subj{SubjToAnalyze(k)};
    
    imagesCurrent = dir([RsaAnalPD '/Maps/' RsaAnalName '_rMap_mask_*_' SubjCurrent '.img']);
    
    for m=1:length(imagesCurrent)
        
        S = strsplit(imagesCurrent(m).name,'_');
        models{m} = strtok(S{5},'.');
        
        images{m}{k,1} = [RsaAnalPD '/Maps/sw' RsaAnalName '_rMap_mask_' models{m} '_' SubjCurrent '.img'];
        
    end
    
end

for m=1:length(models)
    
    dirOutput = [RsaAnalPD '/SPMs/Ttest_Masked_' models{m}];
    if exist(dirOutput,'dir');
        rmdir(dirOutput,'s');
    end
    mkdir(dirOutput);
       
    clear matlabbatch
    load([TempPD '/SimpleTtest.mat']);
    matlabbatch{1}.spm.stats.factorial_design.dir = {dirOutput};
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = images{m};
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {explicitMask};
    save(fullfile(dirOutput,'SimpleTtest.mat'), 'matlabbatch');
    spm_jobman('initcfg')
    spm_jobman('run', matlabbatch);
    
    % Estimate SPM
    cd(dirOutput);
    load([dirOutput '/SPM.mat']);
    spm_spm(SPM);
    
    % Add contrasts
    clear matlabbatch
    matlabbatch{1}.spm.stats.con.spmmat = {[dirOutput '/SPM.mat']};
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = '+ve';
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = '-ve';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = -1;
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.delete = 0;
    save(fullfile(dirOutput,'ContrastManager.mat'), 'matlabbatch');
    spm_jobman('initcfg')
    spm_jobman('run', matlabbatch);
    
end

    