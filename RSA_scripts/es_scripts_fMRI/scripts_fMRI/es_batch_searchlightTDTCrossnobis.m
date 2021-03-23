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
    
    TDTCrossnobisAnalysis_1Subj(SubjCurrent,GLMAnalPD,PreProcPD,TempPD);
    J(k).task = @TDTCrossnobisAnalysis_1Subj; % External function name here
    J(k).n_return_values = 0; % important
    J(k).input_args = {SubjCurrent,GLMAnalPD,PreProcPD,TempPD};
    J(k).depends_on = 0;
    
end

%cbu_qsub(J, S);

%% T-test

version = 'spearmanPartial';

%explicitMask = [];
explicitMask = '/imaging/es03/fMRI_2017/Masks/Bilateral_Temporal_Superior_AAL2.nii';

SubjFirst = Subj{SubjToAnalyze(1)};
imagesSubjFirst = dir([GLMAnalPD '/' SubjFirst '/TDTcrossnobis/' version '/sweffect-map_*.nii']);

clear images
for k=1:length(SubjToAnalyze)
    
    SubjCurrent = Subj{SubjToAnalyze(k)};
        
    for m=1:length(imagesSubjFirst)
        
        S = strsplit(imagesSubjFirst(m).name,'_');
        models{m} = strtok(S{2},'.');
        
        images{m}{k,1} = [GLMAnalPD '/' SubjCurrent '/TDTcrossnobis/' version '/sweffect-map_' models{m} '.nii'];
        
    end
    
end

spm fmri

for m=1:length(models)
    
    dirOutput = [GLMAnalPD '/TDTcrossnobis_' version '/SPMs/Ttest_Masked_' models{m}];
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
    try
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
        matlabbatch{1}.spm.stats.con.consess{3}.fcon.name = '+/-ve';
        matlabbatch{1}.spm.stats.con.consess{3}.fcon.weights = 1;
        matlabbatch{1}.spm.stats.con.consess{3}.fcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.delete = 0;
        save(fullfile(dirOutput,'ContrastManager.mat'), 'matlabbatch');
        spm_jobman('initcfg')
        spm_jobman('run', matlabbatch);
    end
    
end

%% T-test (paired)

version = 'spearman';

%explicitMask = [];
explicitMask = '/imaging/es03/fMRI_2017/Masks/Bilateral_Temporal_Superior_AAL2.nii';

models = {
    {'Syl2StrongM' 'Syl2WeakM'}
    };

clear images
for k=1:length(SubjToAnalyze)
    
    SubjCurrent = Subj{SubjToAnalyze(k)};
        
    for m=1:length(models)
        
        images{m}{k,1} = [GLMAnalPD '/' SubjCurrent '/TDTcrossnobis/' version '/sweffect-map_' models{m}{1} '.nii'];
        images{m}{k,2} = [GLMAnalPD '/' SubjCurrent '/TDTcrossnobis/' version '/sweffect-map_' models{m}{2} '.nii'];
        
    end
    
end

spm fmri

for m=1:length(models)
    
    dirOutput = [GLMAnalPD '/TDTcrossnobis_' version '/SPMs/TtestPaired_Masked_' models{m}{1} 'vs' models{m}{2}];
    if exist(dirOutput,'dir');
        rmdir(dirOutput,'s');
    end
    mkdir(dirOutput);
       
    clear matlabbatch
    load([TempPD '/SimpleTtestPaired.mat']);
    matlabbatch{1}.spm.stats.factorial_design.dir = {dirOutput};
    for k=1:length(SubjToAnalyze)
        matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(k).scans = images{m}(k,:)';
    end
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {explicitMask};
    save(fullfile(dirOutput,'SimpleTtestPaired.mat'), 'matlabbatch');
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
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = '-ve';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.delete = 0;
    save(fullfile(dirOutput,'ContrastManager.mat'), 'matlabbatch');
    spm_jobman('initcfg')
    spm_jobman('run', matlabbatch);
    
end

%% T-test (after computing averages/subtractions of effect images)

%explicitMask = [];
explicitMask = '/imaging/es03/fMRI_2017/Masks/Bilateral_Temporal_Superior_AAL2.nii';

version = 'spearmanPartial';
models = {'EntropyM' 'EntropyMM'};

clear images
for k=1:length(SubjToAnalyze)
    
    SubjCurrent = Subj{SubjToAnalyze(k)};
    
    for m=1:length(models)       
        images{m}{k,1} = [GLMAnalPD '/' SubjCurrent '/TDTcrossnobis/' version '/sweffect-map_' models{m} '.nii'];       
    end
    
    V = spm_vol(images{1}{k,1});
    img1 = spm_read_vols(V);
    V = spm_vol(images{2}{k,1});
    img2 = spm_read_vols(V);
    
    img_new = img2-img1;
    V.fname = [GLMAnalPD '/' SubjCurrent '/TDTcrossnobis/' version '/sweffect-map_EntropyMM-M.nii'];
    spm_write_vol(V,img_new);
    images{3}{k,1} = V.fname;
    
    img_new = (img2+img1)/2;
    V.fname = [GLMAnalPD '/' SubjCurrent '/TDTcrossnobis/' version '/sweffect-map_EntropyMM+M.nii'];
    spm_write_vol(V,img_new);
    images{4}{k,1} = V.fname;
               
end

models = {'EntropyM' 'EntropyMM' 'EntropyMM-M' 'EntropyMM+M'};

%spm fmri

for m=1:length(models)
    
    dirOutput = [GLMAnalPD '/TDTcrossnobis/' version '/SPMs/Ttest_Masked_' models{m}];
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
    try
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
        matlabbatch{1}.spm.stats.con.consess{3}.fcon.name = '+/-ve';
        matlabbatch{1}.spm.stats.con.consess{3}.fcon.weights = 1;
        matlabbatch{1}.spm.stats.con.consess{3}.fcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.delete = 0;
        save(fullfile(dirOutput,'ContrastManager.mat'), 'matlabbatch');
        spm_jobman('initcfg')
        spm_jobman('run', matlabbatch);
    end
    
end

    