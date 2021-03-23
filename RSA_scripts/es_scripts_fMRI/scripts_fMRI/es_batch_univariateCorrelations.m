clearvars

Subj = {'subj1' 'subj2' 'subj3' 'subj4' 'subj5' 'subj6' 'subj7' 'subj8' 'subj9' 'subj10' 'subj11' 'subj12' 'subj13' 'subj14' 'subj15' 'subj16' 'subj17' 'subj18' 'subj19' 'subj20' 'subj21' 'subj22' 'subj23' 'subj24' 'subj25' 'subj26' 'subj27'};
SubjToAnalyze = [3,4,5,6,8,10,11,12,13,15,16,17,18,19,20,21,22,23,24,26,27];
%SubjToAnalyze = [26,27];
PreProcPD = '/imaging/es03/fMRI_2017/PreprocessAnalysis';
GLMAnalPD = '/imaging/es03/fMRI_2017/GLMAnalysisMNISmooth6mmByItem';
TempPD = '/group/language/data/ediz.sohoglu/projects/fMRI_2017/scripts_fMRI/Templates'; % Template PD

%% Parallel computing settings

% addpath /hpc-software/matlab/cbu/
% 
% S = cbu_scheduler();
% S.NumWorkers = length(SubjToAnalyze);
% S.SubmitArguments = '-l mem=16GB -l walltime=10:00:00';

%% Compute correlations

clear J
for k=1:length(SubjToAnalyze)
    
    SubjCurrent = Subj{SubjToAnalyze(k)};
    
    disp(['Correlations for Subject:' SubjCurrent]);
    
    UnivariateCorrelations_1Subj(SubjCurrent,GLMAnalPD,PreProcPD,TempPD);
    J(k).task = @UnivariateCorrelations_1Subj; % External function name here
    J(k).n_return_values = 0; % important
    J(k).input_args = {SubjCurrent,GLMAnalPD,PreProcPD,TempPD};
    J(k).depends_on = 0;
    
end

%cbu_qsub(J, S);

% %% T-test
% 
% %explicitMask = [];
% explicitMask = '/imaging/es03/fMRI_2017/Masks/Bilateral_Temporal_Superior_AAL2.nii';
% 
% version = 'univariateCorrelations';
% models = {'Entropy_M' 'Entropy_MM' 'Entropy_Noise'};
% 
% clear images
% for k=1:length(SubjToAnalyze)
%     
%     SubjCurrent = Subj{SubjToAnalyze(k)};
%     
%     for m=1:length(models)       
%         images{m}{k,1} = [GLMAnalPD '/' SubjCurrent '/' version '/effect-map_' models{m} '.nii'];       
%     end
%     
%     V = spm_vol(images{1}{k,1});
%     img1 = spm_read_vols(V);
%     V = spm_vol(images{2}{k,1});
%     img2 = spm_read_vols(V);
%     
%     img_new = img2-img1;
%     V.fname = [GLMAnalPD '/' SubjCurrent '/' version '/effect-map_Entropy_MM-M.nii'];
%     spm_write_vol(V,img_new);
%     images{4}{k,1} = V.fname;
%     
%     img_new = (img2+img1)/2;
%     V.fname = [GLMAnalPD '/' SubjCurrent '/' version '/effect-map_Entropy_MM+M.nii'];
%     spm_write_vol(V,img_new);
%     images{5}{k,1} = V.fname;
%                
% end
% 
% models = {'Entropy_M' 'Entropy_MM' 'Entropy_Noise' 'Entropy_MM-M' 'Entropy_MM+M'};
% 
% %spm fmri
% 
% for m=1:length(models)
%     
%     dirOutput = [GLMAnalPD '/' version '/SPMs/Ttest_Masked_' models{m}];
%     if exist(dirOutput,'dir');
%         rmdir(dirOutput,'s');
%     end
%     mkdir(dirOutput);
%        
%     clear matlabbatch
%     load([TempPD '/SimpleTtest.mat']);
%     matlabbatch{1}.spm.stats.factorial_design.dir = {dirOutput};
%     matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = images{m};
%     matlabbatch{1}.spm.stats.factorial_design.masking.em = {explicitMask};
%     save(fullfile(dirOutput,'SimpleTtest.mat'), 'matlabbatch');
%     spm_jobman('initcfg')
%     spm_jobman('run', matlabbatch);
%     
%     % Estimate SPM
%     cd(dirOutput);
%     load([dirOutput '/SPM.mat']);
%     try
%         spm_spm(SPM);
%         
%         % Add contrasts
%         clear matlabbatch
%         matlabbatch{1}.spm.stats.con.spmmat = {[dirOutput '/SPM.mat']};
%         matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = '+ve';
%         matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
%         matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
%         matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = '-ve';
%         matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = -1;
%         matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
%         matlabbatch{1}.spm.stats.con.consess{3}.fcon.name = '+/-ve';
%         matlabbatch{1}.spm.stats.con.consess{3}.fcon.weights = 1;
%         matlabbatch{1}.spm.stats.con.consess{3}.fcon.sessrep = 'none';
%         matlabbatch{1}.spm.stats.con.delete = 0;
%         save(fullfile(dirOutput,'ContrastManager.mat'), 'matlabbatch');
%         spm_jobman('initcfg')
%         spm_jobman('run', matlabbatch);
%     end
%     
% end

%% T-test (for mutiple MEG timewindows and hemispheres)

%explicitMask = [];
explicitMask = '/imaging/es03/fMRI_2017/Masks/Bilateral_Temporal_Superior_AAL2.nii';

version = 'univariateCorrelationsMEGfMRI';
models = {'M' 'MM' 'All'};

for hem = 1:2
    
    for t=1:4
        
        clear images
        for k=1:length(SubjToAnalyze)
            
            SubjCurrent = Subj{SubjToAnalyze(k)};
            
            for m=1:length(models)
                images{m}{k,1} = [GLMAnalPD '/' SubjCurrent '/' version '/effect-map_' models{m} '_tWin' num2str(t) '_hem' num2str(hem) '.nii'];
            end
            
            V = spm_vol(images{1}{k,1});
            img1 = spm_read_vols(V);
            V = spm_vol(images{2}{k,1});
            img2 = spm_read_vols(V);
            
            img_new = (img2+img1)/2;
            V.fname = [GLMAnalPD '/' SubjCurrent '/' version '/effect-map_MM+M.nii'];
            spm_write_vol(V,img_new);
            images{4}{k,1} = V.fname;
            
        end
        
        models = {'M' 'MM' 'All' 'MM+M'};
        
        %spm fmri
        
        for m=1:length(models)
            
            dirOutput = [GLMAnalPD '/' version '/SPMs/Ttest_Masked_' models{m} '_tWin' num2str(t) '_hem' num2str(hem)];
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
        
    end
    
end