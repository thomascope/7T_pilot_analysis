clearvars

addpath /group/language/data/ediz.sohoglu/matlab/cvMANOVA_v3
addpath /imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7219
addpath(genpath('/group/language/data/ediz.sohoglu/matlab/rsatoolbox'));

Subj = {'subj1' 'subj2' 'subj3' 'subj4' 'subj5' 'subj6' 'subj7' 'subj8' 'subj9' 'subj10' 'subj11' 'subj12' 'subj13' 'subj14' 'subj15' 'subj16' 'subj17' 'subj18' 'subj19' 'subj20' 'subj21' 'subj22' 'subj23' 'subj24' 'subj25' 'subj26' 'subj27'};
%SubjToAnalyze = [3,4,5,6,8,10,11,12,13,15,16,17,18,19,20,21,22,23,24,26,27];
SubjToAnalyze = 3;
PreProcPD = '/imaging/es03/fMRI_2017/PreprocessAnalysis';
GLMAnalPD = '/imaging/es03/fMRI_2017/GLMAnalysisNative';
TempPD = '/group/language/data/ediz.sohoglu/projects/fMRI_2017/scripts_fMRI/Templates'; % Template PD

%% Parallel computing settings

addpath /hpc-software/matlab/cbu/

S = cbu_scheduler();
S.NumWorkers = length(SubjToAnalyze);
S.SubmitArguments = '-l mem=8GB -l walltime=1440';
%S.SubmitArguments = '-l mem=48GB -l walltime=10080';

%% For Syl 2 GLMs

% Global effects
Cs{1} = [zeros(15,64) zeros(15,32) orth(diff(eye(16))')']'; % Main effect of syllable 2 (Noise+Speech)
Cs{2} = [kron([1/2 1/2],kron([1/2 1/2],orth(diff(eye(16))')')) zeros(15,32) orth(diff(eye(16))')']'; % Main effect of syllable 2 (Clear+Clear,Noise+Speech)
Cs{3} = [kron([1 -1],kron([1/2 1/2],orth(diff(eye(16))')'))]'; % Syllable 2 X Prior congruency (Clear+Clear)
Cs{4} = [kron([1/2 1/2],kron([1 -1],orth(diff(eye(16))')'))]'; % Syllable 2 X Prior strength (Clear+Clear)
Cs{5} = [kron([1 -1],kron([1 -1],orth(diff(eye(16))')'))]'; % Syllable 2 X Prior congruency X Prior strength (Clear+Clear)
Cs{6} = [zeros(15,64) kron([1 -1],orth(diff(eye(16))')')]'; % Syllable 2 X Prior strength (Clear+Noise)

% Simple effects of syllable 2
Cs{7} = [orth(diff(eye(16))')']'; % Syllable 2 (StrongMatching+Clear)
Cs{8} = [zeros(15,16) orth(diff(eye(16))')']'; % Syllable 2 (WeakMatching+Clear)
Cs{9} = [zeros(15,16) zeros(15,16) orth(diff(eye(16))')']'; % Syllable 2 (StrongMismatching+Clear)
Cs{10} = [zeros(15,16) zeros(15,16) zeros(15,16) orth(diff(eye(16))')']'; % Syllable 2 (WeakMismatching+Clear)
Cs{11} = [zeros(15,16) zeros(15,16) zeros(15,16) zeros(15,16) orth(diff(eye(16))')']'; % Syllable 2 (Strong+Noise)
Cs{12} = [zeros(15,16) zeros(15,16) zeros(15,16) zeros(15,16)  zeros(15,16) orth(diff(eye(16))')']'; % Syllable 2 (Weak+Noise)
Cs{13} = [zeros(15,16) zeros(15,16) zeros(15,16) zeros(15,16) zeros(15,16) zeros(15,16) orth(diff(eye(16))')']'; % Syllable 2 (Noise+Speech)

%% For Syl1 GLMs

% Global effects
Cs{1} = [kron([1/2 1/2],orth(diff(eye(32))')')]'; % Main effect of syllable 1
Cs{2} = [kron([1 -1],orth(diff(eye(32))')')]'; % Syllable 1 X Prior strength

% Simple effects of syllable 1
Cs{3} = [kron([1 0],orth(diff(eye(32))')')]'; % Syllable 1 (Strong)
Cs{4} = [kron([0 1],orth(diff(eye(32))')')]'; % Syllable 1 (Weak)

%% For all pairwise comparisons

%[X,a] = indicatorMatrix('allpairs',1:64);
[X,a] = indicatorMatrix('allpairs',1:128);
for c=1:size(X,1); Cs{c} = X(c,:)'; end

%%

clear J
for k=1:length(SubjToAnalyze)
    
    SubjCurrent = Subj{SubjToAnalyze(k)};
    SubjDir = fullfile(GLMAnalPD,SubjCurrent);
    
    cvManovaSearchlight_es(SubjDir, 3, Cs, 0);
    J(k).task = @cvManovaSearchlight_es; % External function name here
    J(k).n_return_values = 0; % important
    J(k).input_args = {SubjDir,3,Cs,0};
    J(k).depends_on = 0;
    
end

%cbu_qsub(J, S);

%% Normalise MANOVA images to MNI template
spm fmri

clear matlabbatch

load([TempPD '/Normalize.mat']);
for k=1:length(SubjToAnalyze)
    
    SubjCurrent = Subj{SubjToAnalyze(k)};
    GLMAnalDir = [GLMAnalPD '/' SubjCurrent];
    
    % Gather images for current subject
    images = cellstr(spm_select('FPList', GLMAnalDir, '^spmDs_.*_P0001.nii$'));
    
    % Write out masks
    images_mask = {};
    for i=1:length(images)
        V = spm_vol(images{i});
        Y = spm_read_vols(V);
        Y(~isnan(Y)) = 1;
        Y(isnan(Y)) = 0;
        images_mask{i,1} = strrep(images{i},'.nii','_nativeSpaceMask.nii');
        spmWriteImage(Y,images_mask{i,1},V.mat);
    end
    
    % Deformation field
    DefFile = cellstr(spm_select('FPList', [PreProcPD '/' SubjCurrent '/Structural/'], '^y.*.nii'));
    
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = [images; images_mask];
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = DefFile;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 1;
    
    save(fullfile(GLMAnalDir,'NormalizeMANOVA.mat'), 'matlabbatch');
    spm_jobman('initcfg')
    spm_jobman('run', matlabbatch);
    
    % Smooth, mask and write out normalised images
    
    % Gather images for current subject
    images = cellstr(spm_select('FPList', GLMAnalDir, '^wspmDs_.*_P0001.nii$'));
    images_mask = cellstr(spm_select('FPList', GLMAnalDir, '^wspmDs.*_nativeSpaceMask.nii$'));
    
    % Mask and smooth normalised data
    mask_threshold = .01;
    for i=1:length(images)
        % Fix normalised mask
        V = spm_vol(images_mask{i});
        fname_mask = V.fname;
        Y_mask = spm_read_vols(V);
        Y_mask(abs(Y_mask)<mask_threshold) = 0;
        Y_mask(isnan(Y_mask)) = 0;
        saveMRImage(Y_mask,fname_mask,V.mat);
        
        % Fix normalised MANOVA-map
        V = spm_vol(images{i});
        fname_image = V.fname;
        Y = spm_read_vols(V);
        Y(abs(Y_mask)<mask_threshold) = 0;
        Y_mask(isnan(Y_mask)) = 0;
        saveMRImage(Y,fname_image,V.mat);
        
        % Smooth and (re)mask normalised MANOVA-map
        fname_smoothed = strrep(images{i},'wspmDs','swspmDs');
        spm_smooth(images{i},fname_smoothed,[6 6 6]);
        V = spm_vol(fname_smoothed);
        Y = spm_read_vols(V);
        Y(Y_mask==0) = NaN;
        saveMRImage(Y,fname_smoothed,V.mat);
    end
    
end

%% Make effect-maps

clear Models
%tmp = modelRDMs;
tmp = modelRDMs2;
Models(1).RDM = tmp.Syl2StrongM;
Models(1).name = 'Syl2StrongM';
Models(2).RDM = tmp.Syl2WeakM;
Models(2).name = 'Syl2WeakM';
Models(3).RDM = tmp.Syl2StrongMM;
Models(3).name = 'Syl2StrongMM';
Models(4).RDM = tmp.Syl2WeakMM;
Models(4).name = 'Syl2WeakMM';

for k=1:length(SubjToAnalyze)
    
    SubjCurrent = Subj{SubjToAnalyze(k)};
    GLMAnalDir = [GLMAnalPD '/' SubjCurrent];
    
    % Gather images for current subject
    images = cellstr(spm_select('FPList', GLMAnalDir, '^spmDs_.*_P0001.nii$'));
    
    % Write out masks
    clear Y_all
    for i=1:length(images)
        V = spm_vol(images{i});
        Y = spm_read_vols(V);
        Y_all(:,i) = Y(:);
    end
    
    for m=1:length(Models)
        %effect = corr(Y_all',vectorizeRDMs(Models(m).RDM(1:128,1:128))','type','Spearman','rows','pairwise');
        effect = Y_all*vectorizeRDMs(Models(m).RDM)';
        effect = reshape(effect,V.dim);
        saveMRImage(effect, [GLMAnalDir '/effect-map_' Models(m).name '.nii'], V.mat);
    end
    
end

%% Normalise effect-maps to MNI template

spm fmri

clear matlabbatch

modelNames = {'Syl2StrongM' 'Syl2WeakM' 'Syl2StrongMM' 'Syl2WeakMM'};

load([TempPD '/Normalize.mat']);
for k=1:length(SubjToAnalyze)
    
    SubjCurrent = Subj{SubjToAnalyze(k)};
    GLMAnalDir = [GLMAnalPD '/' SubjCurrent];
    
    % Gather images for current subject
    clear images
    for m=1:length(modelNames)
        images{m,1} = [GLMAnalDir '/effect-map_' modelNames{m} '.nii'];
    end
    
    % Write out masks
    images_mask = {};
    for i=1:length(images)
        V = spm_vol(images{i});
        Y = spm_read_vols(V);
        Y(~isnan(Y)) = 1;
        Y(isnan(Y)) = 0;
        images_mask{i,1} = strrep(images{i},'.nii','_nativeSpaceMask.nii');
        spmWriteImage(Y,images_mask{i,1},V.mat);
    end
    
    % Deformation field
    DefFile = cellstr(spm_select('FPList', [PreProcPD '/' SubjCurrent '/Structural/'], '^y.*.nii'));
    
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = [images; images_mask];
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = DefFile;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 1;
    
    save(fullfile(GLMAnalDir,'NormalizeMANOVA.mat'), 'matlabbatch');
    spm_jobman('initcfg')
    spm_jobman('run', matlabbatch);
    
    % Smooth, mask and write out normalised images
    
    % Gather images for current subject
    clear images images_mask
    for m=1:length(modelNames)
        images{m,1} = [GLMAnalDir '/weffect-map_' modelNames{m} '.nii'];
        images_mask{m,1} = [GLMAnalDir '/weffect-map_' modelNames{m} '_nativeSpaceMask.nii'];
    end
    
    % Mask and smooth normalised data
    mask_threshold = .01;
    for i=1:length(images)
        % Fix normalised mask
        V = spm_vol(images_mask{i});
        fname_mask = V.fname;
        Y_mask = spm_read_vols(V);
        Y_mask(abs(Y_mask)<mask_threshold) = 0;
        Y_mask(isnan(Y_mask)) = 0;
        saveMRImage(Y_mask,fname_mask,V.mat);
        
        % Fix normalised MANOVA-map
        V = spm_vol(images{i});
        fname_image = V.fname;
        Y = spm_read_vols(V);
        Y(abs(Y_mask)<mask_threshold) = 0;
        Y_mask(isnan(Y_mask)) = 0;
        saveMRImage(Y,fname_image,V.mat);
        
        % Smooth and (re)mask normalised MANOVA-map
        fname_smoothed = strrep(images{i},'weffect-map','sweffect-map');
        spm_smooth(images{i},fname_smoothed,[6 6 6]);
        V = spm_vol(fname_smoothed);
        Y = spm_read_vols(V);
        Y(Y_mask==0) = NaN;
        saveMRImage(Y,fname_smoothed,V.mat);
    end
    
end

%% Average images over conditions

for k=1:length(SubjToAnalyze)
    
    SubjCurrent = Subj{SubjToAnalyze(k)};    
    
    imagesCurrent = dir([GLMAnalPD '/' SubjCurrent '/' 'sweffect-map_*.nii']);
    
    clear data
    for c=1:length(imagesCurrent)
        
        image = [GLMAnalPD '/' SubjCurrent '/' imagesCurrent(c).name];
        V = spm_vol(image);
        data(:,:,:,c) = spm_read_vols(V);
        
    end
    
    data = mean(data,4);
    saveMRImage(data,[GLMAnalPD '/' SubjCurrent '/sweffect-map_mean.nii'],V.mat);
    
end


%% T-test

spm fmri

explicitMask = [];

clear images
for k=1:length(SubjToAnalyze)
    
    SubjCurrent = Subj{SubjToAnalyze(k)};    
    
    imagesCurrent = dir([GLMAnalPD '/' SubjCurrent '/' 'swspmDs_C*_P0001.nii']);
    
    for c=1:length(imagesCurrent)
        
        images{c}{k,1} = [GLMAnalPD '/' SubjCurrent '/' imagesCurrent(c).name];
        
    end
    
end

for c=1:length(images)
    
    dirOutput = sprintf([GLMAnalPD '/Manova/SPMs/Ttest_swspmDs_C%04d'],c);
    if exist(dirOutput,'dir');
        rmdir(dirOutput,'s');
    end
    mkdir(dirOutput);
       
    clear matlabbatch
    load([TempPD '/SimpleTtest.mat']);
    matlabbatch{1}.spm.stats.factorial_design.dir = {dirOutput};
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = images{c};
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