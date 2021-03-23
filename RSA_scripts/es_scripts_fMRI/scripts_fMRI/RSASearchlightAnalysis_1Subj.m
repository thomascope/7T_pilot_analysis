function RSASearchlightAnalysis_1Subj(ProcDataDir,SubjID,RsaPD,TempPD,SpmPD,ScriptsPD)

% Recipe_fMRI_searchlight
%
% Cai Wingfield 11-2009, 2-2010, 3-2010, 8-2010
%__________________________________________________________________________
% Copyright (C) 2010 Medical Research Council

%%%%%%%%%%%%%%%%%%%%
%% Initialisation %%
%%%%%%%%%%%%%%%%%%%%

cd(ScriptsPD); % So that any customised code in project scripts folder is prioritised
toolboxRoot = RsaPD; addpath(genpath(toolboxRoot)); % Catch sight of the toolbox code

%% Run searchlight for each subject specified
    
userOptions = defineUserOptions(SubjID);
if ~exist(userOptions.rootPath,'dir')
    mkdir(userOptions.rootPath);
end
 
%%%%%%%%%%%%%%%%%%%%%%
%% Get names and locations of T images %%
%%%%%%%%%%%%%%%%%%%%%%

userOptions.sourceImage = 'average_spmt';
betas=betaCorrespondence(userOptions);

%%%%%%%%%%%%%%%%%%%%%%
%% Data preparation %%
%%%%%%%%%%%%%%%%%%%%%%

%fullBrainVols = fMRIDataPreparation('SPM', userOptions); % if using beta images
fullBrainVols = fMRIDataPreparation(betas, userOptions); % if using con or T images
binaryMasks_nS = fMRIMaskPreparation(userOptions);

%%%%%%%%%%%%%%%%%%%%%
%% RDM calculation %%
%%%%%%%%%%%%%%%%%%%%%

models = constructModelRDMs(modelRDMs_fromSimulations(), userOptions);
close all

% % Convert to 'weighted contrast' format?
% for m=1:length(models)
%    models(m).RDM(find(models(m).RDM==1)) = 1/length(find(models(m).RDM==1)); 
%    models(m).RDM(find(models(m).RDM==0)) = -1/length(find(models(m).RDM==0)); 
%    models(m).RDM(find(isnan(models(m).RDM))) = 0; 
% end

%%%%%%%%%%%%%%%%%
%% Searchlight %%
%%%%%%%%%%%%%%%%%

%fMRISearchlight(fullBrainVols, binaryMasks_nS, models, 'SPM', userOptions); % if using beta images
fMRISearchlight(fullBrainVols, binaryMasks_nS, models, betas,userOptions); % if using con or T images

%% Normalise RSA images to MNI template

addpath(SpmPD);
spm('FMRI');

load([TempPD '/Normalize.mat']);
    
RsaDir = userOptions.rootPath; % directory with RSA images

% Gather images for current subject
images = cellstr(spm_select('FPList', [RsaDir '/Maps/'], ['^' userOptions.analysisName '.*_rMap_.*' SubjID '.img']));

% Write out masks
images_mask = {};
for i=1:length(images)
    V = spm_vol(images{i});
    Y = spm_read_vols(V);
    Y(~isnan(Y)) = 1;
    Y(isnan(Y)) = 0;
    images_mask{i,1} = strrep(images{i},'rMap','nativeSpaceMask');
    saveMRImage(Y,images_mask{i,1},V.mat);
end

% Normalize
DefFile = cellstr(spm_select('FPList', [ProcDataDir '/' SubjID '/Structural/'], '^y.*.nii'));
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = [images; images_mask];
matlabbatch{1}.spm.spatial.normalise.write.subj.def = DefFile;

save(fullfile(RsaDir,'NormalizeRSA.mat'), 'matlabbatch');
spm_jobman('initcfg')
spm_jobman('run', matlabbatch);

clear images images_mask

% Smooth, mask and write out normalised images

% Gather images for current subject
images = cellstr(spm_select('FPList', [RsaDir '/Maps/'], ['^w' userOptions.analysisName '.*_rMap_.*' SubjID '.img']));
images_mask = cellstr(spm_select('FPList', [RsaDir '/Maps/'], ['^w' userOptions.analysisName '.*_nativeSpaceMask_.*' SubjID '.img']));

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
    
    % Fix normalised r-map
    V = spm_vol(images{i});
    fname_image = V.fname;
    Y = spm_read_vols(V);
    Y(abs(Y_mask)<mask_threshold) = 0;
    Y_mask(isnan(Y_mask)) = 0;
    saveMRImage(Y,fname_image,V.mat);
    
    % Smooth and (re)mask normalised r-map
    fname_smoothed = strrep(images{i},['w' userOptions.analysisName],['sw' userOptions.analysisName]);
    spm_smooth(images{i},fname_smoothed,[6 6 6]);
    V = spm_vol(fname_smoothed);
    Y = spm_read_vols(V);
    Y(Y_mask==0) = NaN;
    saveMRImage(Y,fname_smoothed,V.mat);
end