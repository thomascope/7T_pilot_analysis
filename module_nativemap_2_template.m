function module_nativemap_2_template(GLMDir,downsamp_ratio,StrDir,runagain)
% Normalise effect-maps to MNI template

if ~exist('runagain','var')
    runagain = 0;
end

if ~exist('downsamp_ratio','var')
    downsamp_ratio = 1;
end
% 
% addpath('/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/RSA_scripts/es_scripts_fMRI')
% addpath('/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/RSA_scripts/decoding_toolbox_v3.999')
addpath(genpath('/group/language/data/ediz.sohoglu/matlab/rsatoolbox'));

versionCurrent = 'spearman';

    % Gather images for current subject
if downsamp_ratio == 1
    images = cellstr(spm_select('FPList', [GLMDir '/TDTcrossnobis/' versionCurrent '/'], '^effect-map_.*.nii'));
    images_done = cellstr(spm_select('FPList', [GLMDir '/TDTcrossnobis/' versionCurrent '/'], '^whireseffect-map_.*.nii'));
else
    images = cellstr(spm_select('FPList', [GLMDir '/TDTcrossnobis_downsamp_' num2str(downsamp_ratio) '/' versionCurrent '/'], '^effect-map_.*.nii'));
    images_done = cellstr(spm_select('FPList', [GLMDir '/TDTcrossnobis_downsamp_' num2str(downsamp_ratio) '/' versionCurrent '/'], '^whireseffect-map_.*.nii'));
end
if ~runagain
images = setdiff(images,strrep(images_done,'whireseffect-map','effect-map'));
end

% Write out masks
images_mask = {};
for i=1:length(images)
    V = spm_vol(images{i});
    Y = spm_read_vols(V);
    Y(~isnan(Y)) = 1;
    Y(isnan(Y)) = 0;
    images_mask{i,1} = strrep(images{i},'effect-map','nativeSpaceMask');
     saveMRImage(Y,images_mask{i,1},V.mat);
end

% Normalize
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = [images; images_mask];
matlabbatch{1}.spm.spatial.normalise.write.subj.def =  cellstr([StrDir 'mri/y_structural_csf.nii']);

if downsamp_ratio == 1
    save(fullfile(GLMDir,'TDTcrossnobis',versionCurrent,'NormalizeTDTcrossnobis.mat'), 'matlabbatch');
else
    save(fullfile(GLMDir,['TDTcrossnobis_downsamp_' num2str(downsamp_ratio)],versionCurrent,'NormalizeTDTcrossnobis.mat'), 'matlabbatch');
end
spm_jobman('initcfg')
spm_jobman('run', matlabbatch);

matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'whires';
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
spm_jobman('initcfg')
spm_jobman('run', matlabbatch);


clear images images_mask

% Smooth, mask and write out normalised images

% Gather images for current subject
if downsamp_ratio == 1
    images = cellstr(spm_select('FPList', [GLMDir '/TDTcrossnobis/' versionCurrent '/'], '^weffect-map_.*.nii'));
    images_mask = cellstr(spm_select('FPList', [GLMDir '/TDTcrossnobis/' versionCurrent '/'], '^wnativeSpaceMask_.*.nii'));
else
    images = cellstr(spm_select('FPList', [GLMDir '/TDTcrossnobis_downsamp_' num2str(downsamp_ratio) '/' versionCurrent '/'], '^weffect-map_.*.nii'));
    images_mask = cellstr(spm_select('FPList', [GLMDir '/TDTcrossnobis_downsamp_' num2str(downsamp_ratio) '/' versionCurrent '/'], '^wnativeSpaceMask_.*.nii'));
end

% Mask and smooth normalised data
mask_threshold = .05;
for i=1:length(images)
    % Fix normalised mask
    V = spm_vol(images_mask{i});
    fname_mask = V.fname;
    Y_mask = spm_read_vols(V);
    Y_mask(Y_mask<mask_threshold) = 0;
    Y_mask(isnan(Y_mask)) = 0;
    saveMRImage(Y_mask,fname_mask,V.mat);
    
    % Fix normalised r-map
    V = spm_vol(images{i});
    fname_image = V.fname;
    Y = spm_read_vols(V);
    Y(Y_mask<mask_threshold) = 0;
    Y_mask(isnan(Y_mask)) = 0;
    saveMRImage(Y,fname_image,V.mat);
    
%     % Smooth and (re)mask normalised r-map % I have taken this out to
%     % reduce storage use for now
%     fname_smoothed = strrep(images{i},'weffect-map_','sweffect-map_');
%     spm_smooth(images{i},fname_smoothed,[8 8 8]);
%     V = spm_vol(fname_smoothed);
%     Y = spm_read_vols(V);
%     Y(Y_mask<mask_threshold) = NaN;
%     saveMRImage(Y,fname_smoothed,V.mat);
end

    
