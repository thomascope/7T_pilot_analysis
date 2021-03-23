clearvars
close all
clc

ScriptsPD = '/group/language/data/ediz.sohoglu/projects/fMRI_2017/scripts_fMRI';
PreProcPD = '/imaging/es03/fMRI_2017/PreprocessAnalysis';
TempPD = '/group/language/data/ediz.sohoglu/projects/fMRI_2017/scripts_fMRI/Templates'; % Template PD
RoiPD = '/imaging/es03/fMRI_2017/ROIs';
SpmPD = '/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7219';
MarsbarPD = '/group/language/data/ediz.sohoglu/matlab/spm12/toolbox/marsbar';

Subj = {'subj1' 'subj2' 'subj3' 'subj4' 'subj5' 'subj6' 'subj7' 'subj8' 'subj9' 'subj10' 'subj11' 'subj12' 'subj13' 'subj14' 'subj15' 'subj16' 'subj17' 'subj18' 'subj19' 'subj20' 'subj21' 'subj22' 'subj23' 'subj24' 'subj25' 'subj26' 'subj27'};
SubjToAnalyze = [3,4,5,6,8,10,11,12,13,15,16,17,18,19,20,21,22,23,24,26,27];
%SubjToAnalyze = 5;

addpath(ScriptsPD); % So that any customised code in project scripts folder is prioritised (like spm_my_defaults)
addpath(SpmPD); % So that correct version of SPM is loaded
addpath(MarsbarPD);
spm('FMRI');

for k=1:length(SubjToAnalyze)
    
    SubjCurrent = Subj{SubjToAnalyze(k)};
    SubjDir = fullfile(PreProcPD,SubjCurrent);
    
    clear matlabbatch
    
    load([TempPD '/Normalize.mat']);
    
    % Gather ROIs
    RoiSubjDir = [SubjDir '/ROIs/'];
    if exist(RoiSubjDir,'dir'); rmdir(RoiSubjDir,'s'); end
    mkdir(RoiSubjDir);
    roi_filenames = dir([RoiPD '/*.nii']);
    for f=1:length(roi_filenames)
        copyfile([RoiPD '/' roi_filenames(f).name],[RoiSubjDir roi_filenames(f).name]);
    end
    Rois = cellstr(spm_select('FPList', RoiSubjDir, '.*\.nii'));
    
    % Get bounding box of native space images
    SessEpiDir = fullfile(SubjDir,'Functional','Sess_1');
    EpiFiles = cellstr(spm_select('FPList', SessEpiDir, '^mean.*\.nii$'));
    BB = spm_get_bbox(EpiFiles{1});
    
    % Deformation field
    StrDir = fullfile(SubjDir, 'Structural');
    DefFile = cellstr(spm_select('FPList', StrDir, '^iy.*\.nii$'));
    
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = Rois;
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = DefFile;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [3 3 3.75];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = BB;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 1;
    
    save(fullfile(SubjDir,'NormalizeRois.mat'), 'matlabbatch');
    spm_jobman('initcfg')
    spm_jobman('run', matlabbatch);
    
    for r=1:length(Rois)
        delete(Rois{r});
    end
    
    % Use Marsbar to match image dimensions of native-space ROI to those of functionals
    marsbar('on');
    for r=1:length(Rois)
        [sp,other] = mars_space(EpiFiles{1});
        [filepath,name,ext] = fileparts(Rois{r});
        mars_img2rois(fullfile(filepath,['w' name ext]),'./','tmp')
        roi = maroi('tmp_1_roi.mat');
        save_as_image(roi,fullfile(filepath,['w' name ext]),sp);
        delete('tmp_1_roi.mat');
        
        % For following subject and Right_HG_Syl2GradedNormCrossnobis ROI, normalization to native space results in
        % two non-contiguous clusters. Only the second cluster is needed
        if strcmp(SubjCurrent,'subj5') && ~isempty(strfind(Rois{r},'Right_HG_Syl2GradedNormCrossnobis'))
            roi = maroi('tmp_2_roi.mat');
            save_as_image(roi,fullfile(filepath,['w' name ext]),sp);
            delete('tmp_2_roi.mat');
        end
    end
    
end