function crop_batch(varargin)
% crop_batch.m
% Coregister source images to template and apply to other images.
% crop neck of source image
% This assumes you have subjects under a directory in_parent_dir
% with the structure /sub-21271/ses-20180711/
% Each session directory has files names in 'images' structure to coregister.
% The first image listed is the one coregistered to the target. 
% The remainig images are transformed in the same way. 
% That is it is assumed that images 2:N are already in register
% with the first image.
% Each session directory has no_images_to_crop images to crop 
% in the same way. The remaining images are reoriented but not cropped.
% The final image '_run-03_MT.nii' was added ad hoc. The code could be 
% changed to deal with varying numbers of images but it does not at present


spm('defaults','pet');
spm_jobman('initcfg');

%Change paths to match your setup
addpath('/home/spj24/spm/toolbox/rorden/crop');
in_parent_dir = '/scratch/wbic-beta/spj24/p00477/nifti/derivatives/ants/';

subs={
% 'sub-21271' 'ses-20180711'
% 'sub-21520' 'ses-20180718'
% 'sub-27124' 'ses-20180611'
'sub-27287' 'ses-20180723'
    };

% source_image='mp2rage.nii';

%Only first image is used as a coregistration source
%Amend manually for patients who have 2 MT runs and controls who have 3 
images={
    '_acq-mp2ragesag_out-uni_run-01_MP2RAGE.nii'
    '_acq-mp2ragesag_out-t1_run-01_MP2RAGE.nii'
    '_acq-mp2ragesag_out-inv1_run-01_MP2RAGE.nii'
    '_acq-mp2ragesag_out-inv2_run-01_MP2RAGE.nii'
    '_acq-20lcsnmtoff_run-01_MT.nii'
    '_run-01_MT.nii'
    '_run-02_MT.nii'
%     '_run-03_MT.nii'

    };

no_images_to_crop=4;

for j=1:numel(subs(:,1))
% pafor j=1:numel(subs(:,1))
    sub=subs{j,1};
    sess=subs{j,2};

    vols=fullfile(in_parent_dir,sub,sess,'anat',strcat(sub,'_',sess,images));
    %vols=vertcat(source,other);
    %Make sure to set to use elderly target (4) and crop first 4 images only
    nii_setOrigin12x(vols, 4, true, no_images_to_crop,[-76 -112 -68; 76 76 104]);
end
