%% Set up global variables

%clear all

% make sure EEG modality of SPM software is selected
%spm('EEG');
%spm

% add required paths
addpath(pwd);
% addpath('/group/language/data/ediz.sohoglu/matlab/utilities/');
% addpath('/opt/neuromag/meg_pd_1.2/');

% define paths
rawpathstem = '/imaging/tc02/';
preprocessedpathstem = '/imaging/tc02/PINFA_preprocessed/';

% % define conditions
% conditions = {'Mismatch_4' 'Match_4' 'Mismatch_8' 'Match_8' 'Mismatch_16' 'Match_16'};
% 
% contrast_labels = {'Sum all conditions';'Match-MisMatch'; 'Clear minus Unclear'; 'Gradient difference M-MM'};
% contrast_weights = [1, 1, 1, 1, 1, 1; -1, 1, -1, 1, -1, 1; -1, -1, 0, 0, 1, 1; -1, 1, 0, 0, 1, -1];
conditions = {'Match 3' 'Match 16' 'Mismatch 3' 'Mismatch 16' 'Written'};

% define subjects and blocks (group(cnt) = 1 for controls initial visit, group(cnt) = 2 for patients initial visit, group(cnt) = 3 for controls follow up, group(cnt) = 4 for patients follow up)
cnt = 0;

% cnt = cnt + 1;
% subjects{cnt} = 'P7C01';
% dates{cnt} = '20180523';
% fullid{cnt} = '26934/20180523_U-ID40891';
% basedir{cnt} = 'PINFA';
% blocksin_folders{cnt} = {'Series_020_mp2rage_sag_p3_0.75mm_UNI_Images','Series_021_mp2rage_sag_p3_0.75mm_INV2','Series_010_cmrr_mbep2d_3x2_sparse_238vols','Series_012_cmrr_mbep2d_3x2_sparse_238vols','Series_016_cmrr_mbep2d_3x2_sparse_238vols','Series_018_cmrr_mbep2d_3x2_sparse_238vols','Series_011_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_013_cmrr_mbep2d_3x2_sparse_invPE_SBRef'};
% blocksin{cnt} = {'DATA_0020.nii', 'DATA_0021.nii', 'DATA_0010.nii', 'DATA_0012.nii', 'DATA_0016.nii', 'DATA_0018.nii', 'DATA_0011.nii','DATA_0013.nii'};
% blocksout{cnt} = {'structural', 'INV2', 'Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'};
% minvols(cnt) = 238;
% group(cnt) = 1;
% 
% cnt = cnt + 1;
% subjects{cnt} = 'P7P01';
% dates{cnt} = '20180524';
% fullid{cnt} = '26260/20180524_U-ID40644';
% basedir{cnt} = 'PINFA';
% blocksin_folders{cnt} = {'Series_023_mp2rage_sag_p3_0.75mm_UNI_Images','Series_022_mp2rage_sag_p3_0.75mm_INV2','Series_010_cmrr_mbep2d_3x2_sparse_238vols','Series_012_cmrr_mbep2d_3x2_sparse_238vols','Series_027_cmrr_mbep2d_3x2_sparse_238vols','Series_029_cmrr_mbep2d_3x2_sparse_238vols','Series_011_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_013_cmrr_mbep2d_3x2_sparse_invPE_SBRef'};
% blocksin{cnt} = {'DATA_0023.nii', 'DATA_0022.nii', 'DATA_0010.nii', 'DATA_0012.nii', 'DATA_0027.nii', 'DATA_0029.nii', 'DATA_0011.nii','DATA_0013.nii'};
% blocksout{cnt} = {'structural', 'INV2', 'Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'};
% minvols(cnt) = 238;
% group(cnt) = 2;
% 
% cnt = cnt + 1;
% subjects{cnt} = 'P7P02';
% dates{cnt} = '20180530';
% fullid{cnt} = '26933/20180530_U-ID40643';
% basedir{cnt} = 'PINFA';
% blocksin_folders{cnt} = {'Series_021_mp2rage_sag_p3_0.75mm_UNI_Images','Series_020_mp2rage_sag_p3_0.75mm_INV2','Series_010_cmrr_mbep2d_3x2_sparse_238vols','Series_012_cmrr_mbep2d_3x2_sparse_238vols','Series_016_cmrr_mbep2d_3x2_sparse_238vols','Series_018_cmrr_mbep2d_3x2_sparse_238vols','Series_011_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_013_cmrr_mbep2d_3x2_sparse_invPE_SBRef'};
% blocksin{cnt} = {'DATA_0021.nii', 'DATA_0020.nii', 'DATA_0010.nii', 'DATA_0012.nii', 'DATA_0016.nii', 'DATA_0018.nii', 'DATA_0011.nii','DATA_0013.nii'};
% blocksout{cnt} = {'structural', 'INV2', 'Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'};
% minvols(cnt) = 238;
% group(cnt) = 2;

cnt = cnt + 1;
subjects{cnt} = 'P7C02';
dates{cnt} = '20180611';
fullid{cnt} = '26999/20180611_U-ID41070';
basedir{cnt} = 'PINFA';
blocksin_folders{cnt} = {'Series_023_mp2rage_sag_p3_0.75mm_UNI_Images','Series_022_mp2rage_sag_p3_0.75mm_INV2','Series_010_cmrr_mbep2d_3x2_sparse_238vols','Series_012_cmrr_mbep2d_3x2_sparse_238vols','Series_018_cmrr_mbep2d_3x2_sparse_238vols','Series_020_cmrr_mbep2d_3x2_sparse_238vols','Series_011_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_013_cmrr_mbep2d_3x2_sparse_invPE_SBRef'};
blocksin{cnt} = {'DATA_0023.nii', 'DATA_0022.nii', 'DATA_0010.nii', 'DATA_0012.nii', 'DATA_0018.nii', 'DATA_0020.nii', 'DATA_0011.nii','DATA_0013.nii'};
blocksout{cnt} = {'structural', 'INV2', 'Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'};
minvols(cnt) = 238;
group(cnt) = 1;

cnt = cnt + 1;
subjects{cnt} = 'P7C03';
dates{cnt} = '20180612';
fullid{cnt} = '27000/20180612_U-ID41081';
basedir{cnt} = 'PINFA';
blocksin_folders{cnt} = {'Series_024_mp2rage_sag_p3_0.75mm_UNI_Images','Series_025_mp2rage_sag_p3_0.75mm_INV2','Series_014_cmrr_mbep2d_3x2_sparse_238vols','Series_016_cmrr_mbep2d_3x2_sparse_238vols','Series_020_cmrr_mbep2d_3x2_sparse_238vols','Series_022_cmrr_mbep2d_3x2_sparse_238vols','Series_015_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_017_cmrr_mbep2d_3x2_sparse_invPE_SBRef'};
blocksin{cnt} = {'DATA_0024.nii', 'DATA_0025.nii', 'DATA_0014.nii', 'DATA_0016.nii', 'DATA_0020.nii', 'DATA_0022.nii', 'DATA_0015.nii','DATA_0017.nii'};
blocksout{cnt} = {'structural', 'INV2', 'Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'};
minvols(cnt) = 238;
group(cnt) = 1;