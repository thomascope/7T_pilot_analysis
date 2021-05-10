%% Set up global variables 

%clear all

% add required paths
addpath(pwd);
% addpath('/group/language/data/ediz.sohoglu/matlab/utilities/');
% addpath('/opt/neuromag/meg_pd_1.2/');

% define paths
rawpathstem = '/imaging/mlr/users/tc02/';
preprocessedpathstem = '/imaging/mlr/users/tc02/PINFA_preprocessed_2021/';

% % define conditions
% conditions = {'Mismatch_4' 'Match_4' 'Mismatch_8' 'Match_8' 'Mismatch_16' 'Match_16'};
% 
% contrast_labels = {'Sum all conditions';'Match-MisMatch'; 'Clear minus Unclear'; 'Gradient difference M-MM'};
% contrast_weights = [1, 1, 1, 1, 1, 1; -1, 1, -1, 1, -1, 1; -1, -1, 0, 0, 1, 1; -1, 1, 0, 0, 1, -1];
conditions = {'Match 3' 'Match 16' 'Mismatch 3' 'Mismatch 16' 'Written'};

%% Set up global variables 

cnt = 0;

cnt = cnt + 1; 
subjects{cnt} = 'P7C01'; 
dates{cnt} = '20180523'; 
fullid{cnt} = '26934/20180523_U-ID40891'; 
basedir{cnt} = 'PINFA'; 
blocksin_folders{cnt} = {'Series_020_mp2rage_sag_p3_0.75mm_UNI_Images','Series_021_mp2rage_sag_p3_0.75mm_INV2','Series_010_cmrr_mbep2d_3x2_sparse_238vols','Series_012_cmrr_mbep2d_3x2_sparse_238vols','Series_016_cmrr_mbep2d_3x2_sparse_238vols','Series_018_cmrr_mbep2d_3x2_sparse_238vols','Series_011_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_013_cmrr_mbep2d_3x2_sparse_invPE_SBRef'}; 
blocksin{cnt} = {'DATA_0020.nii','DATA_0021.nii','DATA_0010.nii','DATA_0012.nii','DATA_0016.nii','DATA_0018.nii','DATA_0011.nii','DATA_0013.nii'}; 
blocksout{cnt} = {'structural','INV2','Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'}; 
minvols(cnt) = 238; 
group(cnt) = 1; 

cnt = cnt + 1; 
subjects{cnt} = 'P7P01'; 
dates{cnt} = '20180524'; 
fullid{cnt} = '26260/20180524_U-ID40644'; 
basedir{cnt} = 'PINFA'; 
blocksin_folders{cnt} = {'Series_023_mp2rage_sag_p3_0.75mm_UNI_Images','Series_022_mp2rage_sag_p3_0.75mm_INV2','Series_010_cmrr_mbep2d_3x2_sparse_238vols','Series_012_cmrr_mbep2d_3x2_sparse_238vols','Series_027_cmrr_mbep2d_3x2_sparse_238vols','Series_029_cmrr_mbep2d_3x2_sparse_238vols','Series_011_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_013_cmrr_mbep2d_3x2_sparse_invPE_SBRef'}; 
blocksin{cnt} = {'DATA_0023.nii','DATA_0022.nii','DATA_0010.nii','DATA_0012.nii','DATA_0027.nii','DATA_0029.nii','DATA_0011.nii','DATA_0013.nii'}; 
blocksout{cnt} = {'structural','INV2','Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'}; 
minvols(cnt) = 238; 
group(cnt) = 2; 

cnt = cnt + 1; 
subjects{cnt} = 'P7P02'; 
dates{cnt} = '20180530'; 
fullid{cnt} = '26933/20180530_U-ID40643'; 
basedir{cnt} = 'PINFA'; 
blocksin_folders{cnt} = {'Series_021_mp2rage_sag_p3_0.75mm_UNI_Images','Series_020_mp2rage_sag_p3_0.75mm_INV2','Series_010_cmrr_mbep2d_3x2_sparse_238vols','Series_012_cmrr_mbep2d_3x2_sparse_238vols','Series_016_cmrr_mbep2d_3x2_sparse_238vols','Series_018_cmrr_mbep2d_3x2_sparse_238vols','Series_011_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_013_cmrr_mbep2d_3x2_sparse_invPE_SBRef'}; 
blocksin{cnt} = {'DATA_0021.nii','DATA_0020.nii','DATA_0010.nii','DATA_0012.nii','DATA_0016.nii','DATA_0018.nii','DATA_0011.nii','DATA_0013.nii'}; 
blocksout{cnt} = {'structural','INV2','Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'}; 
minvols(cnt) = 238; 
group(cnt) = 2; 

cnt = cnt + 1; 
subjects{cnt} = 'P7C02'; 
dates{cnt} = '20180611'; 
fullid{cnt} = '26999/20180611_U-ID41070'; 
basedir{cnt} = 'PINFA'; 
blocksin_folders{cnt} = {'Series_023_mp2rage_sag_p3_0.75mm_UNI_Images','Series_022_mp2rage_sag_p3_0.75mm_INV2','Series_010_cmrr_mbep2d_3x2_sparse_238vols','Series_012_cmrr_mbep2d_3x2_sparse_238vols','Series_018_cmrr_mbep2d_3x2_sparse_238vols','Series_020_cmrr_mbep2d_3x2_sparse_238vols','Series_011_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_013_cmrr_mbep2d_3x2_sparse_invPE_SBRef'}; 
blocksin{cnt} = {'DATA_0023.nii','DATA_0022.nii','DATA_0010.nii','DATA_0012.nii','DATA_0018.nii','DATA_0020.nii','DATA_0011.nii','DATA_0013.nii'}; 
blocksout{cnt} = {'structural','INV2','Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'}; 
minvols(cnt) = 238; 
group(cnt) = 1; 

cnt = cnt + 1; 
subjects{cnt} = 'P7C03'; 
dates{cnt} = '20180612'; 
fullid{cnt} = '27000/20180612_U-ID41081'; 
basedir{cnt} = 'PINFA'; 
blocksin_folders{cnt} = {'Series_024_mp2rage_sag_p3_0.75mm_UNI_Images','Series_025_mp2rage_sag_p3_0.75mm_INV2','Series_014_cmrr_mbep2d_3x2_sparse_238vols','Series_016_cmrr_mbep2d_3x2_sparse_238vols','Series_020_cmrr_mbep2d_3x2_sparse_238vols','Series_022_cmrr_mbep2d_3x2_sparse_238vols','Series_015_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_017_cmrr_mbep2d_3x2_sparse_invPE_SBRef'}; 
blocksin{cnt} = {'DATA_0024.nii','DATA_0025.nii','DATA_0014.nii','DATA_0016.nii','DATA_0020.nii','DATA_0022.nii','DATA_0015.nii','DATA_0017.nii'}; 
blocksout{cnt} = {'structural','INV2','Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'}; 
minvols(cnt) = 238; 
group(cnt) = 1; 

cnt = cnt + 1; 
subjects{cnt} = 'P7P03'; 
dates{cnt} = '20180731'; 
fullid{cnt} = '26043/20180731_U-ID41556'; 
basedir{cnt} = 'PINFA'; 
blocksin_folders{cnt} = {'Series_024_mp2rage_sag_p3_0.75mm_UNI_Images','Series_025_mp2rage_sag_p3_0.75mm_INV2','Series_010_cmrr_mbep2d_3x2_sparse_238vols','Series_014_cmrr_mbep2d_3x2_sparse_238vols','Series_020_cmrr_mbep2d_3x2_sparse_238vols','Series_022_cmrr_mbep2d_3x2_sparse_238vols','Series_013_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_015_cmrr_mbep2d_3x2_sparse_invPE_SBRef'}; 
blocksin{cnt} = {'DATA_0024.nii','DATA_0025.nii','DATA_0010.nii','DATA_0014.nii','DATA_0020.nii','DATA_0022.nii','DATA_0013.nii','DATA_0015.nii'}; 
blocksout{cnt} = {'structural','INV2','Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'}; 
minvols(cnt) = 238; 
group(cnt) = 2; 

cnt = cnt + 1; 
subjects{cnt} = 'P7P04'; 
dates{cnt} = '20180810'; 
fullid{cnt} = '26773/20180810_U-ID41293'; 
basedir{cnt} = 'PINFA'; 
blocksin_folders{cnt} = {'Series_021_mp2rage_sag_p3_0.75mm_UNI_Images','Series_020_mp2rage_sag_p3_0.75mm_INV2','Series_010_cmrr_mbep2d_3x2_sparse_238vols','Series_012_cmrr_mbep2d_3x2_sparse_238vols','Series_016_cmrr_mbep2d_3x2_sparse_238vols','Series_018_cmrr_mbep2d_3x2_sparse_238vols','Series_011_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_013_cmrr_mbep2d_3x2_sparse_invPE_SBRef'}; 
blocksin{cnt} = {'DATA_0021.nii','DATA_0020.nii','DATA_0010.nii','DATA_0012.nii','DATA_0016.nii','DATA_0018.nii','DATA_0011.nii','DATA_0013.nii'}; 
blocksout{cnt} = {'structural','INV2','Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'}; 
minvols(cnt) = 238; 
group(cnt) = 2; 

cnt = cnt + 1; 
subjects{cnt} = 'P7P05'; 
dates{cnt} = '20180813'; 
fullid{cnt} = '27379/20180813_U-ID41675'; 
basedir{cnt} = 'PINFA'; 
blocksin_folders{cnt} = {'Series_027_mp2rage_sag_p3_0.75mm_UNI_Images','Series_026_mp2rage_sag_p3_0.75mm_INV2','Series_010_cmrr_mbep2d_3x2_sparse_238vols','Series_014_cmrr_mbep2d_3x2_sparse_238vols','Series_020_cmrr_mbep2d_3x2_sparse_238vols','Series_024_cmrr_mbep2d_3x2_sparse_238vols','Series_013_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_015_cmrr_mbep2d_3x2_sparse_invPE_SBRef'}; 
blocksin{cnt} = {'DATA_0027.nii','DATA_0026.nii','DATA_0010.nii','DATA_0014.nii','DATA_0020.nii','DATA_0024.nii','DATA_0013.nii','DATA_0015.nii'}; 
blocksout{cnt} = {'structural','INV2','Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'}; 
minvols(cnt) = 238; 
group(cnt) = 2; 

cnt = cnt + 1; 
subjects{cnt} = 'P7P06'; 
dates{cnt} = '20180816'; 
fullid{cnt} = '27393/20180816_U-ID41711'; 
basedir{cnt} = 'PINFA'; 
blocksin_folders{cnt} = {'Series_028_mp2rage_sag_p3_0.75mm_UNI_Images','Series_029_mp2rage_sag_p3_0.75mm_INV2','Series_014_cmrr_mbep2d_3x2_sparse_238vols','Series_016_cmrr_mbep2d_3x2_sparse_238vols','Series_020_cmrr_mbep2d_3x2_sparse_238vols','Series_022_cmrr_mbep2d_3x2_sparse_238vols','Series_015_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_017_cmrr_mbep2d_3x2_sparse_invPE_SBRef'}; 
blocksin{cnt} = {'DATA_0028.nii','DATA_0029.nii','DATA_0014.nii','DATA_0016.nii','DATA_0020.nii','DATA_0022.nii','DATA_0015.nii','DATA_0017.nii'}; 
blocksout{cnt} = {'structural','INV2','Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'}; 
minvols(cnt) = 238; 
group(cnt) = 2; 

cnt = cnt + 1; 
subjects{cnt} = 'P7P07'; 
dates{cnt} = '20180911'; 
fullid{cnt} = '23483/20180911_U-ID41986'; 
basedir{cnt} = 'PINFA'; 
blocksin_folders{cnt} = {'Series_020_mp2rage_sag_p3_0.75mm_UNI_Images','Series_021_mp2rage_sag_p3_0.75mm_INV2','Series_010_cmrr_mbep2d_3x2_sparse_238vols','Series_012_cmrr_mbep2d_3x2_sparse_238vols','Series_016_cmrr_mbep2d_3x2_sparse_238vols','Series_018_cmrr_mbep2d_3x2_sparse_238vols','Series_011_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_013_cmrr_mbep2d_3x2_sparse_invPE_SBRef'}; 
blocksin{cnt} = {'DATA_0020.nii','DATA_0021.nii','DATA_0010.nii','DATA_0012.nii','DATA_0016.nii','DATA_0018.nii','DATA_0011.nii','DATA_0013.nii'}; 
blocksout{cnt} = {'structural','INV2','Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'}; 
minvols(cnt) = 238; 
group(cnt) = 2; 

cnt = cnt + 1; 
subjects{cnt} = 'P7C05'; 
dates{cnt} = '20180927'; 
fullid{cnt} = '26566/20180927_U-ID42189'; 
basedir{cnt} = 'PINFA'; 
blocksin_folders{cnt} = {'Series_031_mp2rage_sag_p3_0.75mm_UNI_Images','Series_030_mp2rage_sag_p3_0.75mm_INV2','Series_020_cmrr_mbep2d_3x2_sparse_238vols','Series_022_cmrr_mbep2d_3x2_sparse_238vols','Series_026_cmrr_mbep2d_3x2_sparse_238vols','Series_028_cmrr_mbep2d_3x2_sparse_238vols','Series_021_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_023_cmrr_mbep2d_3x2_sparse_invPE_SBRef'}; 
blocksin{cnt} = {'DATA_0031.nii','DATA_0030.nii','DATA_0020.nii','DATA_0022.nii','DATA_0026.nii','DATA_0028.nii','DATA_0021.nii','DATA_0023.nii'}; 
blocksout{cnt} = {'structural','INV2','Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'}; 
minvols(cnt) = 238; 
group(cnt) = 1; 

cnt = cnt + 1;
subjects{cnt} = 'P7P08';
dates{cnt} = '20180920';
fullid{cnt} = '27567/20180920_U-ID42098';
basedir{cnt} = 'PINFA';
blocksin_folders{cnt} = {'Series_029_mp2rage_sag_p3_0.75mm_UNI_Images','Series_028_mp2rage_sag_p3_0.75mm_INV2','Series_012_cmrr_mbep2d_3x2_sparse_238vols','Series_020_cmrr_mbep2d_3x2_sparse_238vols','Series_024_cmrr_mbep2d_3x2_sparse_238vols','Series_026_cmrr_mbep2d_3x2_sparse_238vols','Series_019_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_021_cmrr_mbep2d_3x2_sparse_invPE_SBRef'};
blocksin{cnt} = {'DATA_0029.nii', 'DATA_0028.nii', 'DATA_0012.nii', 'DATA_0020.nii', 'DATA_0024.nii', 'DATA_0026.nii', 'DATA_0019.nii','DATA_0021.nii'};
blocksout{cnt} = {'structural', 'INV2', 'Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'};
minvols(cnt) = 238;
group(cnt) = 2;

cnt = cnt + 1; 
subjects{cnt} = 'P7P09'; 
dates{cnt} = '20181009'; 
fullid{cnt} = '27614/20181009_U-ID42215'; 
basedir{cnt} = 'PINFA'; 
blocksin_folders{cnt} = {'Series_021_mp2rage_sag_p3_0.75mm_UNI_Images','Series_020_mp2rage_sag_p3_0.75mm_INV2','Series_010_cmrr_mbep2d_3x2_sparse_238vols','Series_012_cmrr_mbep2d_3x2_sparse_238vols','Series_016_cmrr_mbep2d_3x2_sparse_238vols','Series_018_cmrr_mbep2d_3x2_sparse_238vols','Series_011_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_013_cmrr_mbep2d_3x2_sparse_invPE_SBRef'}; 
blocksin{cnt} = {'DATA_0021.nii','DATA_0020.nii','DATA_0010.nii','DATA_0012.nii','DATA_0016.nii','DATA_0018.nii','DATA_0011.nii','DATA_0013.nii'}; 
blocksout{cnt} = {'structural','INV2','Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'}; 
minvols(cnt) = 238; 
group(cnt) = 2; 

cnt = cnt + 1; 
subjects{cnt} = 'P7C06'; 
dates{cnt} = '20181010'; 
fullid{cnt} = '27392/20181010_U-ID41709'; 
basedir{cnt} = 'PINFA'; 
blocksin_folders{cnt} = {'Series_020_mp2rage_sag_p3_0.75mm_UNI_Images','Series_021_mp2rage_sag_p3_0.75mm_INV2','Series_010_cmrr_mbep2d_3x2_sparse_238vols','Series_012_cmrr_mbep2d_3x2_sparse_238vols','Series_016_cmrr_mbep2d_3x2_sparse_238vols','Series_018_cmrr_mbep2d_3x2_sparse_238vols','Series_011_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_013_cmrr_mbep2d_3x2_sparse_invPE_SBRef'}; 
blocksin{cnt} = {'DATA_0020.nii','DATA_0021.nii','DATA_0010.nii','DATA_0012.nii','DATA_0016.nii','DATA_0018.nii','DATA_0011.nii','DATA_0013.nii'}; 
blocksout{cnt} = {'structural','INV2','Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'}; 
minvols(cnt) = 238; 
group(cnt) = 1; 

cnt = cnt + 1; 
subjects{cnt} = 'P7P10'; 
dates{cnt} = '20181022'; 
fullid{cnt} = '27673/20181022_U-ID42394'; 
basedir{cnt} = 'PINFA'; 
blocksin_folders{cnt} = {'Series_021_mp2rage_sag_p3_0.75mm_UNI_Images','Series_020_mp2rage_sag_p3_0.75mm_INV2','Series_010_cmrr_mbep2d_3x2_sparse_238vols','Series_012_cmrr_mbep2d_3x2_sparse_238vols','Series_016_cmrr_mbep2d_3x2_sparse_238vols','Series_018_cmrr_mbep2d_3x2_sparse_238vols','Series_011_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_013_cmrr_mbep2d_3x2_sparse_invPE_SBRef'}; 
blocksin{cnt} = {'DATA_0021.nii','DATA_0020.nii','DATA_0010.nii','DATA_0012.nii','DATA_0016.nii','DATA_0018.nii','DATA_0011.nii','DATA_0013.nii'}; 
blocksout{cnt} = {'structural','INV2','Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'}; 
minvols(cnt) = 238; 
group(cnt) = 2; 

cnt = cnt + 1; 
subjects{cnt} = 'P7C08'; 
dates{cnt} = '20181024'; 
fullid{cnt} = '27685/20181024_U-ID42444'; 
basedir{cnt} = 'PINFA'; 
blocksin_folders{cnt} = {'Series_021_mp2rage_sag_p3_0.75mm_UNI_Images','Series_020_mp2rage_sag_p3_0.75mm_INV2','Series_010_cmrr_mbep2d_3x2_sparse_238vols','Series_012_cmrr_mbep2d_3x2_sparse_238vols','Series_016_cmrr_mbep2d_3x2_sparse_238vols','Series_018_cmrr_mbep2d_3x2_sparse_238vols','Series_011_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_013_cmrr_mbep2d_3x2_sparse_invPE_SBRef'}; 
blocksin{cnt} = {'DATA_0021.nii','DATA_0020.nii','DATA_0010.nii','DATA_0012.nii','DATA_0016.nii','DATA_0018.nii','DATA_0011.nii','DATA_0013.nii'}; 
blocksout{cnt} = {'structural','INV2','Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'}; 
minvols(cnt) = 238; 
group(cnt) = 1; 

cnt = cnt + 1; 
subjects{cnt} = 'P7C09'; 
dates{cnt} = '20181031'; 
fullid{cnt} = '27723/20181031_U-ID42526'; 
basedir{cnt} = 'PINFA'; 
blocksin_folders{cnt} = {'Series_023_mp2rage_sag_p3_0.75mm_UNI_Images','Series_022_mp2rage_sag_p3_0.75mm_INV2','Series_010_cmrr_mbep2d_3x2_sparse_238vols','Series_012_cmrr_mbep2d_3x2_sparse_238vols','Series_018_cmrr_mbep2d_3x2_sparse_238vols','Series_020_cmrr_mbep2d_3x2_sparse_238vols','Series_011_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_013_cmrr_mbep2d_3x2_sparse_invPE_SBRef'}; 
blocksin{cnt} = {'DATA_0023.nii','DATA_0022.nii','DATA_0010.nii','DATA_0012.nii','DATA_0018.nii','DATA_0020.nii','DATA_0011.nii','DATA_0013.nii'}; 
blocksout{cnt} = {'structural','INV2','Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'}; 
minvols(cnt) = 238; 
group(cnt) = 1; 

cnt = cnt + 1; 
subjects{cnt} = 'P7C10'; 
dates{cnt} = '20181102'; 
fullid{cnt} = '27391/20181102_U-ID42554'; 
basedir{cnt} = 'PINFA'; 
blocksin_folders{cnt} = {'Series_020_mp2rage_sag_p3_0.75mm_UNI_Images','Series_021_mp2rage_sag_p3_0.75mm_INV2','Series_010_cmrr_mbep2d_3x2_sparse_238vols','Series_012_cmrr_mbep2d_3x2_sparse_238vols','Series_016_cmrr_mbep2d_3x2_sparse_238vols','Series_018_cmrr_mbep2d_3x2_sparse_238vols','Series_011_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_013_cmrr_mbep2d_3x2_sparse_invPE_SBRef'}; 
blocksin{cnt} = {'DATA_0020.nii','DATA_0021.nii','DATA_0010.nii','DATA_0012.nii','DATA_0016.nii','DATA_0018.nii','DATA_0011.nii','DATA_0013.nii'}; 
blocksout{cnt} = {'structural','INV2','Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'}; 
minvols(cnt) = 238; 
group(cnt) = 1; 

cnt = cnt + 1; 
subjects{cnt} = 'P7P11'; 
dates{cnt} = '20181123'; 
fullid{cnt} = '27722/20181123_U-ID42842'; 
basedir{cnt} = 'PINFA'; 
blocksin_folders{cnt} = {'Series_020_mp2rage_sag_p3_0.75mm_UNI_Images','Series_021_mp2rage_sag_p3_0.75mm_INV2','Series_010_cmrr_mbep2d_3x2_sparse_238vols','Series_012_cmrr_mbep2d_3x2_sparse_238vols','Series_016_cmrr_mbep2d_3x2_sparse_238vols','Series_018_cmrr_mbep2d_3x2_sparse_238vols','Series_011_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_013_cmrr_mbep2d_3x2_sparse_invPE_SBRef'}; 
blocksin{cnt} = {'DATA_0020.nii','DATA_0021.nii','DATA_0010.nii','DATA_0012.nii','DATA_0016.nii','DATA_0018.nii','DATA_0011.nii','DATA_0013.nii'}; 
blocksout{cnt} = {'structural','INV2','Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'}; 
minvols(cnt) = 238; 
group(cnt) = 2; 

cnt = cnt + 1; 
subjects{cnt} = 'P7C11'; 
dates{cnt} = '20181127'; 
fullid{cnt} = '27909/20181127_U-ID42870'; 
basedir{cnt} = 'PINFA'; 
blocksin_folders{cnt} = {'Series_021_mp2rage_sag_p3_0.75mm_UNI_Images','Series_020_mp2rage_sag_p3_0.75mm_INV2','Series_010_cmrr_mbep2d_3x2_sparse_238vols','Series_012_cmrr_mbep2d_3x2_sparse_238vols','Series_016_cmrr_mbep2d_3x2_sparse_238vols','Series_018_cmrr_mbep2d_3x2_sparse_238vols','Series_011_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_013_cmrr_mbep2d_3x2_sparse_invPE_SBRef'}; 
blocksin{cnt} = {'DATA_0021.nii','DATA_0020.nii','DATA_0010.nii','DATA_0012.nii','DATA_0016.nii','DATA_0018.nii','DATA_0011.nii','DATA_0013.nii'}; 
blocksout{cnt} = {'structural','INV2','Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'}; 
minvols(cnt) = 238; 
group(cnt) = 1; 

cnt = cnt + 1; 
subjects{cnt} = 'P7C12'; 
dates{cnt} = '20200924'; 
fullid{cnt} = '31118/20200924_U-ID49997'; 
basedir{cnt} = 'PINFA'; 
blocksin_folders{cnt} = {'Series_023_mp2rage_sag_p3_0.75mm_UNI_Images','Series_022_mp2rage_sag_p3_0.75mm_INV2','Series_010_cmrr_mbep2d_3x2_sparse_238vols','Series_014_cmrr_mbep2d_3x2_sparse_238vols','Series_018_cmrr_mbep2d_3x2_sparse_238vols','Series_020_cmrr_mbep2d_3x2_sparse_238vols','Series_013_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_015_cmrr_mbep2d_3x2_sparse_invPE_SBRef'}; 
blocksin{cnt} = {'DATA_0023.nii','DATA_0022.nii','DATA_0010.nii','DATA_0014.nii','DATA_0018.nii','DATA_0020.nii','DATA_0013.nii','DATA_0015.nii'}; 
blocksout{cnt} = {'structural','INV2','Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'}; 
minvols(cnt) = 238; 
group(cnt) = 1; 

cnt = cnt + 1; 
subjects{cnt} = 'P7P12'; 
dates{cnt} = '20201001'; 
fullid{cnt} = '31143/20201001_U-ID'; 
basedir{cnt} = 'PINFA'; 
blocksin_folders{cnt} = {'Series_023_mp2rage_sag_p3_0.75mm_UNI_Images','Series_022_mp2rage_sag_p3_0.75mm_INV2','Series_010_cmrr_mbep2d_3x2_sparse_238vols','Series_014_cmrr_mbep2d_3x2_sparse_238vols','Series_018_cmrr_mbep2d_3x2_sparse_238vols','Series_020_cmrr_mbep2d_3x2_sparse_238vols','Series_013_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_015_cmrr_mbep2d_3x2_sparse_invPE_SBRef'}; 
blocksin{cnt} = {'DATA_0023.nii','DATA_0022.nii','DATA_0010.nii','DATA_0014.nii','DATA_0018.nii','DATA_0020.nii','DATA_0013.nii','DATA_0015.nii'}; 
blocksout{cnt} = {'structural','INV2','Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'}; 
minvols(cnt) = 238; 
group(cnt) = 2; 

cnt = cnt + 1; 
subjects{cnt} = 'P7P13'; 
dates{cnt} = '20201005'; 
fullid{cnt} = '29902/20201005_U-ID50131'; 
basedir{cnt} = 'PINFA'; 
blocksin_folders{cnt} = {'Series_020_mp2rage_sag_p3_0.75mm_UNI_Images','Series_021_mp2rage_sag_p3_0.75mm_INV2','Series_010_cmrr_mbep2d_3x2_sparse_238vols','Series_012_cmrr_mbep2d_3x2_sparse_238vols','Series_016_cmrr_mbep2d_3x2_sparse_238vols','Series_018_cmrr_mbep2d_3x2_sparse_238vols','Series_011_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_013_cmrr_mbep2d_3x2_sparse_invPE_SBRef'}; 
blocksin{cnt} = {'DATA_0020.nii','DATA_0021.nii','DATA_0010.nii','DATA_0012.nii','DATA_0016.nii','DATA_0018.nii','DATA_0011.nii','DATA_0013.nii'}; 
blocksout{cnt} = {'structural','INV2','Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'}; 
minvols(cnt) = 238; 
group(cnt) = 2; 

cnt = cnt + 1; 
subjects{cnt} = 'P7C13'; 
dates{cnt} = '20201012'; 
fullid{cnt} = '27962/20201012_U-ID50248'; 
basedir{cnt} = 'PINFA'; 
blocksin_folders{cnt} = {'Series_021_mp2rage_sag_p3_0.75mm_UNI_Images','Series_020_mp2rage_sag_p3_0.75mm_INV2','Series_010_cmrr_mbep2d_3x2_sparse_238vols','Series_012_cmrr_mbep2d_3x2_sparse_238vols','Series_016_cmrr_mbep2d_3x2_sparse_238vols','Series_018_cmrr_mbep2d_3x2_sparse_238vols','Series_011_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_013_cmrr_mbep2d_3x2_sparse_invPE_SBRef'}; 
blocksin{cnt} = {'DATA_0021.nii','DATA_0020.nii','DATA_0010.nii','DATA_0012.nii','DATA_0016.nii','DATA_0018.nii','DATA_0011.nii','DATA_0013.nii'}; 
blocksout{cnt} = {'structural','INV2','Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'}; 
minvols(cnt) = 238; 
group(cnt) = 1; 

% Subject excluded - unable to do the task, Pretty minimal, but present,
% normal > written univariate response, and only right parietal clarity/congruency interaction
% cnt = cnt + 1;
% subjects{cnt} = 'P7P14';
% dates{cnt} = '20201008';
% fullid{cnt} = '27567/20201008_U-ID50185';
% basedir{cnt} = 'PINFA';
% blocksin_folders{cnt} = {'Series_027_mp2rage_sag_p3_0.75mm_UNI_Images','Series_026_mp2rage_sag_p3_0.75mm_INV2','Series_016_cmrr_mbep2d_3x2_sparse_238vols','Series_018_cmrr_mbep2d_3x2_sparse_238vols','Series_022_cmrr_mbep2d_3x2_sparse_238vols','Series_024_cmrr_mbep2d_3x2_sparse_238vols','Series_017_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_019_cmrr_mbep2d_3x2_sparse_invPE_SBRef'};
% blocksin{cnt} = {'DATA_0027.nii', 'DATA_0026.nii', 'DATA_0016.nii', 'DATA_0018.nii', 'DATA_0022.nii', 'DATA_0024.nii', 'DATA_0017.nii','DATA_0019.nii'};
% blocksout{cnt} = {'structural', 'INV2', 'Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'};
% minvols(cnt) = 238;
% group(cnt) = 2;

cnt = cnt + 1; 
subjects{cnt} = 'P7P15'; 
dates{cnt} = '20201021'; 
fullid{cnt} = '29217/20201021_U-ID50279'; 
basedir{cnt} = 'PINFA'; 
blocksin_folders{cnt} = {'Series_023_mp2rage_sag_p3_0.75mm_UNI_Images','Series_022_mp2rage_sag_p3_0.75mm_INV2','Series_012_cmrr_mbep2d_3x2_sparse_238vols','Series_014_cmrr_mbep2d_3x2_sparse_238vols','Series_018_cmrr_mbep2d_3x2_sparse_238vols','Series_020_cmrr_mbep2d_3x2_sparse_238vols','Series_013_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_015_cmrr_mbep2d_3x2_sparse_invPE_SBRef'}; 
blocksin{cnt} = {'DATA_0023.nii','DATA_0022.nii','DATA_0012.nii','DATA_0014.nii','DATA_0018.nii','DATA_0020.nii','DATA_0013.nii','DATA_0015.nii'}; 
blocksout{cnt} = {'structural','INV2','Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'}; 
minvols(cnt) = 238; 
group(cnt) = 2; 

cnt = cnt + 1; 
subjects{cnt} = 'P7P16'; 
dates{cnt} = '20201023'; 
fullid{cnt} = '31271/20201023_U-ID50387'; 
basedir{cnt} = 'PINFA'; 
blocksin_folders{cnt} = {'Series_025_mp2rage_sag_p3_0.75mm_UNI_Images','Series_024_mp2rage_sag_p3_0.75mm_INV2','Series_014_cmrr_mbep2d_3x2_sparse_238vols','Series_016_cmrr_mbep2d_3x2_sparse_238vols','Series_020_cmrr_mbep2d_3x2_sparse_238vols','Series_022_cmrr_mbep2d_3x2_sparse_238vols','Series_015_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_017_cmrr_mbep2d_3x2_sparse_invPE_SBRef'}; 
blocksin{cnt} = {'DATA_0025.nii','DATA_0024.nii','DATA_0014.nii','DATA_0016.nii','DATA_0020.nii','DATA_0022.nii','DATA_0015.nii','DATA_0017.nii'}; 
blocksout{cnt} = {'structural','INV2','Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'}; 
minvols(cnt) = 238; 
group(cnt) = 2; 

cnt = cnt + 1; 
subjects{cnt} = 'P7C15'; 
dates{cnt} = '20201022'; 
fullid{cnt} = '26654/20201022_U-ID50276'; 
basedir{cnt} = 'PINFA'; 
blocksin_folders{cnt} = {'Series_020_mp2rage_sag_p3_0.75mm_UNI_Images','Series_021_mp2rage_sag_p3_0.75mm_INV2','Series_010_cmrr_mbep2d_3x2_sparse_238vols','Series_012_cmrr_mbep2d_3x2_sparse_238vols','Series_016_cmrr_mbep2d_3x2_sparse_238vols','Series_018_cmrr_mbep2d_3x2_sparse_238vols','Series_011_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_013_cmrr_mbep2d_3x2_sparse_invPE_SBRef'}; 
blocksin{cnt} = {'DATA_0020.nii','DATA_0021.nii','DATA_0010.nii','DATA_0012.nii','DATA_0016.nii','DATA_0018.nii','DATA_0011.nii','DATA_0013.nii'}; 
blocksout{cnt} = {'structural','INV2','Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'}; 
minvols(cnt) = 238; 
group(cnt) = 1; 

cnt = cnt + 1; 
subjects{cnt} = 'P7C16'; 
dates{cnt} = '20201029'; 
fullid{cnt} = '31233/20201029_U-ID50280'; 
basedir{cnt} = 'PINFA'; 
blocksin_folders{cnt} = {'Series_020_mp2rage_sag_p3_0.75mm_UNI_Images','Series_021_mp2rage_sag_p3_0.75mm_INV2','Series_010_cmrr_mbep2d_3x2_sparse_238vols','Series_012_cmrr_mbep2d_3x2_sparse_238vols','Series_016_cmrr_mbep2d_3x2_sparse_238vols','Series_018_cmrr_mbep2d_3x2_sparse_238vols','Series_011_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_013_cmrr_mbep2d_3x2_sparse_invPE_SBRef'}; 
blocksin{cnt} = {'DATA_0020.nii','DATA_0021.nii','DATA_0010.nii','DATA_0012.nii','DATA_0016.nii','DATA_0018.nii','DATA_0011.nii','DATA_0013.nii'}; 
blocksout{cnt} = {'structural','INV2','Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'}; 
minvols(cnt) = 238; 
group(cnt) = 1; 

cnt = cnt + 1; 
subjects{cnt} = 'P7C17'; 
dates{cnt} = '20201105'; 
fullid{cnt} = '28080/20201105_U-ID50312'; 
basedir{cnt} = 'PINFA'; 
blocksin_folders{cnt} = {'Series_022_mp2rage_sag_p3_0.75mm_UNI_Images','Series_023_mp2rage_sag_p3_0.75mm_INV2','Series_010_cmrr_mbep2d_3x2_sparse_238vols','Series_012_cmrr_mbep2d_3x2_sparse_238vols','Series_016_cmrr_mbep2d_3x2_sparse_238vols','Series_020_cmrr_mbep2d_3x2_sparse_238vols','Series_011_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_013_cmrr_mbep2d_3x2_sparse_invPE_SBRef'}; 
blocksin{cnt} = {'DATA_0022.nii','DATA_0023.nii','DATA_0010.nii','DATA_0012.nii','DATA_0016.nii','DATA_0020.nii','DATA_0011.nii','DATA_0013.nii'}; 
blocksout{cnt} = {'structural','INV2','Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'}; 
minvols(cnt) = 238; 
group(cnt) = 1; 

cnt = cnt + 1; 
subjects{cnt} = 'P7C18'; 
dates{cnt} = '20201106'; 
fullid{cnt} = '28059/20201106_U-ID50598'; 
basedir{cnt} = 'PINFA'; 
blocksin_folders{cnt} = {'Series_021_mp2rage_sag_p3_0.75mm_UNI_Images','Series_020_mp2rage_sag_p3_0.75mm_INV2','Series_010_cmrr_mbep2d_3x2_sparse_238vols','Series_012_cmrr_mbep2d_3x2_sparse_238vols','Series_016_cmrr_mbep2d_3x2_sparse_238vols','Series_018_cmrr_mbep2d_3x2_sparse_238vols','Series_011_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_013_cmrr_mbep2d_3x2_sparse_invPE_SBRef'}; 
blocksin{cnt} = {'DATA_0021.nii','DATA_0020.nii','DATA_0010.nii','DATA_0012.nii','DATA_0016.nii','DATA_0018.nii','DATA_0011.nii','DATA_0013.nii'}; 
blocksout{cnt} = {'structural','INV2','Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'}; 
minvols(cnt) = 238; 
group(cnt) = 1; 

cnt = cnt + 1; 
subjects{cnt} = 'P7C19'; 
dates{cnt} = '20201113'; 
fullid{cnt} = '30131/20201113_U-ID50653'; 
basedir{cnt} = 'PINFA'; 
blocksin_folders{cnt} = {'Series_023_mp2rage_sag_p3_0.75mm_UNI_Images','Series_022_mp2rage_sag_p3_0.75mm_INV2','Series_010_cmrr_mbep2d_3x2_sparse_238vols','Series_014_cmrr_mbep2d_3x2_sparse_238vols','Series_018_cmrr_mbep2d_3x2_sparse_238vols','Series_020_cmrr_mbep2d_3x2_sparse_238vols','Series_013_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_015_cmrr_mbep2d_3x2_sparse_invPE_SBRef'}; 
blocksin{cnt} = {'DATA_0023.nii','DATA_0022.nii','DATA_0010.nii','DATA_0014.nii','DATA_0018.nii','DATA_0020.nii','DATA_0013.nii','DATA_0015.nii'}; 
blocksout{cnt} = {'structural','INV2','Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'}; 
minvols(cnt) = 238; 
group(cnt) = 1; 

cnt = cnt + 1; 
subjects{cnt} = 'P7C20'; 
dates{cnt} = '20201120'; 
fullid{cnt} = '26664/20201120_U-ID50599'; 
basedir{cnt} = 'PINFA'; 
blocksin_folders{cnt} = {'Series_021_mp2rage_sag_p3_0.75mm_UNI_Images','Series_020_mp2rage_sag_p3_0.75mm_INV2','Series_010_cmrr_mbep2d_3x2_sparse_238vols','Series_012_cmrr_mbep2d_3x2_sparse_238vols','Series_016_cmrr_mbep2d_3x2_sparse_238vols','Series_018_cmrr_mbep2d_3x2_sparse_238vols','Series_011_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_013_cmrr_mbep2d_3x2_sparse_invPE_SBRef'}; 
blocksin{cnt} = {'DATA_0021.nii','DATA_0020.nii','DATA_0010.nii','DATA_0012.nii','DATA_0016.nii','DATA_0018.nii','DATA_0011.nii','DATA_0013.nii'}; 
blocksout{cnt} = {'structural','INV2','Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'}; 
minvols(cnt) = 238; 
group(cnt) = 1; 

cnt = cnt + 1; 
subjects{cnt} = 'P7C21'; 
dates{cnt} = '20201126'; 
fullid{cnt} = '25813/20201126_U-ID50655'; 
basedir{cnt} = 'PINFA'; 
blocksin_folders{cnt} = {'Series_021_mp2rage_sag_p3_0.75mm_UNI_Images','Series_020_mp2rage_sag_p3_0.75mm_INV2','Series_010_cmrr_mbep2d_3x2_sparse_238vols','Series_012_cmrr_mbep2d_3x2_sparse_238vols','Series_016_cmrr_mbep2d_3x2_sparse_238vols','Series_018_cmrr_mbep2d_3x2_sparse_238vols','Series_011_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_013_cmrr_mbep2d_3x2_sparse_invPE_SBRef'}; 
blocksin{cnt} = {'DATA_0021.nii','DATA_0020.nii','DATA_0010.nii','DATA_0012.nii','DATA_0016.nii','DATA_0018.nii','DATA_0011.nii','DATA_0013.nii'}; 
blocksout{cnt} = {'structural','INV2','Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'}; 
minvols(cnt) = 238; 
group(cnt) = 1; 

cnt = cnt + 1; 
subjects{cnt} = 'P7C22'; 
dates{cnt} = '20201130'; 
fullid{cnt} = '31388/20201130_U-ID50678'; 
basedir{cnt} = 'PINFA'; 
blocksin_folders{cnt} = {'Series_020_mp2rage_sag_p3_0.75mm_UNI_Images','Series_021_mp2rage_sag_p3_0.75mm_INV2','Series_010_cmrr_mbep2d_3x2_sparse_238vols','Series_012_cmrr_mbep2d_3x2_sparse_238vols','Series_016_cmrr_mbep2d_3x2_sparse_238vols','Series_018_cmrr_mbep2d_3x2_sparse_238vols','Series_011_cmrr_mbep2d_3x2_sparse_238vols_SBRef','Series_013_cmrr_mbep2d_3x2_sparse_invPE_SBRef'}; 
blocksin{cnt} = {'DATA_0020.nii','DATA_0021.nii','DATA_0010.nii','DATA_0012.nii','DATA_0016.nii','DATA_0018.nii','DATA_0011.nii','DATA_0013.nii'}; 
blocksout{cnt} = {'structural','INV2','Run_1','Run_2','Run_3','Run_4','Pos_topup','Neg_topup'}; 
minvols(cnt) = 238; 
group(cnt) = 1; 

%% Add here any bad EPI runs for exclusion, see check_epi_quality.m for a way to look at this from your first level SPM
Bad_run_subj_pairs = {'P7P16' , 4;
'P7C17', 4
'P7C18', 3}; %P7P01 has a frontal dropout on all scans, but probably outside of area of interest, so let him stand.

for i = 1:length(subjects)
    if ~strcmp(subjects{i},Bad_run_subj_pairs(:,1))
        continue
    else
        TF = strcmp(subjects{i},Bad_run_subj_pairs(:,1));
        theseepis = find(strncmp(blocksout{i},'Run',3));
        badepi = theseepis(Bad_run_subj_pairs{TF,2});
        blocksout{i}(badepi) = [];
        blocksin{i}(badepi) = [];
        blocksin_folders{i}(badepi) = [];
    end
end

