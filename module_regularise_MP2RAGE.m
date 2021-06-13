function module_regularise_MP2RAGE(UNI_file,INV1_file,INV2_file,Out_file)

% Module to implement regularisation as per
%
% O?Brien, et al, 2014. 
% Robust T1-Weighted Structural Brain Imaging and Morphometry at 7T Using 
% MP2RAGE. v
% PLOS ONE 9, e99676. doi:10.1371/journal.pone.0099676
%
% which allows the creation of MP2RAGE T1w images without the strong 
% background noise in air regions. 
%
% although in the original paper the method only worked on raw multichannel
% data, here that contrain has been overcome and the correction can be
% implemented if both SOS images of the two inversion times exist and a
% MP2RAGE T1w image that has been calculated directly from the multichannel
% data as initially proposed in Marques et al, Neuroimage, 2009
%

addpath(genpath('./MP2RAGE_denoise/'))

MP2RAGE.filenameUNI=UNI_file; % standard MP2RAGE T1w image;
MP2RAGE.filenameINV1=INV1_file;% Inversion Time 1 M-P2RAGE T1w image;
MP2RAGE.filenameINV2=INV2_file;% Inversion Time 2 MP2RAGE T1w image;
MP2RAGE.filenameOUT=Out_file;% image without background noise;

regularisation = 5;
RobustCombination(MP2RAGE,regularisation);