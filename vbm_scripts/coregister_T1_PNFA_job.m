%-----------------------------------------------------------------------
% Job saved on 06-Sep-2016 13:01:50 by cfg_util (rev $Rev: 6460 $)
% spm SPM - SPM12 (6685)
%-----------------------------------------------------------------------
matlabbatch{1}.spm.spatial.coreg.estimate.ref = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/p00259/elderly/mean_T1_elderly_demented_brain.nii,1'};
matlabbatch{1}.spm.spatial.coreg.estimate.source = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/p00259/con_10075/mprage_125/con_10075_mprage_125.nii,1'};
matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
