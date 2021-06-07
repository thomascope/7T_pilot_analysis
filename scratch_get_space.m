for tissue_class = 1:3
    M = spm_get_space(scanname);
    copyfile(['c' num2str(tissue_class) scanname],['test_coreg_c' num2str(tissue_class) scanname]);
    spm_get_space(['test_coreg_c' num2str(tissue_class) scanname],M);
    all_coreged_vols = [all_coreged_vols; 'test_coreg_c' num2str(tissue_class) scanname];
end
Vo = spm_imcalc(all_coreged_vols,'test_coreg_brainmask.nii','i1+i2+i3>0');
temp = spm_read_vols(spm_vol(Vo));
filled_mask = imfill(temp,'holes');
Vo.fname = 'test_coreg_filled_brainmask.nii';
spm_write_vol(Vo,filled_mask);
spm_imcalc(char({scanname;'test_coreg_filled_brainmask.nii'}),'test_coreg_filled_skullstripped.nii','i1.*i2');

