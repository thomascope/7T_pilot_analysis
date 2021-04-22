%-----------------------------------------------------------------------
% Job saved on 06-Sep-2016 18:34:51 by cfg_util (rev $Rev: 6460 $)
% spm SPM - SPM12 (6685)
%-----------------------------------------------------------------------
matlabbatch{1}.spm.tools.dartel.mni_norm.template = '<UNDEFINED>';
matlabbatch{1}.spm.tools.dartel.mni_norm.data.subjs.flowfields = '<UNDEFINED>';
matlabbatch{1}.spm.tools.dartel.mni_norm.data.subjs.images = {};
matlabbatch{1}.spm.tools.dartel.mni_norm.vox = [1 1 1];
matlabbatch{1}.spm.tools.dartel.mni_norm.bb = [-90 -126 -72
                                               90 90 108];
matlabbatch{1}.spm.tools.dartel.mni_norm.preserve = 1;
matlabbatch{1}.spm.tools.dartel.mni_norm.fwhm = [8 8 8];
matlabbatch{2}.spm.stats.factorial_design.dir = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/Stats/factorial_full_group_vbm'};
matlabbatch{2}.spm.stats.factorial_design.des.t2.scans1 = '<UNDEFINED>';
matlabbatch{2}.spm.stats.factorial_design.des.t2.scans2 = '<UNDEFINED>';
matlabbatch{2}.spm.stats.factorial_design.des.t2.dept = 0;
matlabbatch{2}.spm.stats.factorial_design.des.t2.variance = 1;
matlabbatch{2}.spm.stats.factorial_design.des.t2.gmsca = 0;
matlabbatch{2}.spm.stats.factorial_design.des.t2.ancova = 0;
matlabbatch{2}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{2}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{2}.spm.stats.factorial_design.masking.tm.tmr.rthresh = 0.8;
matlabbatch{2}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{2}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{2}.spm.stats.factorial_design.globalc.g_mean = 1;
matlabbatch{2}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{2}.spm.stats.factorial_design.globalm.glonorm = 2;
matlabbatch{3}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{3}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{4}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{4}.spm.stats.con.consess{1}.tcon.name = 'Controls>PNFA';
matlabbatch{4}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
matlabbatch{4}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{4}.spm.stats.con.consess{2}.tcon.name = 'PNFA>Controls';
matlabbatch{4}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
matlabbatch{4}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{4}.spm.stats.con.delete = 0;
matlabbatch{5}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{5}.spm.stats.results.conspec(1).titlestr = '';
matlabbatch{5}.spm.stats.results.conspec(1).contrasts = 1;
matlabbatch{5}.spm.stats.results.conspec(1).threshdesc = 'FWE';
matlabbatch{5}.spm.stats.results.conspec(1).thresh = 0.05;
matlabbatch{5}.spm.stats.results.conspec(1).extent = 0;
matlabbatch{5}.spm.stats.results.conspec(1).conjunction = 1;
matlabbatch{5}.spm.stats.results.conspec(1).mask.none = 1;
matlabbatch{5}.spm.stats.results.conspec(2).titlestr = '';
matlabbatch{5}.spm.stats.results.conspec(2).contrasts = 2;
matlabbatch{5}.spm.stats.results.conspec(2).threshdesc = 'FWE';
matlabbatch{5}.spm.stats.results.conspec(2).thresh = 0.05;
matlabbatch{5}.spm.stats.results.conspec(2).extent = 0;
matlabbatch{5}.spm.stats.results.conspec(2).conjunction = 1;
matlabbatch{5}.spm.stats.results.conspec(2).mask.none = 1;
matlabbatch{5}.spm.stats.results.units = 1;
matlabbatch{5}.spm.stats.results.print = 'ps';
matlabbatch{5}.spm.stats.results.write.none = 1;
