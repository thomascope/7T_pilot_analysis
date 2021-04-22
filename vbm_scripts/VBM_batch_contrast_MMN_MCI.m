%-----------------------------------------------------------------------
% Job saved on 08-Sep-2016 10:08:43 by cfg_util (rev $Rev: 6460 $)
% spm SPM - SPM12 (6685)
%-----------------------------------------------------------------------
matlabbatch{1}.spm.stats.con.spmmat = '<UNDEFINED>';
% Group order {'matched_HCs' 'pca' 'bvFTD' 'pnfa'};
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Controls > pca';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'Controls > bvFTD';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [1 0 -1];
matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'Controls > nfvPPA';
matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = [1 0 0 -1];
matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'Controls < pca';
matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = [-1 1];
matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = 'Controls < bvFTD';
matlabbatch{1}.spm.stats.con.consess{5}.tcon.weights = [-1 0 1];
matlabbatch{1}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{6}.tcon.name = 'Controls < nfvPPA';
matlabbatch{1}.spm.stats.con.consess{6}.tcon.weights = [-1 0 0 1];
matlabbatch{1}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{7}.tcon.name = 'pca > bvFTD';
matlabbatch{1}.spm.stats.con.consess{7}.tcon.weights = [0 1 -1];
matlabbatch{1}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{8}.tcon.name = 'pca < bvFTD';
matlabbatch{1}.spm.stats.con.consess{8}.tcon.weights = [0 -1 1];
matlabbatch{1}.spm.stats.con.consess{8}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{9}.tcon.name = 'pca > nfvPPA';
matlabbatch{1}.spm.stats.con.consess{9}.tcon.weights = [0 1 0 -1];
matlabbatch{1}.spm.stats.con.consess{9}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{10}.tcon.name = 'pca < nfvPPA';
matlabbatch{1}.spm.stats.con.consess{10}.tcon.weights = [0 -1 0 1];
matlabbatch{1}.spm.stats.con.consess{10}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{11}.tcon.name = 'bvFTD > nfvPPA';
matlabbatch{1}.spm.stats.con.consess{11}.tcon.weights = [0 0 1 -1];
matlabbatch{1}.spm.stats.con.consess{11}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{12}.tcon.name = 'bvFTD < nfvPPA';
matlabbatch{1}.spm.stats.con.consess{12}.tcon.weights = [0 0 -1 1];
matlabbatch{1}.spm.stats.con.consess{12}.tcon.sessrep = 'none';
% Added MCI - do these contrasts last so the visualisation scripts require
% minimal edit
matlabbatch{1}.spm.stats.con.consess{13}.tcon.name = 'Controls > ADMCI';
matlabbatch{1}.spm.stats.con.consess{13}.tcon.weights = [1 0 0 0 -1];
matlabbatch{1}.spm.stats.con.consess{13}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{14}.tcon.name = 'Controls < ADMCI';
matlabbatch{1}.spm.stats.con.consess{14}.tcon.weights = [-1 0 0 0 1];
matlabbatch{1}.spm.stats.con.consess{14}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{15}.tcon.name = 'pca > ADMCI';
matlabbatch{1}.spm.stats.con.consess{15}.tcon.weights = [0 1 0 0 -1];
matlabbatch{1}.spm.stats.con.consess{15}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{16}.tcon.name = 'pca < ADMCI';
matlabbatch{1}.spm.stats.con.consess{16}.tcon.weights = [0 -1 0 0 1];
matlabbatch{1}.spm.stats.con.consess{16}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{17}.tcon.name = 'bvFTD > ADMCI';
matlabbatch{1}.spm.stats.con.consess{17}.tcon.weights = [0 0 1 0 -1];
matlabbatch{1}.spm.stats.con.consess{17}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{18}.tcon.name = 'bvFTD < ADMCI';
matlabbatch{1}.spm.stats.con.consess{18}.tcon.weights = [0 0 -1 0 1];
matlabbatch{1}.spm.stats.con.consess{18}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{19}.tcon.name = 'nfvPPA > ADMCI';
matlabbatch{1}.spm.stats.con.consess{19}.tcon.weights = [0 0 0 1 -1];
matlabbatch{1}.spm.stats.con.consess{19}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{20}.tcon.name = 'nfvPPA < ADMCI';
matlabbatch{1}.spm.stats.con.consess{20}.tcon.weights = [0 0 0 -1 1];
matlabbatch{1}.spm.stats.con.consess{20}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 0;
