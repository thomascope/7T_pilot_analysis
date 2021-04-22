% List of open inputs
% Coregister: Estimate: Source Image - cfg_files

T1_list = {
'./p00259/con_16672/mprage_125/con_01_16672_mprage_125.nii'
'./p00259/ppa_23503/mprage_125/ppa_23503_mprage_125.nii'
'./p00259/con_16255/mprage_125/con_16255_mprage_125.nii'
'./p00259/con_14863/mprage_125/con_14863_mprage_125.nii'
'./p00259/con_14844/mprage_125/con_14844_mprage_125.nii'
'./p00259/ppa_20218/mprage_125/ppa_20218_mprage_125.nii'
'./p00259/con_10075/mprage_125/con_10075_mprage_125.nii'
'./p00259/ppa_19727/mprage_125/ppa_19727_mprage_125.nii'
'./p00259/con_15581/mprage_125/con_15581_mprage_125.nii'
'./p00259/con_16254/mprage_125/con_16254_mprage_125.nii'
'./p00259/ppa_23512/mprage_125/ppa_23512_mprage_125.nii'
'./p00259/con_14847/mprage_125/con_14847_mprage_125.nii'
'./p00259/ppa_23416/mprage_125/ppa_23416_mprage_125.nii'
'./p00259/con_14848/mprage_125/con_14848_mprage_125.nii'
'./p00259/con_19617/mprage_125/con_19617_mprage_125.nii'
'./p00259/con_14855/mprage_125/con_14855_mprage_125.nii'
'./p00259/con_14827/mprage_125/con_14827_mprage_125.nii'
'./p00259/con_14842/mprage_125/con_14842_mprage_125.nii'
'./p00259/con_14854/mprage_125/con_14854_mprage_125.nii'
'./p00259/con_14826/mprage_125/con_14826_mprage_125.nii'
'./p00259/con_16258/mprage_125/con_16258_mprage_125.nii'
'./p00259/ppa_23483/mprage_125/ppa_23483_mprage_125.nii'
'./p00259/con_14869/mprage_125/con_14869_mprage_125.nii'
'./p00259/con_23705/mprage_125/con_23705_mprage_125.nii'
'./p00259/con_14843/mprage_125/con_14843_mprage_125.nii'
'./p00259/con_14867/mprage_125/con_14867_mprage_125.nii'
};

nrun = length(T1_list); % enter the number of runs here
jobfile = {'/imaging/tc02/vespa/scans/PNFA_VBM/tom/coregister_T1_blanked_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(1, nrun);

for crun = 1:nrun
    inputs{1, crun} = cellstr(T1_list{crun}); % for co_register_estimate
end

coregisterworkedcorrectly = zeros(1,nrun);
jobs = repmat(jobfile, 1, 1);

cbupool(nrun)

parfor crun = 1:nrun
    spm('defaults', 'PET');
    spm_jobman('initcfg')
    try
        spm_jobman('run', jobs, inputs{:,crun});
        coregisterworkedcorrectly(crun) = 1;
    catch
        coregisterworkedcorrectly(crun) = 0;
    end
end

% matlabpool close