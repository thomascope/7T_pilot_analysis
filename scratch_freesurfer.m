parfor crun = 1:nrun
    rawdatafolder = [rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))}];
    scanname = blocksin{crun}{find(strcmp(blocksout{crun},'structural'))};
    %rawdatafolder = [preprocessedpathstem subjects{crun}];
    %scanname = 'structural_csf.nii';
    overwrite = 1
    setenv('FREESURFER_HOME','/home/tc02/freesurfer'); %Freesurfer 7.1.0
    !csh /home/tc02/freesurfer/SetUpFreeSurfer.csh
    
    module_freesurfer_hires(preprocessedpathstem,subjects,crun,rawdatafolder,scanname,scriptdir,overwrite)
end