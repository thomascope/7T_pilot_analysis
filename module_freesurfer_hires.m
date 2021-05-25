function module_freesurfer_hires(preprocessedpathstem,subjects,crun,rawdatafolder,scanname,scriptdir,overwrite)
%Recon-all with hires flag and mris_inflate -n 30 on already skullstripped
%structural with csf retained
% this_scan(crun) = cellstr([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '/' blocksin{crun}{find(strcmp(blocksout{crun},'structural'))}]);
% this_scan(crun) = cellstr([preprocessedpathstem subjects{crun} '/structural_csf.nii']);

this_scan(crun) = cellstr([rawdatafolder '/' scanname]);

pathstem = [preprocessedpathstem subjects{crun}];
setenv('RAW_DATA_FOLDER',rawdatafolder);
if ~exist([pathstem '/freesurfer'])
    mkdir([pathstem '/freesurfer']);
elseif exist([pathstem '/freesurfer'])&&overwrite
    rmdir([pathstem '/freesurfer'],'s')
    mkdir([pathstem '/freesurfer']);
end
setenv('SUBJECTS_DIR',[pathstem '/freesurfer']);
setenv('FSF_OUTPUT_FORMAT','nii');

cmd = ['recon-all -all -s ' subjects{crun} ' -hires -i ' this_scan{crun} ' -expert ' scriptdir 'expert_hires.opts'];
fprintf(['Submitting the following command: ' cmd]);
system(cmd);
