function module_freesurfer_hires(preprocessedpathstem,subjects,crun,rawdatafolder,scanname,scriptdir,overwrite,skullstripped)
%Recon-all with hires flag and mris_inflate -n 30 on already skullstripped
%structural with csf retained
% this_scan(crun) = cellstr([rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))} '/' blocksin{crun}{find(strcmp(blocksout{crun},'structural'))}]);
% this_scan(crun) = cellstr([preprocessedpathstem subjects{crun} '/structural_csf.nii']);

this_scan(crun) = cellstr([rawdatafolder '/' scanname]);
if skullstripped
   this_subjects_dir = [preprocessedpathstem '/freesurfer_skullstripped/']; 
else
    this_subjects_dir = [preprocessedpathstem '/freesurfer/']; 
end

%pathstem = [preprocessedpathstem subjects{crun}];
setenv('RAW_DATA_FOLDER',rawdatafolder);
if ~exist(this_subjects_dir)
    mkdir(this_subjects_dir);
elseif exist([this_subjects_dir subjects{crun}])&&overwrite
    rmdir([this_subjects_dir subjects{crun}],'s')
    %mkdir([this_subjects_dir subjects{crun}]);
end
setenv('SUBJECTS_DIR',this_subjects_dir);
setenv('FSF_OUTPUT_FORMAT','nii');

if skullstripped == 0
    %cmd = ['recon-all -all -s ' subjects{crun} ' -hires -i ' this_scan{crun} ' -notal-check -cw256 -expert ' scriptdir 'expert_hires.opts'];
    cmd = ['recon-all -all -s ' subjects{crun} ' -hires -i ' this_scan{crun} ' -notal-check -cw256'];
    fprintf(['Submitting the following command: ' cmd]);
    system(cmd);
else
    %error('Skullstripped version not yet implemented - see recipe below')
    % Q. I have already skull-stripped data. Can I submit it to recon-all?
    % A: If your skull-stripped volume does not have the cerebellum, then no. If it does, then yes, however you will have to run the data a bit differently.
    % First you must run only -autorecon1 like this:
    % recon-all -autorecon1 -noskullstrip -s <subjid>
    cmd = ['recon-all -autorecon1 -noskullstrip -s ' subjects{crun} ' -hires -i ' this_scan{crun} ' -notal-check -cw256'];
    fprintf(['Submitting the following first stage command: ' cmd]);
    system(cmd);
    
    % Then you will have to make a symbolic link or copy T1.mgz to brainmask.auto.mgz and a link from brainmask.auto.mgz to brainmask.mgz. Finally, open this brainmask.mgz file and check that it looks okay (there is no skull, cerebellum is intact; use the sample subject bert that comes with your FreeSurfer installation to make sure it looks comparable). From there you can run the final stages of recon-all:
    % recon-all -autrecon2 -autorecon3 -s <subjid>
    data_dir = [this_subjects_dir subjects{crun}];
    copyfile([data_dir '/mri/T1.mgz'],[data_dir '/mri/brainmask.auto.mgz'])
    copyfile([data_dir '/mri/T1.mgz'],[data_dir '/mri/brainmask.mgz'])
    cmd = ['recon-all -autrecon2 -autorecon3 -noskullstrip -s ' subjects{crun} ' -hires -notal-check -cw256'];
    fprintf(['Submitting the following second stage command: ' cmd]);
    system(cmd);
end