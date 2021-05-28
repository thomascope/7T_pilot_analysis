function module_freesurfer_hires(preprocessedpathstem,subjects,crun,rawdatafolder,scanname,scriptdir,overwrite,skullstripped)
%Recon-all with hires flag

this_scan(crun) = cellstr([rawdatafolder '/' scanname]);
if skullstripped == 1
   this_subjects_dir = [preprocessedpathstem '/freesurfer_skullstripped/']; 
elseif skullstripped == 0
    this_subjects_dir = [preprocessedpathstem '/freesurfer/']; 
elseif skullstripped == 2
    this_subjects_dir = [preprocessedpathstem '/freesurfer_masked/']; 
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

if skullstripped == 0 % Let Freesurfer skullstrip - default but fails more than a third of the time with PINFA 7T dataset
    %cmd = ['recon-all -all -s ' subjects{crun} ' -hires -i ' this_scan{crun} ' -notal-check -cw256 -expert ' scriptdir 'expert_hires.opts'];
    cmd = ['recon-all -all -s ' subjects{crun} ' -hires -i ' this_scan{crun} ' -notal-check -cw256'];
    fprintf(['Submitting the following command: ' cmd]);
    system(cmd);
elseif skullstripped == 1 % Use the SPM skullstripped image with csf as the input - works every time but has holes in the white matter that interfere with temporal lobe.
    % Q. I have already skull-stripped data. Can I submit it to recon-all?
    % A: If your skull-stripped volume does not have the cerebellum, then no. If it does, then yes, however you will have to run the data a bit differently.
    % First you must run only -autorecon1 like this:
    % recon-all -autorecon1 -noskullstrip -s <subjid>
    cmd = ['recon-all -autorecon1 -noskullstrip -s ' subjects{crun} ' -hires -i ' this_scan{crun} ' -notal-check -cw256'];
    fprintf(['Submitting the following first stage command: ' cmd]);
    system(cmd);
    
    % Then you will have to make a symbolic link or copy T1.mgz to brainmask.auto.mgz and a link from brainmask.auto.mgz to brainmask.mgz. Finally, open this brainmask.mgz file and check that it looks okay (there is no skull, cerebellum is intact; use the sample subject bert that comes with your FreeSurfer installation to make sure it looks comparable). From there you can run the final stages of recon-all:
    % recon-all -autorecon2 -autorecon3 -s <subjid>
    data_dir = [this_subjects_dir subjects{crun}];
    copyfile([data_dir '/mri/T1.mgz'],[data_dir '/mri/brainmask.auto.mgz'])
    copyfile([data_dir '/mri/T1.mgz'],[data_dir '/mri/brainmask.mgz'])
    cmd = ['recon-all -autorecon2 -autorecon3 -noskullstrip -s ' subjects{crun} ' -hires -notal-check -cw256'];
    fprintf(['Submitting the following second stage command: ' cmd]);
    system(cmd);
elseif skullstripped == 2 % Create an imfilled mask for freesurfer from the SPM segmentation, which must have been run previously, and use this to skullstrip the input
    % Q. I have already skull-stripped data. Can I submit it to recon-all?
    % A: If your skull-stripped volume does not have the cerebellum, then no. If it does, then yes, however you will have to run the data a bit differently.
    % First you must run only -autorecon1 like this:
    % recon-all -autorecon1 -noskullstrip -s <subjid>
    cmd = ['recon-all -autorecon1 -noskullstrip -s ' subjects{crun} ' -hires -i ' this_scan{crun} ' -notal-check -cw256 -bigventricles'];
    fprintf(['Submitting the following first stage command: ' cmd]);
    system(cmd);
    
    % Then you will have to make a symbolic link or copy T1.mgz to brainmask.auto.mgz and a link from brainmask.auto.mgz to brainmask.mgz. Finally, open this brainmask.mgz file and check that it looks okay (there is no skull, cerebellum is intact; use the sample subject bert that comes with your FreeSurfer installation to make sure it looks comparable). From there you can run the final stages of recon-all:
    % recon-all -autorecon2 -autorecon3 -s <subjid>
    data_dir = [this_subjects_dir subjects{crun}];
    copyfile([data_dir '/mri/T1.mgz'],[data_dir '/mri/brainmask.auto.mgz'])
    copyfile([data_dir '/mri/T1.mgz'],[data_dir '/mri/brainmask.mgz'])
    cmd = ['recon-all -autorecon2 -autorecon3 -noskullstrip -s ' subjects{crun} ' -hires -notal-check -cw256 -bigventricles'];
    fprintf(['Submitting the following second stage command: ' cmd]);
    system(cmd);
else
    error('Unknown skullstrip argument specified')
end