parfor crun = 1:nrun
   coreg_and_crop_images = 0; %Necessary, but can skip if you've done before
%     rawdatafolder = [rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))}];
%     scanname = blocksin{crun}{find(strcmp(blocksout{crun},'structural'))};
%     skullstripped = 0;
    rawdatafolder = [preprocessedpathstem subjects{crun}];
    scanname = 'structural_csf.nii';
    skullstripped = 1;
    overwrite = 0
    setenv('FREESURFER_HOME','/home/tc02/freesurfer'); %Freesurfer 7.1.0
    !csh /home/tc02/freesurfer/SetUpFreeSurfer.csh
    
    this_scan = cellstr([rawdatafolder '/' scanname]);
    if coreg_and_crop_images
        addpath('./crop')
        crop_images(this_scan) % Vital to do rough Talaraiching first as Freesurfer can't manage it.
    end
    scanname = ['p' scanname];
    
    module_freesurfer_hires(preprocessedpathstem,subjects,crun,rawdatafolder,scanname,scriptdir,overwrite,skullstripped)
    
end

% Now manually view outputs and decide which ones are bad
output_folder = [preprocessedpathstem '/freesurfer/'];
bad_subjects = {    
    'P7P01'
    'P7P02'
    'P7C02'
    'P7P06'
    'P7P08'
    'P7P09'
    'P7C08'
    'P7P11'
    'P7C11'
    'P7C12'
    'P7P13'
    'P7P16'
    'P7C16'
    'P7C17'
    'P7C22'
    };

if isempty(bad_subjects)
    for crun = 1:nrun
        cmd = ['freeview -v ' output_folder subjects{crun} '/mri/T1.mgz ' output_folder subjects{crun} '/mri/wm.mgz ' output_folder subjects{crun} '/mri/brainmask.mgz ' output_folder subjects{crun} '/mri/aseg.mgz:colormap=lut:opacity=0.2 -f ' output_folder subjects{crun} '/surf/lh.white:edgecolor=blue ' output_folder subjects{crun} '/surf/lh.pial:edgecolor=red ' output_folder subjects{crun} '/surf/rh.white:edgecolor=blue ' output_folder subjects{crun} '/surf/rh.pial:edgecolor=red'];
        system(cmd)
        isgood(crun) = input('Does this subject look good? y or n.');
    end
    bad_subjects = subjects(isgood=='n')';
end

