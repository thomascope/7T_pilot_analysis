
coreg_and_crop_images = 0; %Necessary, but can skip if you've done before
overwrite = 1; % If you want to re-do a step - otherwise will crash if data already exist.
flags.which = 1; % For reslicing
view_scans = 0; % Obviously doesn't work in parallel
parfor crun = 1:nrun
    setenv('FREESURFER_HOME','/home/tc02/freesurfer'); %Freesurfer 7.1.0
    !csh /home/tc02/freesurfer/SetUpFreeSurfer.csh
    setenv('FSF_OUTPUT_FORMAT','nii');
    if skullstripped == 0 % Let Freesurfer skullstrip - fails more than a third of the time
        rawdatafolder = [rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))}];
        scanname = blocksin{crun}{find(strcmp(blocksout{crun},'structural'))};
        this_scan = cellstr([rawdatafolder '/' scanname]);
        if coreg_and_crop_images
            addpath('./crop')
            crop_images(this_scan) % Vital to do rough Talaraiching first as Freesurfer can't manage it.
        end
        scanname = ['p' scanname];
        
    elseif skullstripped == 1 % Use the SPM skullstripped image with csf as the input - works every time but has holes in the white matter that interfere with temporal lobe.
        rawdatafolder = [preprocessedpathstem subjects{crun}];
        scanname = 'structural_csf.nii';
        this_scan = cellstr([rawdatafolder '/' scanname]);
        if coreg_and_crop_images
            addpath('./crop')
            crop_images(this_scan) % Vital to do rough Talaraiching first as Freesurfer can't manage it.
        end
        scanname = ['p' scanname];
        
    elseif skullstripped == 2 % Create an imfilled mask for freesurfer from the SPM segmentation, which must have been run previously
        rawdatafolder = [rawpathstem basedir{crun} '/' fullid{crun} '/' blocksin_folders{crun}{find(strcmp(blocksout{crun},'structural'))}];
        scanname = blocksin{crun}{find(strcmp(blocksout{crun},'structural'))};
        this_scan = cellstr([rawdatafolder '/' scanname]);
        if coreg_and_crop_images
            addpath('./crop')
            crop_images(this_scan) % Vital to do rough Talaraiching first as Freesurfer can't manage it.
        end
        
        this_dir = pwd;
        cd(rawdatafolder)
        all_coreged_vols = [];
        for tissue_class = 1:3
            vol1 = spm_vol(char(this_scan));
            vol2 = spm_vol(['c' num2str(tissue_class) scanname]);
            new_vol = vol2;
            new_vol.fname = ['coreg_' vol2.fname];
            new_vol.mat = vol1.mat;
            spm_write_vol(new_vol,spm_read_vols(vol2));
            all_coreged_vols = [all_coreged_vols; 'coreg_c' num2str(tissue_class) scanname];
        end
        spm_reslice(char([this_scan; all_coreged_vols]),flags);
        rs = char(ones(size(all_coreged_vols,1),1) * 'r');
        all_resliced_coreged_vols = [rs, all_coreged_vols];
        Vo = spm_imcalc(all_resliced_coreged_vols,'coreg_brainmask.nii','i1+i2+i3>0');
        Vo = spm_imcalc(all_resliced_coreged_vols,'coreg_filled_brainmask.nii','imfill(i1+i2+i3>0,''holes'')');
        % % Doing it this way led to the images being written to one side -
        % % I don't understand why
        %         temp = spm_read_vols(spm_vol(Vo));
        %         filled_mask = imfill(temp,'holes');
        %         spm_write_vol(spm_vol('coreg_brainmask.nii'),filled_mask);
        %         spm_imcalc(char({scanname;'coreg_brainmask.nii'}),'coreg_filled_skullstripped.nii','i1.*i2');
        Vo = spm_imcalc(char([all_resliced_coreged_vols;this_scan]),'coreg_filled_skullstripped.nii','imfill(i1+i2+i3>0,''holes'').*i4');
        scanname = 'coreg_filled_skullstripped.nii';
        if view_scans
            %         % Optional checks - highly recommended as this step has gone all
            %         % kinds of wrong for reasons I don't understand
            spm_check_registration(char([this_scan; all_coreged_vols; all_resliced_coreged_vols; 'coreg_brainmask.nii'; 'coreg_filled_brainmask.nii'; 'coreg_filled_skullstripped.nii']))
            system(['freeview -v coreg_filled_skullstripped.nii'])
        end
        cd(this_dir)
    end
    module_freesurfer_hires(preprocessedpathstem,subjects,crun,rawdatafolder,scanname,scriptdir,overwrite,skullstripped)
    
end

% Now manually view outputs and decide which ones are bad
bad_subjects = {};
% output_folder = [preprocessedpathstem '/freesurfer/']; % For skullstripped == 0
% bad_subjects = {
%     'P7P01'
%     'P7P02'
%     'P7C02'
%     'P7P06'
%     'P7P08'
%     'P7P09'
%     'P7C08'
%     'P7P11'
%     'P7C11'
%     'P7C12'
%     'P7P13'
%     'P7P16'
%     'P7C16'
%     'P7C17'
%     'P7C22'
%     }; % With raw data
% 
% bad_subjects = {}; %None with either of the other methods, but temporal
% %lobes often wrong with skullstripped - masked is better
% output_folder = [preprocessedpathstem '/freesurfer_skullstripped/']; % For skullstripped == 1
output_folder = [preprocessedpathstem '/freesurfer_masked/']; % For skullstripped == 2
if isempty(bad_subjects)
    for crun = 1:nrun
        cmd = ['freeview -v ' output_folder subjects{crun} '/mri/T1.mgz ' output_folder subjects{crun} '/mri/wm.mgz ' output_folder subjects{crun} '/mri/brainmask.mgz ' output_folder subjects{crun} '/mri/aseg.mgz:colormap=lut:opacity=0.2 -f ' output_folder subjects{crun} '/surf/lh.white:edgecolor=blue ' output_folder subjects{crun} '/surf/lh.pial:edgecolor=red ' output_folder subjects{crun} '/surf/rh.white:edgecolor=blue ' output_folder subjects{crun} '/surf/rh.pial:edgecolor=red'];
        system(cmd)
        isgood(crun) = input('Does this subject look good? y or n.');
    end
    bad_subjects = subjects(isgood=='n')';
end

