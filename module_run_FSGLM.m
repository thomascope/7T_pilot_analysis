function module_run_FSGLM(datafolder, group_names, studyname, precached, view_data)
% Runs the FSGLM using previously specified FSGD and contrasts (see
% module_make_FSGD.m)
setenv('SUBJECTS_DIR',datafolder);
if ~exist([datafolder '/glm'])
    mkdir([datafolder '/glm'])
end
this_dir = pwd;
cd([datafolder '/glm'])
hemispheres = {'lh','rh'};
if~precached
      
      % First resample each hemisphere into average space
   parfor i = 1:2
       setenv('SUBJECTS_DIR',datafolder);
       cmd = ['mris_preproc --fsgd ' datafolder '/FSGD/' studyname '.fsgd --target fsaverage --hemi ' hemispheres{i} ' --meas thickness --out  ' hemispheres{i} '.' studyname '.thickness.00.mgh']
       system(cmd)
   end
   % Next smooth the data - here with 10mm kernel, cortex only
   parfor i = 1:2
       setenv('SUBJECTS_DIR',datafolder);
       cmd = ['mri_surf2surf --hemi ' hemispheres{i} ' --s fsaverage --sval ' hemispheres{i} '.' studyname '.thickness.00.mgh --fwhm 10 --cortex --tval ' hemispheres{i} '.' studyname '.thickness.10.mgh']
       system(cmd)
   end
  
   
else
    error('Not yet written for precached data, but this could save time')
end

% Next fit the GLM contrast previously specified
parfor i = 1:2
    setenv('SUBJECTS_DIR',datafolder);
    cmd = ['mri_glmfit --y  ' hemispheres{i} '.' studyname '.thickness.10.mgh --fsgd ' datafolder '/FSGD/' studyname '.fsgd doss --C ' datafolder '/Contrast/' group_names{1} '-' group_names{2} '.mtx --surf fsaverage ' hemispheres{i} ' --cortex --glmdir ' hemispheres{i} '.' group_names{1} '-' group_names{2} '.glmdir --eres-save']
    system(cmd)
end
% Optional view uncorrected map at generous threshold
if view_data
    for i = 1:2
        cmd = ['freeview -f $SUBJECTS_DIR/fsaverage/surf/' hemispheres{i} '.inflated:annot=aparc.annot:annot_outline=1:overlay=' hemispheres{i} '.' group_names{1} '-' group_names{2} '.glmdir/' group_names{1} '-' group_names{2} '/sig.mgh:overlay_threshold=2,5 -viewport 3d'];
        system(cmd)
    end
end
% Now do clusterwise statistics
parfor i = 1:2
    cmd = ['mri_glmfit-sim --glmdir ' hemispheres{i} '.' group_names{1} '-' group_names{2} '.glmdir --perm 1000 4 pos --cwp 0.05 --bg 1 --overwrite']
    %cmd = ['mri_glmfit-sim --glmdir ' hemispheres{i} '.' group_names{1} '-' group_names{2} '.glmdir --perm 100 4 pos --cwp 0.99 --bg 1 --overwrite'] % Quick for code check
    system(cmd)
end
if view_data
    for i = 1:2
        cmd = ['freeview -f $SUBJECTS_DIR/fsaverage/surf/' hemispheres{i} '.inflated:overlay=' hemispheres{i} '.' group_names{1} '-' group_names{2} '.glmdir/' group_names{1} '-' group_names{2} '/perm.th40.pos.sig.cluster.mgh:overlay_threshold=2,5:annot=' hemispheres{i} '.' group_names{1} '-' group_names{2} '.glmdir/' group_names{1} '-' group_names{2} '/perm.th40.pos.sig.ocn.annot -viewport 3d -layout 1']
        system(cmd)
    end
end

cd(this_dir)