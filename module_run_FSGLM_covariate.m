function module_run_FSGLM_covariate(datafolder, group_names, studyname, fsstats_already_done, view_data, smoothing_kernel, covariates)
% Runs the FSGLM using previously specified FSGD and contrasts (see
% module_make_FSGD.m)
setenv('SUBJECTS_DIR',datafolder);
if ~exist([datafolder '/glm'])
    mkdir([datafolder '/glm'])
end
this_dir = pwd;
cd([datafolder '/glm'])
hemispheres = {'lh','rh'};
glm_types = {
    'grouped','Grouped'
    'ungrouped','Ungrouped'
    'patients','Patients'
    };

if ~fsstats_already_done
    % First resample each hemisphere into average space
    parfor l = 1:2*(size(glm_types,1))
        i=mod(l,2)+1;
        k = ceil(l/2);
        setenv('SUBJECTS_DIR',datafolder);
        cmd = ['mris_preproc --fsgd ' datafolder '/FSGD/' studyname '_' glm_types{k,1} '.fsgd --target fsaverage --hemi ' hemispheres{i} ' --meas thickness --out  ' hemispheres{i} '.' studyname '_' glm_types{k,1} '.thickness.00.mgh']
        system(cmd)
    end
    
    % Next smooth the data - here with 10mm kernel, cortex only
    parfor l = 1:2*(size(glm_types,1))
        i=mod(l,2)+1;
        k = ceil(l/2);
        setenv('SUBJECTS_DIR',datafolder);
        cmd = ['mri_surf2surf --hemi ' hemispheres{i} ' --s fsaverage --sval ' hemispheres{i} '.' studyname '_' glm_types{k,1} '.thickness.00.mgh --fwhm ' num2str(smoothing_kernel) ' --cortex --tval ' hemispheres{i} '.' studyname '_' glm_types{k,1} '.thickness.' num2str(smoothing_kernel) '.mgh']
        
        system(cmd)
    end
end

% Next fit the GLM contrasts previously specified
for j = 1:size(covariates,2)
    parfor l = 1:2*(size(glm_types,1))
        i=mod(l,2)+1;
        k = ceil(l/2);
        setenv('SUBJECTS_DIR',datafolder);
        cmd = ['mri_glmfit --y  ' hemispheres{i} '.' studyname '_' glm_types{k,1} '.thickness.' num2str(smoothing_kernel) '.mgh --fsgd ' datafolder '/FSGD/' studyname '_' glm_types{k,1} '.fsgd doss --C ' datafolder '/Contrast/' glm_types{k,2} '_Covariate_' num2str(j) '.mtx --surf fsaverage ' hemispheres{i} ' --cortex --glmdir ' hemispheres{i} '.' glm_types{k,2} '_Covariate_' num2str(j) '.' num2str(smoothing_kernel) '.glmdir --eres-save']
        system(cmd)
    end
    % Optional view uncorrected map at generous threshold
    if view_data
        for l = 1:2*(size(glm_types,1))
            i=mod(l,2)+1;
            k = ceil(l/2);
            setenv('SUBJECTS_DIR',datafolder);
            cmd = ['freeview -f $SUBJECTS_DIR/fsaverage/surf/' hemispheres{i} '.inflated:annot=aparc.annot:annot_outline=1:overlay=' hemispheres{i} '.' glm_types{k,2} '_Covariate_' num2str(j) '.' num2str(smoothing_kernel) '.glmdir/' glm_types{k,2} '_Covariate_' num2str(j) '/sig.mgh:overlay_threshold=2,5 -viewport 3d'];
            system(cmd)
        end
    end
    % Now do clusterwise statistics
    parfor l = 1:2*(size(glm_types,1))
        i=mod(l,2)+1;
        k = ceil(l/2);
        setenv('SUBJECTS_DIR',datafolder);
        nperm = 10000;
        %nperm = 100; %Quick for code check
        cmd = ['mri_glmfit-sim --glmdir ' hemispheres{i} '.' glm_types{k,2} '_Covariate_' num2str(j) '.' num2str(smoothing_kernel) '.glmdir --perm ' nperm ' 4 pos --cwp 0.99 --bg 1 --overwrite'] % Quick for code check
        system(cmd)
    end
    if view_data
        for l = 1:2*(size(glm_types,1))
            i=mod(l,2)+1;
            k = ceil(l/2);
            setenv('SUBJECTS_DIR',datafolder);
            cmd = ['freeview -f $SUBJECTS_DIR/fsaverage/surf/' hemispheres{i} '.inflated:overlay=' hemispheres{i} '.' glm_types{k,2} '_Covariate_' num2str(j) '.' num2str(smoothing_kernel) '.glmdir/' glm_types{k,2} '_Covariate_' num2str(j) '/perm.th40.pos.sig.cluster.mgh:overlay_threshold=2,5:annot=' hemispheres{i} '.' glm_types{k,2} '_Covariate_' num2str(j) '.' num2str(smoothing_kernel) '.glmdir/' glm_types{k,2} '_Covariate_' num2str(j) '/perm.th40.pos.sig.ocn.annot -viewport 3d -layout 1'];
            system(cmd)
        end
    end
end

cd(this_dir)