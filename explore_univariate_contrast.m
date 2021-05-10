function explore_univariate_contrast(subjects,preprocessedpathstem,this_smooth,all_conditions)
% 
%     all_conditions = {
%         'con_0005.nii','Match > Mismatch';
%         'con_0010.nii','Mismatch > Match';
%         'con_0015.nii','Normal > Written';
%         'con_0020.nii','Written > Normal';
%         'con_0025.nii','Normal > Silence';
%         'con_0030.nii','Clear > Unclear';
%         'con_0035.nii','Unclear > Clear';
%         'con_0040.nii','Clarity Congruency Interaction'};
    expected_sessions = 4;
   
    nrun = size(all_conditions,1); % enter the number of runs here
    %jobfile = {'/group/language/data/thomascope/vespa/SPM12version/Standalone preprocessing pipeline/tc_source/batch_forwardmodel_job_noheadpoints.m'};
    
    this_scan = {};
    this_t_scan = {};
    firstlevel_folder = ['stats4_' num2str(this_smooth)];
    
    jobfile = {'/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/module_explore_univariate_job.m'};
    jobs = repmat(jobfile, 1, nrun);
    inputs = {};
    
    for this_condition = 1:nrun
        group1_mrilist = {}; %NB: Patient MRIs, so here group 2 (sorry)
        group2_mrilist = {};
        
        for crun = 1:size(subjects,2)
            this_spm_temp = load([preprocessedpathstem subjects{crun} filesep firstlevel_folder filesep 'SPM.mat']);
            if length(this_spm_temp.SPM.Sess)~=expected_sessions
                disp([subjects{crun} ' has ' num2str(length(this_spm_temp.SPM.Sess)) ' sessions when ' num2str(expected_sessions) ' expected. Check this is what you want'])
                con_num = str2num(all_conditions{this_condition,1}(5:8));
                new_con_num = (con_num/(expected_sessions+1))*(length(this_spm_temp.SPM.Sess)+1);
                disp(['Replacing contrast ' num2str(con_num,'%04.f') ' with ' num2str(new_con_num,'%04.f')])
            else
                new_con_num = str2num(all_conditions{this_condition,1}(5:8));
            end
                inputs{crun}{1, this_condition} = cellstr([preprocessedpathstem subjects{crun} filesep firstlevel_folder filesep 'SPM.mat']);
                inputs{crun}{2, this_condition} = new_con_num;
        end

    end
    
    spm('defaults', 'fMRI');
    spm_jobman('initcfg')
    for crun = 1:size(subjects,2)
        for this_condition = 1:nrun
            spm_jobman('run', jobs{this_condition}, inputs{crun}{:,this_condition});
            pause
        end
    end
end
