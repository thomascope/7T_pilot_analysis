function secondlevelworkedcorrectly = module_searchlight_interaction_hires(GLMDir,subjects,group,age_lookup,outpath,downsamp_ratio,conds_top,conds_bottom,cond_names)
% Normalise effect-maps to MNI template

do_smoothed_maps = 0;  % if you want to do smoothed maps change this
if ~exist('downsamp_ratio','var')
    downsamp_ratio = 1;
end

versionCurrent = 'spearman';

% Gather images for current subject
if downsamp_ratio == 1
    images = cellstr(spm_select('FPList', [GLMDir '/TDTcrossnobis/' versionCurrent '/'], '^whireseffect-map_.*.nii'));
else
    images = cellstr(spm_select('FPList', [GLMDir '/TDTcrossnobis_downsamp_' num2str(downsamp_ratio) '/' versionCurrent '/'], '^whireseffect-map_.*.nii'));
end

nrun = length(conds_top); % enter the number of runs here - if want to do smoothed as well
%jobfile = {'/group/language/data/thomascope/vespa/SPM12version/Standalone preprocessing pipeline/tc_source/batch_forwardmodel_job_noheadpoints.m'};

for i = 1:nrun
    for this_top = 1:length(conds_top{i})
        if downsamp_ratio == 1
            these_tops{i}(this_top) = cellstr(spm_select('FPList', [GLMDir '/TDTcrossnobis/' versionCurrent '/'], ['^whireseffect-map_' conds_top{i}{this_top} '.nii']));
            assert(logical(exist(these_tops{i}{this_top},'file')),['^whireseffect-map_' conds_top{i}{this_top} '.nii not found']);
        else
            these_tops{i}(this_top) = cellstr(spm_select('FPList', [GLMDir '/TDTcrossnobis_downsamp_' num2str(downsamp_ratio) '/' versionCurrent '/'], ['^whireseffect-map_' conds_top{i}{this_top} '.nii']));
            assert(logical(exist(these_tops{i}{this_top},'file')),['^whireseffect-map_' conds_top{i}{this_top} '.nii not found']);
        end
    end
    for this_bottom = 1:length(conds_bottom{i})
        if downsamp_ratio == 1
            these_bottoms{i}(this_bottom) = cellstr(spm_select('FPList', [GLMDir '/TDTcrossnobis/' versionCurrent '/'], ['^whireseffect-map_' conds_bottom{i}{this_bottom} '.nii']));
            assert(logical(exist(these_tops{i}{this_top},'file')),['^whireseffect-map_' conds_bottom{i}{this_bottom} '.nii not found']);
        else
            these_bottoms{i}(this_bottom) = cellstr(spm_select('FPList', [GLMDir '/TDTcrossnobis_downsamp_' num2str(downsamp_ratio) '/' versionCurrent '/'], ['^whireseffect-map_' conds_bottom{i}{this_bottom} '.nii']));
            assert(logical(exist(these_tops{i}{this_top},'file')),['^whireseffect-map_' conds_bottom{i}{this_bottom} '.nii not found']);
        end
    end
end
    
% Now organise inputs to suit Rik Henson's batch_spm_anova script for
% repeated measures
% The only complicated bit is organising imgfiles correctly, so here's an example:
%
%  imgfiles{1}{1} = ['mydir/grp1_sub1_con1.nii'; 'mydir/grp1_sub1_con2.nii'];
%  imgfiles{1}{2} = ['mydir/grp1_sub2_con1.nii'; 'mydir/grp1_sub2_con2.nii'];
%  imgfiles{2}{1} = ['mydir/grp2_sub1_con1.nii'; 'mydir/grp2_sub1_con2.nii'];
%  imgfiles{2}{2} = ['mydir/grp2_sub2_con1.nii'; 'mydir/grp2_sub2_con2.nii'];
%
% Actually, the organisation of the later-added user_regs is also a bit complex:
%
%  user_regs{1} = [[1:4]' rand(4,1)];       % 2 regressors for group 1 
%  user_regs{2} = [[1:4]' rand(4,1)];       % 2 regressors for group 2

this_scan = {};

myWrapper = @(x) exist(x, 'file');
clear S
for this_contrast = 1:nrun
    S{this_contrast}.imgfiles{1} = {}; %NB: Patient MRIs, so here group 2 (sorry)
    S{this_contrast}.user_regs{1} = [];
    S{this_contrast}.imgfiles{2} = {};
    S{this_contrast}.user_regs{2} = [];
    
    S{this_contrast}.outdir = [outpath filesep 'interactions' filesep cond_names{this_contrast}];
    
    group1_mrilist = {}; 
    group1_ages = [];
    group2_mrilist = {};
    group2_ages = [];
    
     for crun = 1:size(subjects,2)
        this_age = age_lookup.Age(strcmp(age_lookup.Study_ID,subjects{crun}));
        these_images=[strrep(these_tops{this_contrast},subjects{1},subjects{crun})'; strrep(these_bottoms{this_contrast},subjects{1},subjects{crun})'];
        
        if group(crun) == 1 % Controls
            S{this_contrast}.imgfiles{2}{end+1} = these_images;
            S{this_contrast}.user_regs{2}(end+1) = this_age;
            assert(all(cellfun(myWrapper,S{this_contrast}.imgfiles{2}{end})),['Missing input for subject ' num2str(crun)])
            S{this_contrast}.imgfiles{2}{end} = char(S{this_contrast}.imgfiles{2}{end});
        elseif group(crun) == 2 % Patients
            S{this_contrast}.imgfiles{1}{end+1} = these_images;
            S{this_contrast}.user_regs{1}(end+1) = this_age;
            assert(all(cellfun(myWrapper,S{this_contrast}.imgfiles{1}{end})),['Missing input for subject ' num2str(crun)])
            S{this_contrast}.imgfiles{1}{end} = char(S{this_contrast}.imgfiles{1}{end});
        end
     end
    S{this_contrast}.user_regs{1}=repmat(S{this_contrast}.user_regs{1}',size(S{this_contrast}.imgfiles{1}{end},1),1);
    S{this_contrast}.user_regs{2}=repmat(S{this_contrast}.user_regs{2}',size(S{this_contrast}.imgfiles{2}{end},1),1);
    S{this_contrast}.mask = '/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/control_majority_unsmoothed_mask_c1_thr0.05_cons0.8.img';
    
    
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Controls > chance';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [0 1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'Patients > chance';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [1 0];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'Controls > Patients';
matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [-1 1];
matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'Patients > Controls';
matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [1 -1];
matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = 'Patients + Controls';
matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [1 1];
matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{6}.tcon.name = 'Controls < chance';
matlabbatch{3}.spm.stats.con.consess{6}.tcon.weights = [0 -1];
matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{7}.tcon.name = 'Patients < chance';
matlabbatch{3}.spm.stats.con.consess{7}.tcon.weights = [-1 0];
matlabbatch{3}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{8}.tcon.name = 'Patients + Controls Negative';
matlabbatch{3}.spm.stats.con.consess{8}.tcon.weights = [-1 -1];
matlabbatch{3}.spm.stats.con.consess{8}.tcon.sessrep = 'none';
    
main_contrast = [ones(1,length(conds_top{this_contrast}))/length(conds_top{this_contrast}) -ones(1,length(conds_bottom{this_contrast}))/length(conds_bottom{this_contrast}) zeros(size(S{this_contrast}.user_regs{1},2))];

S{this_contrast}.contrasts{1}.c = [0*main_contrast 1*main_contrast];
S{this_contrast}.contrasts{1}.type = 'T';
S{this_contrast}.contrasts{1}.name =  'Controls > chance';

S{this_contrast}.contrasts{2}.c = [1*main_contrast 0*main_contrast];
S{this_contrast}.contrasts{2}.type = 'T';
S{this_contrast}.contrasts{2}.name =  'Patients > chance';

S{this_contrast}.contrasts{3}.c = [-1*main_contrast 1*main_contrast];
S{this_contrast}.contrasts{3}.type = 'T';
S{this_contrast}.contrasts{3}.name =  'Controls > Patients';

S{this_contrast}.contrasts{4}.c = [1*main_contrast -1*main_contrast];
S{this_contrast}.contrasts{4}.type = 'T';
S{this_contrast}.contrasts{4}.name =  'Patients > Controls';

S{this_contrast}.contrasts{5}.c = [1*main_contrast 1*main_contrast];
S{this_contrast}.contrasts{5}.type = 'T';
S{this_contrast}.contrasts{5}.name =  'Patients + Controls';

S{this_contrast}.contrasts{6}.c = [0*main_contrast -1*main_contrast];
S{this_contrast}.contrasts{6}.type = 'T';
S{this_contrast}.contrasts{6}.name =  'Controls < chance';

S{this_contrast}.contrasts{7}.c = [-1*main_contrast 0*main_contrast];
S{this_contrast}.contrasts{7}.type = 'T';
S{this_contrast}.contrasts{7}.name =  'Patients < chance';

S{this_contrast}.contrasts{8}.c = [-1*main_contrast -1*main_contrast];
S{this_contrast}.contrasts{8}.type = 'T';
S{this_contrast}.contrasts{8}.name =  'Patients + Controls Negative';
end

secondlevelworkedcorrectly = zeros(1,nrun);

parfor this_contrast = 1:nrun
    if exist(fullfile(S{this_contrast}.outdir,'SPM.mat'),'file')
        disp([fullfile(S{this_contrast}.outdir,'SPM.mat') ' already exists, delete it if you want to re-make, otherwise moving on.'])
    else
        try
            try mkdir(S{this_contrast}.outdir); end
            spm('defaults', 'fMRI');
            batch_spm_anova(S{this_contrast})
            secondlevelworkedcorrectly(this_contrast) = 1;
        catch
            secondlevelworkedcorrectly(this_contrast) = 0;
        end
    end
end