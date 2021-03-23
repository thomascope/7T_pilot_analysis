function []=GLMAnalysis_1Subj(ProcDataDir,GLMAnalDir,TempPD,SubjID,NumSess,SpmPD,ScriptsPD)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load SPM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(ScriptsPD); % So that any customised code in project scripts folder is prioritised (like spm_my_defaults)
addpath(SpmPD); % So that correct version of SPM is loaded
spm('FMRI');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model Specify
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load([TempPD '/GLMSpecify.mat']);

matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(GLMAnalDir);

matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2.5;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 32;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 16;

matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];


for i=1:NumSess
    SessEpiDir = fullfile(ProcDataDir,'Functional',['Sess_' num2str(i)]);
    SessFunctionals = cellstr(spm_select('FPList', SessEpiDir, '^uf.*\.nii$'));
    assertSizeAtLeast(1, SessFunctionals, ['Specify: Not expected number of files' num2str(i)]);
   
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).scans = SessFunctionals;
    RegFile=cellstr(spm_select('FPList',GLMAnalDir,  ['Regressors_Sess_' num2str(i) '_' SubjID '.mat']));
    assertSizeAtLeast(1, RegFile, ['Regressor file not found for subject_' SubjID]);
    
    MotionPhysioRegFile = cellstr(spm_select('FPList',GLMAnalDir,  ['MotionPhysio_Sess_' num2str(i) '_' SubjID '.mat']));
    assertSizeAtLeast(1,MotionPhysioRegFile, ['Motion Physio Regressor file not found for subject_' SubjID]);
    
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond = struct([]);
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi = RegFile;
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).regress = struct([]);
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi_reg = MotionPhysioRegFile;
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).hpf = 128;
end

save(fullfile(GLMAnalDir,['GLMSpecify_' SubjID '.mat']), 'matlabbatch');
spm_jobman('initcfg')
spm_jobman('run', matlabbatch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model Estimate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear matlabbatch

load([TempPD '/ModelEstimate.mat']);

SPMDir=GLMAnalDir;
SPMFile=cellstr(spm_select('FPList',SPMDir,  'SPM.mat'));
assertSizeIs(1, SPMFile, 'SPM file not found');
matlabbatch{1}.spm.stats.fmri_est.spmmat=SPMFile;

spm_jobman('initcfg')
spm_jobman('run', matlabbatch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contrast Manager
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear matlabbatch

load([TempPD '/ContrastManager.mat']);

% Load regressor names because there are so many (will use this to specify contrasts)
load(fullfile(GLMAnalDir,  ['Regressors_Sess_' num2str(1) '_' SubjID '.mat']));

% Use rows of identity matrix to make contrast of each regressor (versus implicit baseline)
contrastsAll = eye(length(names));

% Initialise separate T contrast counter
Tcount = 0;

% T contrasts
for regIndex = 1:length(names)
    Tcount = Tcount+1;
    Tcons(Tcount).c = contrastsAll(regIndex,:);
    Tcons(Tcount).name = names{regIndex};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SPMDir=GLMAnalDir;
SPMFile=cellstr(spm_select('FPList',SPMDir,  'SPM.mat'));
assertSizeIs(1, SPMFile, 'SPM file not found');

matlabbatch{1}.spm.stats.con.spmmat=SPMFile;
matlabbatch{1}.spm.stats.con.delete=1;

for k=1:length(Tcons)
    matlabbatch{1}.spm.stats.con.consess{k}.tcon.name=Tcons(k).name;
    matlabbatch{1}.spm.stats.con.consess{k}.tcon.convec=Tcons(k).c;
    matlabbatch{1}.spm.stats.con.consess{k}.tcon.sessrep='replsc';
end

save(fullfile(GLMAnalDir,['ContrastManager_' SubjID '.mat']), 'matlabbatch');
spm_jobman('initcfg')
spm_jobman('run', matlabbatch);