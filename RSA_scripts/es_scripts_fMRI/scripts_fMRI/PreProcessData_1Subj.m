function [] =PreProcessData_1Subj(SubjDir,FuncDir,StrDir,TempPD,SpmPD,ScriptsPD)

NumSess = length(FuncDir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load SPM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(ScriptsPD); % So that any customised code in project scripts folder is prioritised (like spm_my_defaults)
addpath(SpmPD); % So that correct version of SPM is loaded
spm('FMRI');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create fieldmap vdm files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load([TempPD '/FieldMaps.mat']);

% default values
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.et =[10 12.46];
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.maskbrain=1;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.blipdir =-1;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.tert =37;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.epifm=0;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsval.ajm=0;
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.template = {fullfile(SpmPD,'FieldMap/T1.nii')};
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.matchvdm=1;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.writeunwarped=0;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.matchanat=0;

PhaseImg = cellstr(spm_select('FPList', fullfile(SubjDir,'FieldMaps','Phase'), '.nii'));
MagnImg = cellstr(spm_select('FPList', fullfile(SubjDir,'FieldMaps','Magnitude'), '.nii'));

% Select just the first of two mag images
MagnImg = MagnImg(1);

assertSizeIs(1, PhaseImg, 'Expected one phase map');
assertSizeIs(1, MagnImg, 'Expected one magnitude map');

matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.phase =PhaseImg;
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.magnitude =MagnImg;

% EPI for each session
matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.session = struct();
for i = 1:NumSess
    % Identify first EPI of session
    SessEpiDir = FuncDir{i};
    FirstEpi = cellstr(spm_select('FPList', SessEpiDir, '.*\.nii$'));
    FirstEpi = FirstEpi(1);
    assertSizeIs(1, FirstEpi, ['Fieldmap: expected one EPI for session' num2str(i)]);   
    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.session(i).epi = FirstEpi;
end
save(fullfile(SubjDir,'FieldMapSpecify.mat'), 'matlabbatch');
spm_jobman('initcfg')
spm_jobman('run', matlabbatch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Realign and unwarp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear matlabbatch

load([TempPD '/RealignAndUnwarp.mat']);
matlabbatch{1}.spm.spatial.realignunwarp.data = struct();

for i = 1:NumSess
    % EPIs for session
    SessEpiDir = FuncDir{i};
    EpiFiles = cellstr(spm_select('FPList', SessEpiDir, '^f.*\.nii$'));
    assertSizeAtLeast(1, EpiFiles, ['Realign: expected EPIs for session' num2str(i)]);
    matlabbatch{1}.spm.spatial.realignunwarp.data(i).scans = EpiFiles;
    clear SessEpiDir
    
    % Phase map for session
    PhaseDir = fullfile(SubjDir,'FieldMaps','Phase');
    VdmFile = cellstr(spm_select('FPList', PhaseDir, ['.*session' num2str(i) '\.nii$']));
    assertSizeIs(1, VdmFile, ['Realign: expected phasemap for session' num2str(i)]);
    matlabbatch{1}.spm.spatial.realignunwarp.data(i).pmscan = VdmFile;
end

save(fullfile(SubjDir,'RealignAndUnwarp.mat'), 'matlabbatch');
spm_jobman('initcfg')
spm_jobman('run', matlabbatch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coregistration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear matlabbatch

load([TempPD '/Coregister.mat']);

StrFile = cellstr(spm_select('FPList', StrDir, '^s.*\.nii$'));
EpiDir = FuncDir{1};
MeanFunctional = cellstr(spm_select('FPList',EpiDir, '^mean.*\.nii$'));
assertSizeIs(1, MeanFunctional, 'Coregistration: expected mean functional');

matlabbatch{1}.spm.spatial.coreg.estimate.source = StrFile;
matlabbatch{1}.spm.spatial.coreg.estimate.ref = MeanFunctional;

% % Specify 'Other images' to coregister (in this case, all the individual functionals). Helps when default Coreg fails?
% IndividualFunctional = {};
% for i=1:NumSess
%     EpiDir = FuncDir{i};
%     IndividualFunctional = [IndividualFunctional; cellstr(spm_select('FPList',EpiDir, '^f.*\.nii$'))];
% end
% matlabbatch{1}.spm.spatial.coreg.estimate.source = MeanFunctional;
% matlabbatch{1}.spm.spatial.coreg.estimate.other = IndividualFunctional;
% matlabbatch{1}.spm.spatial.coreg.estimate.ref = StrFile;

save(fullfile(SubjDir,'Coregister.mat'), 'matlabbatch');
spm_jobman('initcfg')
spm_jobman('run', matlabbatch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear matlabbatch

load([TempPD '/Segment.mat']);

StrFile = cellstr(spm_select('FPList', StrDir, '^s.*\.nii$'));
assertSizeIs(1, StrFile, 'Segmentation: expected structural');
matlabbatch{1}.spm.spatial.preproc.channel.vols = StrFile;
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 1;
matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];
matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {fullfile(SpmPD,'/tpm/TPM.nii,1')};
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {fullfile(SpmPD,'/tpm/TPM.nii,2')};
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {fullfile(SpmPD,'/tpm/TPM.nii,3')};
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {fullfile(SpmPD,'/tpm/TPM.nii,4')};
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {fullfile(SpmPD,'/tpm/TPM.nii,5')};
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {fullfile(SpmPD,'/tpm/TPM.nii,6')};

save(fullfile(SubjDir,'Segment.mat'), 'matlabbatch');
spm_jobman('initcfg')
spm_jobman('run', matlabbatch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalize Functionals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear matlabbatch

load([TempPD '/Normalize.mat']);

% Gather functionals from all sessions
AllSessFunctionals = {};
for i = 1:NumSess
    SessEpiDir = FuncDir{i};
    SessFunctionals = cellstr(spm_select('FPList', SessEpiDir, '^uf.*\.nii$')); 
    assertSizeAtLeast(1, SessFunctionals, ['Normalize functionals: expected EPIs for session' num2str(i)]);
    AllSessFunctionals = [AllSessFunctionals ; SessFunctionals];
end

% Deformation field
DefFile = cellstr(spm_select('FPList', StrDir, '^y.*\.nii$'));
assertSizeIs(1, DefFile, 'Normalize functionals: expected seg file');

matlabbatch{1}.spm.spatial.normalise.write.subj.resample = AllSessFunctionals;
matlabbatch{1}.spm.spatial.normalise.write.subj.def = DefFile;

save(fullfile(SubjDir,'NormalizeFunctional.mat'), 'matlabbatch');
spm_jobman('initcfg')
spm_jobman('run', matlabbatch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalize Structural
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear matlabbatch

load([TempPD '/Normalize.mat']);

% Structural

StrFile = cellstr(spm_select('FPList',StrDir, '^ms.*\.nii$'));
assertSizeIs(1, StrFile, 'Normalize structural: expected structural');

% Seg file from segmentation step
DefFile = cellstr(spm_select('FPList', StrDir, '^y.*\.nii$'));
assertSizeIs(1, DefFile, 'Normalize functionals: expected seg file');

matlabbatch{1}.spm.spatial.normalise.write.subj.resample = StrFile;
matlabbatch{1}.spm.spatial.normalise.write.subj.def = DefFile;

save(fullfile(SubjDir,'NormalizeStructural.mat'), 'matlabbatch');
spm_jobman('initcfg')
spm_jobman('run', matlabbatch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Smooth Functionals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear matlabbatch SessFunctionals

load([TempPD '/Smooth.mat']);

% Gather functionals from all sessions
AllSessFunctionals = {};
for i = 1:NumSess
    SessEpiDir = FuncDir{i};
    SessFunctionals = cellstr(spm_select('FPList', SessEpiDir, '^wuf.*\.nii$'));
    assertSizeAtLeast(1, SessFunctionals, ['Smooth: expected EPIs for session' num2str(i)]);
    
    AllSessFunctionals = [AllSessFunctionals; SessFunctionals ];
end

matlabbatch{1}.spm.spatial.smooth.data = AllSessFunctionals;
matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];

save(fullfile(SubjDir,'Smooth.mat'), 'matlabbatch');
spm_jobman('initcfg')
spm_jobman('run', matlabbatch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rename graphics .ps output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find most recent SPM graphics output file
file = dir('spm*.ps');
% work out subject ID from SubjDir variable
startInd = strfind(SubjDir,'/');
SubjID = SubjDir(startInd(end)+1:end);
% prepend subject ID to graphics file name
movefile(file(1).name,[SubjID '_' file(1).name]);
