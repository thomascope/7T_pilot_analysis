function dartel_template_to_mni()
%Use the Template_6_2mni.mat file generated from the 'Normalise to MNI
%space routine to put Template_6 in MNI space

gmmult=5000; %4000 or 5000
wmmult=7500; %6500 or 7500


%Any representitive warped image that you want the template to match. It
%will have the same bounding box.
%target = '/work/imagingC/Rowe/p00259/longitudinal_vbm/subs/pt001_22033/mprage_125_sess_01/s_0_wc1avg_pt001_22033_mprage_125_sess_01.nii';
target = '/imaging/tc02/vespa/scans/PNFA_VBM/tom/Stats/copy_of_full_vbm_for_Bayes_null_H0_testing/FWE_05_GM_PNFA.nii';

%load('/work/imagingC/Rowe/p00259/longitudinal_vbm/subs/pt001_22033/mprage_125_sess_01/Template_psp10_6_2mni.mat');
load('/imaging/tc02/vespa/scans/PNFA_VBM/tom/p00259/ppa_20218/mprage_125/Template_6_2mni.mat');

M = mni.affine;

P='/imaging/tc02/vespa/scans/PNFA_VBM/tom/p00259/ppa_20218/mprage_125/Template_6_2mni_mat.nii';

Nii = nifti(P);
Nii.mat = M;
Nii.mat_intent = 4;

create(Nii);

prefix = 'r';

flags.mask = 0;
flags.mean = 0;
flags.interp = 4;
flags.which = 1; % reslice images 2:n
flags.wrap = [ 0 0 0];
flags.prefix = prefix;  

source=P;

Ps=deblank(char(target, source));

spm_reslice(Ps,flags);

%Create an average T1.
[pth,nam,ext,num] = spm_fileparts(P);
% Pi=fullfile(pth,strcat(prefix,nam,ext));
gm=fullfile(pth,strcat(prefix,nam,ext,',1'));
wm=fullfile(pth,strcat(prefix,nam,ext,',2'));

clear flags

flags.dmtx = 0;
flags.mask = 0;
flags.interp = 0;
flags.dtype =  spm_type('int16'); % 

% f=sprintf('i1*%d+i2*%d',gmmult,wmmult);
f=sprintf('i1');
Po=fullfile(pth,strcat(prefix,nam,'_GM',ext));

Vo = spm_imcalc(char(gm), Po, f ,flags);
