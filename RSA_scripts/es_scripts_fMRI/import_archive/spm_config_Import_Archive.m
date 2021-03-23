function opts = spm_config_Import_Archive
% SPM5 Configuration file for Import_Archive toolbox
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Methods & Physics Groups
% $Id: spm_config_Import_Archive.m 148 2008-10-06 10:37:41Z guillaume $

%--------------------------------------------------------------------------

data.type = 'files';
data.name = 'Archive files';
data.tag  = 'data';
data.filter = '^.*\.tar(\.gz)?$';
data.num  = [1 Inf];
data.help = {'Select archive files.'};

%--------------------------------------------------------------------------

outdir.type = 'files';
outdir.name = 'Output directory';
outdir.tag  = 'outdir';
outdir.filter = 'dir';
outdir.num  = 1;
outdir.help = {'Select a directory where files are written.'};

%--------------------------------------------------------------------------

opts.type = 'branch';
opts.name = 'Import Archive';
opts.tag  = 'imparch';
opts.val  = {data,outdir};
opts.prog = @imparch;
opts.help = {[...
'This routine import a list of archive files (.tar.gz) into NIfTI files.']};

%__________________________________________________________________________
function imparch(job)

addpath(fullfile(spm('Dir'),'toolbox','Import_Archive'),'-begin');

Import_Archive(strvcat(job.data),job.outdir{1});
