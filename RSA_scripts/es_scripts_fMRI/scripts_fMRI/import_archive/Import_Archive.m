function Import_Archive(archive,pth)
% Import a list of archive files into NIfTI files
% FORMAT Import_Archive(ARCHIVE,PTH) imports archive files ARCHIVE (can be
% .tar, .tar.gz) in output directory PTH.
%__________________________________________________________________________
%
% This function uses MATLAB function <untar>, available from MATLAB 7.0.4.
% If not, it tries to use system commands (tar on Unix and provided tar.exe
% gzip.exe on Windows).
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Methods & Physics Groups
% $Id: Import_Archive.m 149 2008-10-06 10:59:00Z guillaume $

SVNid = '$Rev: 149 $';

spm('FnBanner',mfilename,SVNid);
Finter = spm('FigName','Import archive...'); 

%-Retrieve input arguments
%--------------------------------------------------------------------------
if nargin < 1
    [archive,sts] = spm_select([1 Inf],'^.*\.tar(\.gz)?$','Select archive files');
    if ~sts, return; end
end
archive = cellstr(archive);

if nargin < 2
    [pth,sts] = spm_select([1],'dir','Select output directory');
    if ~sts, return; end
end

spm('Pointer','Watch');

%-Loop over archives
%--------------------------------------------------------------------------
for i=1:numel(archive)
    
    %-Loop overt archives
    %----------------------------------------------------------------------
    [p,n,e] = fileparts(archive{i});
    if strcmpi(e,'.gz')
        [p,n] = fileparts(fullfile(p,n));
    end
    fprintf('\r%-40s: %30s',sprintf('Archive %s',n),' ');
    
    %-Create output subdirectory
    %----------------------------------------------------------------------
    fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...create output directory');
    outdir = fullfile(pth,n);
    if ~exist(outdir,'dir')
        [sts,msg] = mkdir(pth,n);
        if ~sts, error(msg); end
    end
    
    %-Create temporary directory
    %----------------------------------------------------------------------
    fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...create temporary directory');
    %my_tempdir = tempname;
    [b, my_tempdir] = fileparts(tempname);
    my_tempdir = fullfile(outdir,['__SAFE_TO_DELETE__' my_tempdir]);
    if ~exist(my_tempdir,'dir')
        [sts,msg] = mkdir(my_tempdir);
        if ~sts, error(msg); end
    end
    
    %-Extract archive in a temporary directory
    %----------------------------------------------------------------------
    fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...extract archive');
    try
        farchive = feval(archive_extract_fcn,archive{i},my_tempdir);
    catch
        fprintf('\nDirectory "%s" should be manually deleted.\n',my_tempdir);
        % Could use onCleanup for Matlab 7.6
        %  cleanup = onCleanup(@()(deletetemporarydir(my_tempdir)));
        %  cleanup.delete;
        rethrow(lasterror);
    end
    
    %-Reconstruct images
    %----------------------------------------------------------------------
    fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...reconstruct images');
    f_ima = spm_select('FPList',my_tempdir,'^.*\.ima$');
    if ~(isempty(f_ima) || ... % work-around for stupid bug with FPList
            (size(f_ima,1) == 1 && size(f_ima,2) <= length(my_tempdir)+1))
        cwd = pwd;
        cd(outdir);
        
        fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...read DICOM headers');
        hdr = spm_dicom_headers(f_ima,true);
        
        fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...convert DICOM images');
        spm_dicom_convert(hdr,'all','flat','img');
        
        cd(cwd);
        ia_msg = 'done';
    else
        ia_msg = 'skipped';
    end
    
    %-Delete temporary files and directory
    %----------------------------------------------------------------------
    fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...delete temporary files');
    for j=1:numel(farchive)
        spm_unlink(farchive{j});
    end
    [sts, msg] = rmdir(my_tempdir); % could use the mode option
    if ~sts, warning(msg); end
    
    fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),['...' ia_msg]);
    
end

fprintf('\n');
spm('FigName','Import archive: done',Finter); spm('Pointer','Arrow')


%==========================================================================
function f = archive_extract_fcn

if usejava('jvm') 
    % Matlab untar function available from Matlab 7.0.4 (R14SP2)
    f = @untar;
else
    if ispc
        f = @windows_untar;
    elseif isunix || ismac
        f = @unix_untar;
    else
        error('Please use Matlab 7.0.4 or greater with Java VM.');
    end
end

%==========================================================================
function f = unix_untar(tarfilename,outputdir)

[p,n,e] = fileparts(tarfilename);
if strcmpi(e,'.gz')
    opt = 'xvfz';
else
    opt = 'xvf';
end

[s, w] = unix(['tar ' opt ' "' tarfilename '" -C "' outputdir '"']);
if s
    error('[unix_untar] Data extraction failed. Please check your input files.');
end

f = strread(w,'%s','delimiter',sprintf('\n'));
for i=1:length(f), f{i} = spm_select('CPath',f{i},outputdir); end

%==========================================================================
function f = windows_untar(tarfilename,outputdir)

tar_exe  = fullfile(fileparts(which(mfilename)),'tar.exe');
gzip_exe = fullfile(fileparts(which(mfilename)),'gzip.exe');

[p,n,e] = fileparts(tarfilename);
if strcmpi(e,'.gz')
    [s, w] = system(['"' gzip_exe '" -d "' tarfilename '"']);
    if s
        error('[win_gzip] Data extraction failed. Please check your input files.');
    end
    tarfilename = fullfile(p,n);
end

cwd = pwd;
cd(outputdir); % -C directive does not seem to be working
[s, w] = system(['"' tar_exe '" xvf "' tarfilename '"']);
cd(cwd);
if s
    error('[win_tar] Data extraction failed. Please check your input files.');
end

f = strread(w,'%s','delimiter',sprintf('\r\n'));
for i=1:length(f), f{i} = spm_select('CPath',f{i},outputdir); end
