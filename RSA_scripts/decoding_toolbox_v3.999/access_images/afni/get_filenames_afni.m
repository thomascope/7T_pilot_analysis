% function varargout = get_filenames_afni(varargin)
% 
% Function to get file names using afni_matlab. Supports both wildcard and
% regular expressions for file names (not for paths).
%
% INPUT:
%       file_path (1 x n string)
%       optional: filenames with wildcards (e.g. *.BRIK or rf*.BRIK), or 
%                 filenames with regular expressions starting with 
%                 REGEXP: (e.g. 'REGEXP:^rf.*\.BRIK$')
%
% OUTPUT: filenames as n x m char array (n = number of files) with space padding at end
%
% Example calls:
%       fnames = get_filenames_afni('/home/resultsdir','*.BRIK') 
%           will select all nifti files in the provided folder
%
%       fnames = get_filenames_afni('/home/resultsdir/*.BRIK')
%           same result as above
%
%       fnames = get_filenames_afni('/home/resultsdir','REGEXP:^results.*\.BRIK$')
%           will select all nifti files in the provided folder starting
%           with results and ending with .BRIK

function varargout = get_filenames_afni(varargin)


if nargin == 1
    [fpath,fname,fext] = fileparts(varargin{1});
    varargin{1} = fpath;
    varargin{2} = [fname fext];
end

if length(varargin) >= 2
    if strfind(varargin{2},'REGEXP:') == 1
        fname_regexp = varargin{2}(8:end);
    else
        fname_regexp = wildcard2regexp(varargin{2});
    end
else
    fname_regexp = '.*';
end

h = dir(varargin{1});
% exclude directories
h([h.isdir]) = [];
n_files = length(h);
fnames = cell(n_files,1);
keepind = false(n_files,1);
for i_file = 1:n_files
    keepind(i_file) = ~isempty(regexp(h(i_file).name,fname_regexp, 'once'));
    if keepind(i_file)
        fnames{i_file} = fullfile(varargin{1},h(i_file).name);
    end
end

varargout = {char(fnames(keepind))};