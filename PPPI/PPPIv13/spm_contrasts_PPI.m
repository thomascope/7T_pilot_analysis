function SPM=spm_contrasts_PPI(SPM,Subject,SPMContrasts,Weighted,Method,GroupDir)
%Estimates PPI contrasts
%   SPM is the SPM file for PPI or the SPM structure for PPI after
%   estimation
%   Contrasts is a structure of the contrasts to be applied to the PPI
%   estimation, if not specified, computes all possible combinations
%   Weighted is numeric cutoff for the minimum time of an event/epoch/block
%   to use run averaging as opposed to trial averaging. Default is 0, so
%   all events will be analyzed using run averaging.
%   Method is a string that is either 'trad' or 'cond'.
%
%
% License:
%   Copyright (c) 2011-16, Donald G. McLaren and Aaron Schultz
%   All rights reserved.
%
%    Redistribution, with or without modification, is permitted provided that the following conditions are met:
%    1. Redistributions must reproduce the above copyright
%        notice, this list of conditions and the following disclaimer in the
%        documentation and/or other materials provided with the distribution.
%    2. All advertising materials mentioning features or use of this software must display the following acknowledgement:
%        This product includes software developed by the Harvard Aging Brain Project.
%    3. Neither the Harvard Aging Brain Project nor the
%        names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
%    4. You are not permitted under this Licence to use these files
%        commercially. Use for which any financial return is received shall be defined as commercial use, and includes (1) integration of all 	
%        or part of the source code or the Software into a product for sale or license by or on behalf of Licensee to third parties or (2) use 	
%        of the Software or any derivative of it for research with the final aim of developing software products for sale or license to a third 	
%        party or (3) use of the Software or any derivative of it for research with the final aim of developing non-software products for sale 
%        or license to a third party, or (4) use of the Software to provide any service to an external organisation for which payment is received.
%
%   THIS SOFTWARE IS PROVIDED BY DONALD G. MCLAREN (mclaren@nmr.mgh.harvard.edu) AND AARON SCHULTZ (aschultz@nmr.mgh.harvard.edu)
%   ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
%   FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
%   SPECIAL, EXEMPLARY,/usr/pubsw/common/scripts/fmri/Utilities_DGM/PPPI OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
%   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR 
%   TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
%   USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
%   DAMAGE.
%
%   Last modified on 2/14/2012 by Donald G. McLaren
%   (mclaren@nmr.mgh.harvard.edu)
%   GRECC, Bedford VAMC
%   Department of Neurology, Massachusetts General Hospital and Havard
%   Medical School
%
%   2/14/2012:
%   Updated the prefix portion to work with using multiple prefixes and
%   prefix structures. This will now allow the user to specify individual
%   runs.
%
%   7/7/2015:
%   Removed move commands
%% Program begins here
%Check input arguments
if ~isstruct(SPM)
    if exist(SPM,'file')==2
        load(SPM)
    else
        disp('PPI not estimated')
        return;
    end
end

cd(SPM.swd)
if nargin==1
    error('Subject must be specified')
end
if nargin==2 || isempty(SPMContrasts)
    try
        SPMContrasts=defContrasts(SPM,Weighted);
    catch
        SPMContrasts=defContrasts(SPM,0);
    end
end
if nargin<4
    Weighted=0;
end
if ~isempty(SPMContrasts) && iscellstr(SPMContrasts)
    try
        SPMContrasts=defContrasts(SPM,Weighted,0,SPMContrasts);
    catch
        SPMContrasts=defContrasts(SPM,0,0,SPMContrasts);
    end
end
if ~isfield(SPMContrasts,'Weighted')
    SPMContrasts(1).Weighted=Weighted;
end

%Configure Contrasts
ind=zeros(length(SPMContrasts),1);
display('Generate Contrast Vectors')
for ii = 1:length(SPMContrasts)
    if ~isfield(SPMContrasts(ii),'c') || isempty(SPMContrasts(ii).c)
        if isempty(SPMContrasts(ii).Weighted)
            SPMContrasts(ii).Weighted=Weighted;
        end
        if strcmpi(Method,'trad')
            try
                SPMContrasts(ii).c=createVec(SPMContrasts(ii).left,SPMContrasts(ii).right,SPM,SPMContrasts(ii).Weighted,SPMContrasts(ii).Prefix)';
            catch
                SPMContrasts(ii).c=createVec(SPMContrasts(ii).left,SPMContrasts(ii).right,SPM,SPMContrasts(ii).Weighted)';
            end
        else
            if isfield(SPMContrasts(ii),'Prefix') && ~isempty(SPMContrasts(ii).Prefix) && isfield(SPMContrasts(ii),'Contrail') && ~isempty(SPMContrasts(ii).Contrail)
                Prefix=checkprefix(SPMContrasts(ii).Prefix);
                try
                    SPMContrasts(ii).c=createVec(SPMContrasts(ii).left,SPMContrasts(ii).right,SPM,SPMContrasts(ii).Weighted,Prefix,SPMContrasts(ii).Contrail,SPMContrasts(ii).MinEvents,SPMContrasts(ii).MinEventsPer)';
                catch
                    SPMContrasts(ii).c=createVec(SPMContrasts(ii).left,SPMContrasts(ii).right,SPM,SPMContrasts(ii).Weighted,Prefix,SPMContrasts(ii).Contrail,SPMContrasts(ii).MinEvents)';
                end
                  
            elseif isfield(SPMContrasts(ii),'Prefix') && ~isempty(SPMContrasts(ii).Prefix)
                Prefix=checkprefix(SPMContrasts(ii).Prefix);
                try
                    SPMContrasts(ii).c=createVec(SPMContrasts(ii).left,SPMContrasts(ii).right,SPM,SPMContrasts(ii).Weighted,Prefix,[],SPMContrasts(ii).MinEvents,SPMContrasts(ii).MinEventsPer)';
                catch
                    SPMContrasts(ii).c=createVec(SPMContrasts(ii).left,SPMContrasts(ii).right,SPM,SPMContrasts(ii).Weighted,Prefix,[],SPMContrasts(ii).MinEvents)';
                end
            elseif isfield(SPMContrasts(ii),'Contrail') && ~isempty(SPMContrasts(ii).Contrail)
                try
                    SPMContrasts(ii).c=createVec(SPMContrasts(ii).left,SPMContrasts(ii).right,SPM,SPMContrasts(ii).Weighted,'PPI_',SPMContrasts(ii).Contrail,SPMContrasts(ii).MinEvents,SPMContrasts(ii).MinEventsPer)';
                catch
                    SPMContrasts(ii).c=createVec(SPMContrasts(ii).left,SPMContrasts(ii).right,SPM,SPMContrasts(ii).Weighted,'PPI_',SPMContrasts(ii).Contrail,SPMContrasts(ii).MinEvents)';
                end
            else
                try
                    SPMContrasts(ii).c=createVec(SPMContrasts(ii).left,SPMContrasts(ii).right,SPM,SPMContrasts(ii).Weighted,'PPI_',[],SPMContrasts(ii).MinEvents,SPMContrasts(ii).MinEventsPer)';
                catch
                    SPMContrasts(ii).c=createVec(SPMContrasts(ii).left,SPMContrasts(ii).right,SPM,SPMContrasts(ii).Weighted,'PPI_',[],SPMContrasts(ii).MinEvents)';
                end
            end
        end
        
        if mean(SPMContrasts(ii).c==0)~=1; ind(ii)=1; end
        if isempty(SPMContrasts(ii).name)
            if isfield(SPMContrasts(ii),'Prefix') && ~isempty(SPMContrasts(ii).Prefix) && isfield(SPMContrasts(ii),'Contrail') && ~isempty(SPMContrasts(ii).Contrail)
                tmp=[SPMContrasts(ii).Prefix '_PPI_' SPMContrasts(ii).left 'minus' SPMContrasts(ii).right '_' SPMContrasts(ii).Contrail];
            elseif isfield(SPMContrasts(ii),'Prefix') && ~isempty(SPMContrasts(ii).Prefix)
                tmp=[SPMContrasts(ii).Prefix '_PPI_' SPMContrasts(ii).left 'minus' SPMContrasts(ii).right];
            elseif isfield(SPMContrasts(ii),'Contrail') && ~isempty(SPMContrasts(ii).Contrail)
                tmp=['PPI_' SPMContrasts(ii).left 'minus' SPMContrasts(ii).right '_' SPMContrasts(ii).Contrail];
            else
                tmp=['PPI_' SPMContrasts(ii).left 'minus' SPMContrasts(ii).right];
            end
            tmp=sprintf('%s_',tmp{:});
            tmp=tmp(1:end-1);
            SPMContrasts(ii).name=tmp;
        else
            if iscellstr(SPMContrasts(ii).name)
                SPMContrasts(ii).name=['PPI_' cell2mat(SPMContrasts(ii).name)];
            else
                SPMContrasts(ii).name=['PPI_' SPMContrasts(ii).name];
            end
        end
        if isempty(SPMContrasts(ii).STAT)
            SPMContrasts(ii).STAT='T';
        end
    end
end
SPMContrasts=SPMContrasts(ind==1);
for ii = 1:length(SPMContrasts)
    disp('Generate xCon')
    xCon(ii) = spm_FcUtil('Set',SPMContrasts(ii).name,SPMContrasts(ii).STAT,'c',SPMContrasts(ii).c,SPM.xX.xKXs);
end

%Compute Contrasts
disp('Generate Contrasts')
try
    init=length(SPM.xCon);
catch
    init=0;
end
if init~=0
    SPM.xCon(init+1:init+length(xCon)) = xCon;
else
    SPM.xCon = xCon;
end
SPM = spm_contrasts(SPM,init+1:length(SPM.xCon));
% Move contrasts
for ii=(1+init):numel(SPM.xCon)
    disp('Moving Contrast Images');
    DD=['_' date];
    f1 = SPM.xCon(ii).Vcon.fname;
    if isempty(strfind(f1,'.nii'))
        f2 = [SPM.xCon(ii).Vcon.fname(1:end-3) 'hdr'];
        ext='.img';
    else
        ext='.nii';
    end
    [junk fname]=fileparts(SPM.xCon(ii).Vcon.fname); fname=fname(6:end); clear junk
    if length(SPM.xCon(ii).name)==4 && strcmp(fname,SPM.xCon(ii).name);
        disp(['Contrast: ' SPM.xCon(ii).Vcon.fname ' was not moved.'])
        continue
    end
    err=0;
    if ~exist([f1(1:4) SPM.xCon(ii).name '_' Subject ext],'file')
        try
            movefile(f1, [f1(1:4) SPM.xCon(ii).name '_' Subject ext],'f');
            try movefile(f2, [f2(1:4) SPM.xCon(ii).name '_' Subject '.hdr'],'f'); end
        catch
            disp(['error moving contrast '  SPM.xCon(ii).name])
            err=1;
        end
        if exist('GroupDir','var') && exist(GroupDir,'dir')==7 && err==0
            copyfile([f1(1:4) SPM.xCon(ii).name '_' Subject ext],GroupDir,'f');
            try copyfile([f2(1:4) SPM.xCon(ii).name '_' Subject '.hdr'],GroupDir,'f'); end
        end
    else
        try
            movefile([f1(1:4) SPM.xCon(ii).name '_' Subject ext], [f1(1:4) SPM.xCon(ii).name DD '_' Subject  ext],'f');
            try movefile([f2(1:4) SPM.xCon(ii).name '_' Subject '.hdr'], [f2(1:4) SPM.xCon(ii).name DD '_' Subject  '.hdr'],'f'); end
            movefile(f1, [f1(1:4) SPM.xCon(ii).name '_' Subject ext],'f');
            try movefile(f2, [f2(1:4) SPM.xCon(ii).name '_' Subject '.hdr'],'f'); end
            for kk = 1:numel(SPM.xCon)
                if strcmp(SPM.xCon(kk).Vcon.fname,[f1(1:4) SPM.xCon(ii).name '_' Subject ext]) > 0
                    SPM.xCon(kk).Vcon.fname=[f1(1:4) SPM.xCon(ii).name DD '_' Subject ext];
                    SPM.xCon(kk).name=[SPM.xCon(ii).name DD];
                    break
                end
            end
           
        catch
            disp(['error moving contrast '  SPM.xCon(ii).name])
            err=1;
        end
        if exist('GroupDir','var') && exist(GroupDir,'dir')==7 && err==0
            copyfile([f1(1:4) SPM.xCon(ii).name '_' Subject ext],GroupDir,'f');
            try copyfile([f2(1:4) SPM.xCon(ii).name '_' Subject '.hdr'],GroupDir,'f'); end
        end
    end
    if ~strcmp(SPM.xCon(ii).Vcon.fname,[f1(1:4) SPM.xCon(ii).name '_' Subject ext])
        SPM.xCon(ii).Vcon.fname=[f1(1:4) SPM.xCon(ii).name '_' Subject ext];
    end
    f1 = SPM.xCon(ii).Vspm.fname;
    if isempty(strfind(f1,'.nii'))
        f2 = [SPM.xCon(ii).Vspm.fname(1:end-3) 'hdr'];
        ext='.img';
    else
        ext='.nii';
    end
    err=0;
    if ~exist([f1(1:4) '_' SPM.xCon(ii).name '_' Subject ext],'file')
       try
            movefile(f1, [f1(1:4) '_' SPM.xCon(ii).name '_' Subject ext],'f');
            try movefile(f2, [f2(1:4) '_' SPM.xCon(ii).name '_' Subject '.hdr'],'f'); end
       catch
           disp(['error moving contrast '  SPM.xCon(ii).name])
           err=1;
       end
       if exist('GroupDir','var') && exist(GroupDir,'dir')==7 && err==0
            copyfile([f1(1:4) '_' SPM.xCon(ii).name '_' Subject ext],GroupDir,'f');
            try copyfile([f2(1:4) '_' SPM.xCon(ii).name '_' Subject '.hdr'],GroupDir,'f'); end %#ok<*TRYNC>
        end
    else
        try
            movefile([f1(1:4) '_' SPM.xCon(ii).name '_' Subject ext], [f1(1:4) '_' SPM.xCon(ii).name DD '_' Subject ext],'f');
            try movefile([f2(1:4) '_' SPM.xCon(ii).name '_' Subject '.hdr'], [f2(1:4) '_' SPM.xCon(ii).name DD '_' Subject '.hdr'],'f'); end
            movefile(f1, [f1(1:4) '_' SPM.xCon(ii).name '_' Subject ext],'f');
            try movefile(f2, [f2(1:4) '_' SPM.xCon(ii).name '_' Subject '.hdr'],'f'); end
            for kk = 1:numel(SPM.xCon)
                if strcmp(SPM.xCon(kk).Vspm.fname,[f1(1:4) SPM.xCon(ii).name '_' Subject ext]) > 0
                    SPM.xCon(kk).Vspm.fname=[f1(1:4) SPM.xCon(ii).name DD '_' Subject ext];
                    break
                end
            end
       catch
           disp(['error moving contrast '  SPM.xCon(ii).name])
           err=1;
        end
        if exist('GroupDir','var') && exist(GroupDir,'dir')==7 && err==0
                copyfile([f1(1:4) '_' SPM.xCon(ii).name '_' Subject ext],GroupDir,'f');
                try copyfile([f2(1:4) '_' SPM.xCon(ii).name '_' Subject '.hdr'],GroupDir,'f'); end
        end
    end 
    if ~strcmp(SPM.xCon(ii).Vspm.fname,[f1(1:4) '_' SPM.xCon(ii).name '_' Subject ext])
        SPM.xCon(ii).Vspm.fname=[f1(1:4) '_' SPM.xCon(ii).name '_' Subject ext];
    end 
end
save SPM.mat SPM
    
function Prefix=checkprefix(Prefixin)
    if iscellstr(Prefixin)
        Prefix=cell(size(Prefixin));
        for ii=1:numel(Prefixin)
            Prefix{ii}=[Prefixin{ii} '.*PPI_'];
        end
    elseif isstruct(Prefixin)
        Prefix.Left=cell(size(Prefixin.Left));
        for ii=1:numel(Prefixin.Left)
            try
                Prefix.Left{ii}=[Prefixin.Left{ii} '.*PPI_'];
            catch
                Prefix.Left{ii}=[Prefixin.Left '.*PPI_'];
            end
        end
        Prefix.Right=cell(size(Prefixin.Right));
        for ii=1:numel(Prefixin.Right)
            try
                Prefix.Right{ii}=[Prefixin.Right{ii} '.*PPI_'];
            catch
                Prefix.Right{ii}=[Prefixin.Right '.*PPI_'];
            end
        end
    elseif ischar(Prefixin)
        Prefix=[Prefixin '.*PPI_'];
    else
        error('Prefix is not valid. Must be a cellstr or structure')
    end
return

    