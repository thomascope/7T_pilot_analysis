function P = spm_mesh_project_tc(M, dat, sampling_distance, method, varargin)
% Project volumetric data onto a mesh
% FORMAT P = spm_mesh_project(M, dat, method)
% M        - a patch structure, a handle to a patch 
%            or a [nx3] vertices array
% dat      - a structure array [1xm] with fields dim, mat, XYZ and t 
%            (see spm_render.m)
%            or a structure array [1xm] with fields mat and dat
%            or a char array/cellstr of image filenames
% method   - interpolation method {'nn'}
% varargin - other parameters required by the interpolation method
%
% P        - a [mxn] curvature vector
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_project.m 3409 2009-09-18 15:40:25Z guillaume $

if ishandle(M)
    V = get(M,'Vertices');
elseif isstruct(M)
    V = M.vertices;
else
    V = M;
end

if nargin < 4, method = 'nn'; end
if ~strcmpi(method,'nn')
    error('Only Nearest Neighbours interpolation is available.');
end

if ischar(dat), dat = cellstr(dat); end
P = zeros(length(dat),size(V,1));
for i=1:numel(dat)
    if iscellstr(dat)
        v      = spm_vol(dat{i});
        Y      = spm_read_vols(v);
        mat    = v.mat;
    elseif isfield(dat,'dat')
        Y      = dat(i).dat;
        mat    = dat(i).mat;
    else
        Y      = zeros(dat(i).dim(1:3)');
        OFF    = dat(i).XYZ(1,:) + dat(i).dim(1)*(dat(i).XYZ(2,:)-1 + dat(i).dim(2)*(dat(i).XYZ(3,:)-1));
        Y(OFF) = dat(i).t; % .* (dat(i).t > 0);
        mat    = dat(i).mat;
    end
    XYZ        = double(inv(mat)*[V';ones(1,size(V,1))]);
    
    disp(['Sampling volume from mesh at a distance of ' num2str(sampling_distance) ' voxels. This can take some time'])

    all_data_locations = dat(i).XYZ(1:3,~isnan(dat(i).t)&dat(i).t~=0);
    for this_loc = 1:size(XYZ,2)
        if mod(this_loc,100) == 1
            disp(['Sampling vertex ' num2str(this_loc) ' of ' num2str(size(XYZ,2))])
        end

        these_distances = pdist2(XYZ(1:3,this_loc)',all_data_locations');
        these_points = these_distances<sampling_distance;
        if sum(these_points)~=0
            these_data_locations = all_data_locations(:,these_points);
            these_data = zeros(1,size(these_data_locations,2));
            for this_voxel = 1:size(these_data_locations,2)
                these_data(this_voxel) = Y(these_data_locations(1,this_voxel),these_data_locations(2,this_voxel),these_data_locations(3,this_voxel));
            end
            P(i,this_loc) = max(these_data);
        else
            P(i,this_loc) = NaN;
        end
    end

    
%     disp(['Sampling volume from mesh at a distance of ' num2str(sampling_distance) 'mm. This can take some time'])
%     
%     for this_loc = 1:size(XYZ,2)
%         if mod(this_loc,100) == 1
%             disp(['Sampling vertex ' num2str(this_loc) ' of ' num2str(size(XYZ,2))])
%         end
%         
%         
%         
%         these_distances = pdist2(XYZ(1:3,this_loc)',dat(i).XYZ(1:3,~isnan(dat(i).t))');
%         these_points = these_distances<sampling_distance;
%         if ~isempty(these_data)
%             these_data =  spm_sample_vol(Y,dat(i).XYZ(1,these_points),dat(i).XYZ(2,these_points),dat(i).XYZ(2,these_points),0);
%             P(i,this_loc) = max(these_data);
%         else
%             P(i,this_loc) = NaN;
%         end
%     end
    
    
end
