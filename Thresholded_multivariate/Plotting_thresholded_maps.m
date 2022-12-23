%Display all thresholded maps of interest (i.e. in this folder), scaled to
%their data range.
clear cfg
% cfg.plots = [1:2];
cfg.plots = [1];
% cfg.symmetricity = 'symmetrical'; % Works well for jet, but not hot or cool
% cfg.normalise = 1;
% cfg.threshold = [5 40];
cfg.inflate = 25;
addpath(genpath(['../atlas_Neuromorphometrics/']))
cfg.normalise = 0;
% cfg.rendfile = './lh.inflated1.surf.gii';
cfg.rendfile = './BrainMesh_ICBM152.lh.gii'; % Different overlay brains

all_scans = dir('./univariate/*.nii');

for sampling_distance = [3];
    cfg.sampling_distance = sampling_distance; % tc modification - When mesh is lower resolution than overlay image this can be problematic, now displays largest value within a Euclidean distance specified in number of voxels. Can be slow for high resolution images
    for i = 1:length(all_scans)
        
%         this_scan_data = spm_read_vols(spm_vol(all_scans(i).name));
%         upper_range = max(max(max(this_scan_data)));
%         lower_range = min(min(min(abs(this_scan_data))));
upper_range = 12.72; %Univariate maximum
%lower_range = 3.3748; % p<0.001 unc
lower_range = 0; % data already thresholded
        
        % spm_smooth(all_scans(i).name,['s5x_' all_scans(i).name],[5 0 0])
        % this_smoothed_scan_data = spm_read_vols(spm_vol(['s5x_' all_scans(i).name]));
        % upper_range_smooth = max(max(max(this_scan_data)));
        % lower_range_smooth = min(min(min(abs(this_scan_data))));
        
        %cfg.threshold = [3.3748 6.5]; %p=0.001
        cfg.threshold = [lower_range upper_range]; %scaled to data
        % cfg.threshold = [lower_range_smooth upper_range_smooth]; %scaled to data
        
        jp_spm8_surfacerender2_version_tc(['./univariate/' all_scans(i).name],'hot',cfg)
        
        savepath = ['./rendered_images/' all_scans(i).name(1:end-4) '_' num2str(lower_range,'%.2f') '_' num2str(upper_range,'%.2f') '_' num2str(cfg.sampling_distance) 'voxelsampling'];
        savepath = strrep(savepath,'>','-');
        %eval(['export_fig ' savepath '.png -transparent -m2'])
        eval(['export_fig ' savepath '.png -transparent'])
        close all
    end
end

cfg.rendfile = './BrainMesh_ICBM152.gii'; % For right hemisphere - twice as many vertices so twoce as long to compute.
cfg.plots = [2];
for sampling_distance = [3];
    cfg.sampling_distance = sampling_distance; % tc modification - When mesh is lower resolution than overlay image this can be problematic, now displays largest value within a Euclidean distance specified in number of voxels. Can be slow for high resolution images
    for i = 1:length(all_scans)
        
%         this_scan_data = spm_read_vols(spm_vol(all_scans(i).name));
%         upper_range = max(max(max(this_scan_data)));
%         lower_range = min(min(min(abs(this_scan_data))));
upper_range = 12.72; %Univariate maximum
%lower_range = 3.3748; % p<0.001 unc
lower_range = 0; % data already thresholded
        
        % spm_smooth(all_scans(i).name,['s5x_' all_scans(i).name],[5 0 0])
        % this_smoothed_scan_data = spm_read_vols(spm_vol(['s5x_' all_scans(i).name]));
        % upper_range_smooth = max(max(max(this_scan_data)));
        % lower_range_smooth = min(min(min(abs(this_scan_data))));
        
        %cfg.threshold = [3.3748 6.5]; %p=0.001
        cfg.threshold = [lower_range upper_range]; %scaled to data
        % cfg.threshold = [lower_range_smooth upper_range_smooth]; %scaled to data
        
        jp_spm8_surfacerender2_version_tc(['./univariate/' all_scans(i).name],'hot',cfg)
        
        savepath = ['./rendered_images/' all_scans(i).name(1:end-4) '_' num2str(lower_range,'%.2f') '_' num2str(upper_range,'%.2f') '_' num2str(cfg.sampling_distance) 'voxelsampling_rh'];
        savepath = strrep(savepath,'>','-');
        %eval(['export_fig ' savepath '.png -transparent -m2'])
        eval(['export_fig ' savepath '.png -transparent'])
        close all
    end
end

cfg.rendfile = './BrainMesh_ICBM152.lh.gii'; % Back to lh only
cfg.plots = [1];
all_scans = dir('./multivariate_segments/*.nii');

for sampling_distance = [3];
    cfg.sampling_distance = sampling_distance; % tc modification - When mesh is lower resolution than overlay image this can be problematic, now displays largest value within a Euclidean distance specified in number of voxels. Can be slow for high resolution images
    for i = 1:length(all_scans)
        
        this_scan_data = spm_read_vols(spm_vol(['./multivariate_segments/' all_scans(i).name]));
         upper_range = max(max(max(this_scan_data)));
%         lower_range = min(min(min(abs(this_scan_data))));
upper_range = 9.84 % For matching shared segment scales
lower_range = 0; % data already thresholded
        
        % spm_smooth(all_scans(i).name,['s5x_' all_scans(i).name],[5 0 0])
        % this_smoothed_scan_data = spm_read_vols(spm_vol(['s5x_' all_scans(i).name]));
        % upper_range_smooth = max(max(max(this_scan_data)));
        % lower_range_smooth = min(min(min(abs(this_scan_data))));
        
        %cfg.threshold = [3.3748 6.5]; %p=0.001
        cfg.threshold = [lower_range upper_range]; %scaled to data
        % cfg.threshold = [lower_range_smooth upper_range_smooth]; %scaled to data
        
        jp_spm8_surfacerender2_version_tc(['./multivariate_segments/' all_scans(i).name],'hot',cfg)
        
        savepath = ['./rendered_images/' all_scans(i).name(1:end-4) '_' num2str(lower_range,'%.2f') '_' num2str(upper_range,'%.2f') '_' num2str(cfg.sampling_distance) 'voxelsampling'];
        savepath = strrep(savepath,'>','-');
        %eval(['export_fig ' savepath '.png -transparent -m2'])
        eval(['export_fig ' savepath '.png -transparent'])
        close all
    end
end

all_scans = dir('./multivariate_predictions/*.nii');

for sampling_distance = [3];
    cfg.sampling_distance = sampling_distance; % tc modification - When mesh is lower resolution than overlay image this can be problematic, now displays largest value within a Euclidean distance specified in number of voxels. Can be slow for high resolution images
    for i = 1:length(all_scans)
        
        this_scan_data = spm_read_vols(spm_vol(['./multivariate_predictions/' all_scans(i).name]));
         upper_range = max(max(max(this_scan_data)));
%         lower_range = min(min(min(abs(this_scan_data))));
upper_range = 5.97; % For matching prediction scales
lower_range = 0; % data already thresholded
        
        % spm_smooth(all_scans(i).name,['s5x_' all_scans(i).name],[5 0 0])
        % this_smoothed_scan_data = spm_read_vols(spm_vol(['s5x_' all_scans(i).name]));
        % upper_range_smooth = max(max(max(this_scan_data)));
        % lower_range_smooth = min(min(min(abs(this_scan_data))));
        
        %cfg.threshold = [3.3748 6.5]; %p=0.001
        cfg.threshold = [lower_range upper_range]; %scaled to data
        % cfg.threshold = [lower_range_smooth upper_range_smooth]; %scaled to data
        
        jp_spm8_surfacerender2_version_tc(['./multivariate_predictions/' all_scans(i).name],'hot',cfg)
        
        savepath = ['./rendered_images/' all_scans(i).name(1:end-4) '_' num2str(lower_range,'%.2f') '_' num2str(upper_range,'%.2f') '_' num2str(cfg.sampling_distance) 'voxelsampling'];
        savepath = strrep(savepath,'>','-');
        %eval(['export_fig ' savepath '.png -transparent -m2'])
        eval(['export_fig ' savepath '.png -transparent'])
        close all
    end
end

cfg.rendfile = './BrainMesh_ICBM152.gii'; % Now rh
cfg.plots = [2];
all_scans = dir('./multivariate_segments/*.nii');

for sampling_distance = [3];
    cfg.sampling_distance = sampling_distance; % tc modification - When mesh is lower resolution than overlay image this can be problematic, now displays largest value within a Euclidean distance specified in number of voxels. Can be slow for high resolution images
    for i = 1:length(all_scans)
        
        this_scan_data = spm_read_vols(spm_vol(['./multivariate_segments/' all_scans(i).name]));
         upper_range = max(max(max(this_scan_data)));
%         lower_range = min(min(min(abs(this_scan_data))));
upper_range = 9.84 % For matching shared segment scales
lower_range = 0; % data already thresholded
        
        % spm_smooth(all_scans(i).name,['s5x_' all_scans(i).name],[5 0 0])
        % this_smoothed_scan_data = spm_read_vols(spm_vol(['s5x_' all_scans(i).name]));
        % upper_range_smooth = max(max(max(this_scan_data)));
        % lower_range_smooth = min(min(min(abs(this_scan_data))));
        
        %cfg.threshold = [3.3748 6.5]; %p=0.001
        cfg.threshold = [lower_range upper_range]; %scaled to data
        % cfg.threshold = [lower_range_smooth upper_range_smooth]; %scaled to data
        
        jp_spm8_surfacerender2_version_tc(['./multivariate_segments/' all_scans(i).name],'hot',cfg)
        
        savepath = ['./rendered_images/' all_scans(i).name(1:end-4) '_' num2str(lower_range,'%.2f') '_' num2str(upper_range,'%.2f') '_' num2str(cfg.sampling_distance) 'voxelsampling_rh'];
        savepath = strrep(savepath,'>','-');
        %eval(['export_fig ' savepath '.png -transparent -m2'])
        eval(['export_fig ' savepath '.png -transparent'])
        close all
    end
end

all_scans = dir('./multivariate_predictions/*.nii');

for sampling_distance = [3];
    cfg.sampling_distance = sampling_distance; % tc modification - When mesh is lower resolution than overlay image this can be problematic, now displays largest value within a Euclidean distance specified in number of voxels. Can be slow for high resolution images
    for i = 1:length(all_scans)
        
        this_scan_data = spm_read_vols(spm_vol(['./multivariate_predictions/' all_scans(i).name]));
         upper_range = max(max(max(this_scan_data)));
%         lower_range = min(min(min(abs(this_scan_data))));
upper_range = 5.97; % For matching prediction scales
lower_range = 0; % data already thresholded
        
        % spm_smooth(all_scans(i).name,['s5x_' all_scans(i).name],[5 0 0])
        % this_smoothed_scan_data = spm_read_vols(spm_vol(['s5x_' all_scans(i).name]));
        % upper_range_smooth = max(max(max(this_scan_data)));
        % lower_range_smooth = min(min(min(abs(this_scan_data))));
        
        %cfg.threshold = [3.3748 6.5]; %p=0.001
        cfg.threshold = [lower_range upper_range]; %scaled to data
        % cfg.threshold = [lower_range_smooth upper_range_smooth]; %scaled to data
        
        jp_spm8_surfacerender2_version_tc(['./multivariate_predictions/' all_scans(i).name],'hot',cfg)
        
        savepath = ['./rendered_images/' all_scans(i).name(1:end-4) '_' num2str(lower_range,'%.2f') '_' num2str(upper_range,'%.2f') '_' num2str(cfg.sampling_distance) 'voxelsampling_rh'];
        savepath = strrep(savepath,'>','-');
        %eval(['export_fig ' savepath '.png -transparent -m2'])
        eval(['export_fig ' savepath '.png -transparent'])
        close all
    end
end

all_scans = dir('./physio_physio/*.nii');

for sampling_distance = [3];
    cfg.sampling_distance = sampling_distance; % tc modification - When mesh is lower resolution than overlay image this can be problematic, now displays largest value within a Euclidean distance specified in number of voxels. Can be slow for high resolution images
    for i = 1:length(all_scans)
        
        this_scan_data = spm_read_vols(spm_vol(['./physio_physio/' all_scans(i).name]));
         upper_range = max(max(max(this_scan_data)));
%         lower_range = min(min(min(abs(this_scan_data))));
upper_range = 9.17; % For matching prediction scales
lower_range = 0; % data already thresholded
        
        % spm_smooth(all_scans(i).name,['s5x_' all_scans(i).name],[5 0 0])
        % this_smoothed_scan_data = spm_read_vols(spm_vol(['s5x_' all_scans(i).name]));
        % upper_range_smooth = max(max(max(this_scan_data)));
        % lower_range_smooth = min(min(min(abs(this_scan_data))));
        
        %cfg.threshold = [3.3748 6.5]; %p=0.001
        cfg.threshold = [lower_range upper_range]; %scaled to data
        % cfg.threshold = [lower_range_smooth upper_range_smooth]; %scaled to data
        
        jp_spm8_surfacerender2_version_tc(['./physio_physio/' all_scans(i).name],'hot',cfg)
        
        savepath = ['./rendered_images/' all_scans(i).name(1:end-4) '_' num2str(lower_range,'%.2f') '_' num2str(upper_range,'%.2f') '_' num2str(cfg.sampling_distance) 'voxelsampling'];
        savepath = strrep(savepath,'>','-');
        %eval(['export_fig ' savepath '.png -transparent -m2'])
        eval(['export_fig ' savepath '.png -transparent'])
        close all
    end
end


% Now plot PPI with ROIs in blue

for sampling_distance = [3];
    cfg.sampling_distance = sampling_distance; % tc modification - When mesh is lower resolution than overlay image this can be problematic, now displays largest value within a Euclidean distance specified in number of voxels. Can be slow for high resolution images
    cfg.overwrite = 1;
    upper_range = 9.17; % For matching prediction scales
    lower_range = 0; % data already thresholded
    
    cfg.threshold = [lower_range upper_range];
    cfg.threshold2 = [0 1]; %Binary ROI
    jp_spm8_surfacerender2_version_tc(['./ROI_scans/PPI-PrG-STG-noseed.nii'],'hot',cfg,['./ROI_scans/Left_Precentral_Univariate_Interaction_combined.nii'],repmat([0 0 0.5],256,1))
    
    savepath = ['./rendered_images/PPI-PrG-STG-noseed_' num2str(lower_range,'%.2f') '_' num2str(upper_range,'%.2f') '_' num2str(cfg.sampling_distance) 'voxelsampling'];
    savepath = strrep(savepath,'>','-');
    %eval(['export_fig ' savepath '.png -transparent -m2'])
    eval(['export_fig ' savepath '.png -transparent'])
    close all
    
    jp_spm8_surfacerender2_version_tc(['./physio_physio/PPI-Precentral-STG-Cluster.nii'],'hot',cfg,['./ROI_scans/Left_Precentral_Univariate_Interaction_combined.nii'],repmat([0 0 0.5],256,1))
    
    savepath = ['./rendered_images/PPI-PrG-STG-noseed_' num2str(lower_range,'%.2f') '_' num2str(upper_range,'%.2f') '_' num2str(cfg.sampling_distance) 'voxelsampling'];
    savepath = strrep(savepath,'>','-');
    %eval(['export_fig ' savepath '.png -transparent -m2'])
    eval(['export_fig ' savepath '.png -transparent'])
    close all
    
    jp_spm8_surfacerender2_version_tc(['./physio_physio/PPI-STG-Precentral-cluster.nii'],'hot',cfg,['./ROI_scans/Left_STG_Univariate3mm_15>3.nii'],repmat([0 0 0.5],256,1))
    
    savepath = ['./rendered_images/PPI-STG-PrG-noseed_' num2str(lower_range,'%.2f') '_' num2str(upper_range,'%.2f') '_' num2str(cfg.sampling_distance) 'voxelsampling'];
    savepath = strrep(savepath,'>','-');
    %eval(['export_fig ' savepath '.png -transparent -m2'])
    eval(['export_fig ' savepath '.png -transparent'])
    close all
    
    spm_imcalc(char({'./ROI_scans/Left_STG_Univariate3mm_15>3.nii';'./ROI_scans/Left_Precentral_Univariate_Interaction_combined.nii'}),['./ROI_scans/PrG+STG.nii'],['i1>0|i2>0'])
    jp_spm8_surfacerender2_version_tc(['./physio_physio/PPI-PrecentralxSTG-Negative-Cluster.nii'],'hot',cfg,['./ROI_scans/PrG+STG.nii'],repmat([0 0 0.5],256,1))
    
    savepath = ['./rendered_images/PPI-Interaction-noseed_' num2str(lower_range,'%.2f') '_' num2str(upper_range,'%.2f') '_' num2str(cfg.sampling_distance) 'voxelsampling'];
    savepath = strrep(savepath,'>','-');
    %eval(['export_fig ' savepath '.png -transparent -m2'])
    eval(['export_fig ' savepath '.png -transparent'])
    close all
        
    cfg.overwrite = 0;
end

% 
% %Now create ROI figures for illustration
% im1 = './ROI_scans/Left_STG_Univariate3mm_15>3.nii';
% im2 = './ROI_scans/Left_Precentral_Univariate_Interaction_combined.nii';
% cfg.threshold = [0 1];
% 
% jp_spm8_surfacerender2_version_tc(im1,'hot',cfg,im2,'cool')
