function all_data = module_extract_searchlight_volume(GLMDir,subjects,secondlevelnames,downsamp_ratio,masknames,group)
%Function for extracting the single subject RSA results from within a
%searchlight volume

addpath ./plotting

do_smoothed_maps = 0;  % if you want to do smoothed maps change this
if ~exist('downsamp_ratio','var')
    downsamp_ratio = 1;
end

versionCurrent = 'spearman';

% Gather images for first subject
images = {};
for m = 1:length(secondlevelnames)
    if downsamp_ratio == 1
        images{m} = spm_select('FPList', [GLMDir '/TDTcrossnobis/' versionCurrent '/'], ['^whireseffect-map_' secondlevelnames{m} '.nii']);
    else
        images{m} = spm_select('FPList', [GLMDir '/TDTcrossnobis_downsamp_' num2str(downsamp_ratio) '/' versionCurrent '/'], ['^whireseffect-map_' secondlevelnames{m} '.nii']);
    end
end

% masknames = {
%     %     'Left_Precentral_Univariate_Interaction_combined';
%     %   'Left_STG_Univariate3mm_15>3'
%     %'Left_Frontal_Univariate_MM>M'
%     %'Left_Temporal_Univariate_MM>M'
%     'Left_Angular_Univariate_Interaction_combined'
%     };
% secondlevelnames = {
%     'Clarity Congruency Interaction Positive';
%     };

all_data = [];
for this_mask = 1:length(masknames)
    for this_second_level = 1:length(secondlevelnames)
        this_scan = {};
        Y = spm_read_vols(spm_vol(['./atlas_Neuromorphometrics/' masknames{this_mask} '.nii']),1);
        indx = find(Y>0);
        [x,y,z] = ind2sub(size(Y),indx);
        XYZ = [x y z]';
        
        for crun = 1:size(subjects,2)
            disp(['Working on subject ' subjects{crun}]);
            
            condition_name = strsplit(images{this_second_level},'whireseffect-map_');
            condition_name = condition_name{2}(1:end-4);
            
            this_scan(end+1) = cellstr(strrep(images{this_second_level},subjects{1},subjects{crun}));
            
            all_data(crun,this_second_level,this_mask) = mean(spm_get_data(this_scan(end), XYZ),2);
        end
    end
    
    figure
    barweb([mean(all_data(group==1,:,this_mask));mean(all_data(group==2,:,this_mask))],[std(all_data(group==1,:,this_mask))/sqrt(sum(group==1));std(all_data(group==2,:,this_mask))/sqrt(sum(group==1))],[],{'Controls','Patients'},['Extracted searchlight RSAs in ' masknames{this_mask}],[],'Spearman correlation',[],[],secondlevelnames) ;
    ylim([(min(mean(all_data(:,:,this_mask)))-4*max(std(all_data(:,:,this_mask))/sqrt(sum(group==2)))),(max(mean(all_data(:,:,this_mask)))+4*max(std(all_data(:,:,this_mask))/sqrt(sum(group==2))))])
    
    
end