function module_extract_avneuralRDM_cluster(GLMDir,downsamp_ratio,mask_names)
%For taking already calculated crossnobis distances and doing RSA

if ~exist('downsamp_ratio','var')
    downsamp_ratio = 1;
end

addpath('/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/RSA_scripts/es_scripts_fMRI')
addpath('/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/RSA_scripts/decoding_toolbox_v3.999')
addpath(genpath('/group/language/data/ediz.sohoglu/matlab/rsatoolbox'));

%Define input data location
if downsamp_ratio == 1
    cfg.results.dir = fullfile(GLMDir,'TDTcrossnobis');
else
    cfg.results.dir = fullfile(GLMDir,['TDTcrossnobis_downsamp_' num2str(downsamp_ratio)]);
end

version = 'spearman'; % how to assess accuracy of model RDMs (pearson, spearman, weighted average)

outputDir = fullfile(cfg.results.dir,version,'extracted_neuralRDMs');
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

load(fullfile(cfg.results.dir,'res_other_average.mat'));
data = results.other_average.output;
notempty_data = find(~cellfun(@isempty,results.other_average.output));

for this_mask = 1:length(mask_names)
V = spm_vol(fullfile(GLMDir,'mask.nii')); % Searchlight mask
mask = spm_read_vols(V);
mask_index = results.mask_index;

V2 = spm_vol(fullfile(GLMDir,['../' mask_names{this_mask}])); % Searchlight mask
extraction_mask = spm_read_vols(V2);
mask_index = results.mask_index;

assert(all(all(V.mat==V2.mat)),'Extraction mask and searchlight data are not in the same space')
voxels_to_average = extraction_mask&mask;

RDMs_to_average{this_mask} = [];
for vx=1:numel(data)
    if voxels_to_average(mask_index(vx)) && ~isempty(data{vx})
        if ~isempty(RDMs_to_average{this_mask})
            RDMs_to_average{this_mask}(:,:,end+1) = data{vx};
        else
            RDMs_to_average{this_mask}(:,:,1) = data{vx};
        end
    end
end

average_RDM = squeeze(mean(RDMs_to_average{this_mask},3));
save(fullfile(outputDir,['average_RDM_' mask_names{this_mask}(1:end-4)]), 'average_RDM')


end
% 
%     modelRDM = vectorizeRDMs(models{m})';
%     effectMap = NaN(size(mask));
%     for vx=1:numel(data)
%         neuralRDM = vectorizeRDMs(data{vx})';
%         if isempty(neuralRDM)
%             continue
%         end
%         notempty = vx;
%         if ~isempty(strfind(version,'pearson'))
%             effectMap(mask_index(vx)) = fisherTransform(corr(modelRDM,neuralRDM,'type','Pearson','Rows','pairwise'));
%         elseif ~isempty(strfind(version,'spearman'))
%             effectMap(mask_index(vx)) = fisherTransform(corr(modelRDM,neuralRDM,'type','Spearman','Rows','pairwise'));
%         elseif ~isempty(strfind(version,'average'))
%             %effectMap(mask_index(vx)) = mean(neuralRDM(find(~isnan(modelRDM)),:),1);
%             effectMap(mask_index(vx)) = mean(neuralRDM(find(modelRDM==1),:),1);
%         end
%         if ~mod(vx,100)
%             disp(['Processing voxel ' num2str(vx) ' of ' num2str(numel(data))])
%         end
%     end
%     dims = size(effectMap);
%     downsamped_effectMap = effectMap(1:downsamp_ratio:dims(1),1:downsamp_ratio:dims(2),1:downsamp_ratio:dims(3));
%     downsamped_V.mat = V.mat;
%     downsamped_V.mat(1:3,1:3)=downsamped_V.mat(1:3,1:3)*downsamp_ratio;
%     
%     saveMRImage(downsamped_effectMap,fullfile(outputDir,['xzcxeffect-map_' this_model_name{m} '.nii']),downsamped_V.mat);
% end
% end
% 
% %% Make effect-maps (by correlating neural RDMs to model RDMs)
% 
% clear models
% 
% basemodels.vowels = zeros(16,16);
% basemodels.vowels(1:17:end) = 1;
% basemodels.vowels(2:68:end) = 1/3;
% basemodels.vowels(3:68:end) = 1/3;
% basemodels.vowels(4:68:end) = 1/3;
% basemodels.vowels(17:68:end) = 1/3;
% basemodels.vowels(19:68:end) = 1/3;
% basemodels.vowels(20:68:end) = 1/3;
% basemodels.vowels(33:68:end) = 1/3;
% basemodels.vowels(34:68:end) = 1/3;
% basemodels.vowels(36:68:end) = 1/3;
% basemodels.vowels(49:68:end) = 1/3;
% basemodels.vowels(50:68:end) = 1/3;
% basemodels.vowels(51:68:end) = 1/3;
% basemodels.vowels = 1-basemodels.vowels;
% 
% %Squares based on all shared features
% 
% basemodels.shared_segments = zeros(16,16);
% basemodels.shared_segments(1:17:end) = 1;
% basemodels.shared_segments(2:68:end) = 2/3;
% basemodels.shared_segments(3:68:end) = 2/3;
% basemodels.shared_segments(4:68:end) = 1/3;
% basemodels.shared_segments(17:68:end) = 2/3;
% basemodels.shared_segments(19:68:end) = 1/3;
% basemodels.shared_segments(20:68:end) = 2/3;
% basemodels.shared_segments(33:68:end) = 2/3;
% basemodels.shared_segments(34:68:end) = 1/3;
% basemodels.shared_segments(36:68:end) = 2/3;
% basemodels.shared_segments(49:68:end) = 1/3;
% basemodels.shared_segments(50:68:end) = 2/3;
% basemodels.shared_segments(51:68:end) = 2/3;
% 
% basemodels.shared_segments(1,16) = 1/3;
% basemodels.shared_segments(1,14) = 1/3;
% basemodels.shared_segments(16,1) = 1/3;
% basemodels.shared_segments(14,1) = 1/3;
% basemodels.shared_segments(3,16) = 1/3;
% basemodels.shared_segments(3,14) = 1/3;
% basemodels.shared_segments(16,3) = 1/3;
% basemodels.shared_segments(14,3) = 1/3;
% 
% basemodels.shared_segments(5,9) = 1/3;
% basemodels.shared_segments(7,9) = 1/3;
% basemodels.shared_segments(9,5) = 1/3;
% basemodels.shared_segments(9,7) = 1/3;
% basemodels.shared_segments(5,11) = 1/3;
% basemodels.shared_segments(7,11) = 1/3;
% basemodels.shared_segments(11,5) = 1/3;
% basemodels.shared_segments(11,7) = 1/3;
% 
% basemodels.shared_segments(6,10) = 1/3;
% basemodels.shared_segments(8,10) = 1/3;
% basemodels.shared_segments(10,6) = 1/3;
% basemodels.shared_segments(10,8) = 1/3;
% basemodels.shared_segments(6,12) = 1/3;
% basemodels.shared_segments(8,12) = 1/3;
% basemodels.shared_segments(12,6) = 1/3;
% basemodels.shared_segments(12,8) = 1/3;
% 
% basemodels.shared_segments(15,3) = 1/3;
% basemodels.shared_segments(15,4) = 1/3;
% basemodels.shared_segments(16,4) = 1/3;
% basemodels.shared_segments(16,3) = 1/3;
% basemodels.shared_segments(3,16) = 1/3;
% basemodels.shared_segments(4,16) = 1/3;
% basemodels.shared_segments(4,15) = 1/3;
% basemodels.shared_segments(3,15) = 1/3;
% 
% basemodels.shared_segments = 1-basemodels.shared_segments;
% 
% basemodelNames = {'vowels','shared_segments'};
% 
% % Set the label names to the regressor names which you want to use for
% % your similarity analysis, e.g.
% %labelnames = {'Strong+M_Set1_Item1','Strong+M_Set1_Item2','Strong+M_Set1_Item3','Strong+M_Set1_Item4','Strong+M_Set2_Item1','Strong+M_Set2_Item2','Strong+M_Set2_Item3','Strong+M_Set2_Item4','Strong+M_Set3_Item1','Strong+M_Set3_Item2','Strong+M_Set3_Item3','Strong+M_Set3_Item4','Strong+M_Set4_Item1','Strong+M_Set4_Item2','Strong+M_Set4_Item3','Strong+M_Set4_Item4','Strong+M_Set5_Item1','Strong+M_Set5_Item2','Strong+M_Set5_Item3','Strong+M_Set5_Item4','Strong+M_Set6_Item1','Strong+M_Set6_Item2','Strong+M_Set6_Item3','Strong+M_Set6_Item4','Strong+M_Set7_Item1','Strong+M_Set7_Item2','Strong+M_Set7_Item3','Strong+M_Set7_Item4','Strong+M_Set8_Item1','Strong+M_Set8_Item2','Strong+M_Set8_Item3','Strong+M_Set8_Item4','Weak+M_Set1_Item1','Weak+M_Set1_Item2','Weak+M_Set1_Item3','Weak+M_Set1_Item4','Weak+M_Set2_Item1','Weak+M_Set2_Item2','Weak+M_Set2_Item3','Weak+M_Set2_Item4','Weak+M_Set3_Item1','Weak+M_Set3_Item2','Weak+M_Set3_Item3','Weak+M_Set3_Item4','Weak+M_Set4_Item1','Weak+M_Set4_Item2','Weak+M_Set4_Item3','Weak+M_Set4_Item4','Weak+M_Set5_Item1','Weak+M_Set5_Item2','Weak+M_Set5_Item3','Weak+M_Set5_Item4','Weak+M_Set6_Item1','Weak+M_Set6_Item2','Weak+M_Set6_Item3','Weak+M_Set6_Item4','Weak+M_Set7_Item1','Weak+M_Set7_Item2','Weak+M_Set7_Item3','Weak+M_Set7_Item4','Weak+M_Set8_Item1','Weak+M_Set8_Item2','Weak+M_Set8_Item3','Weak+M_Set8_Item4','Strong+MM_Set1_Item1','Strong+MM_Set1_Item2','Strong+MM_Set1_Item3','Strong+MM_Set1_Item4','Strong+MM_Set2_Item1','Strong+MM_Set2_Item2','Strong+MM_Set2_Item3','Strong+MM_Set2_Item4','Strong+MM_Set3_Item1','Strong+MM_Set3_Item2','Strong+MM_Set3_Item3','Strong+MM_Set3_Item4','Strong+MM_Set4_Item1','Strong+MM_Set4_Item2','Strong+MM_Set4_Item3','Strong+MM_Set4_Item4','Strong+MM_Set5_Item1','Strong+MM_Set5_Item2','Strong+MM_Set5_Item3','Strong+MM_Set5_Item4','Strong+MM_Set6_Item1','Strong+MM_Set6_Item2','Strong+MM_Set6_Item3','Strong+MM_Set6_Item4','Strong+MM_Set7_Item1','Strong+MM_Set7_Item2','Strong+MM_Set7_Item3','Strong+MM_Set7_Item4','Strong+MM_Set8_Item1','Strong+MM_Set8_Item2','Strong+MM_Set8_Item3','Strong+MM_Set8_Item4','Weak+MM_Set1_Item1','Weak+MM_Set1_Item2','Weak+MM_Set1_Item3','Weak+MM_Set1_Item4','Weak+MM_Set2_Item1','Weak+MM_Set2_Item2','Weak+MM_Set2_Item3','Weak+MM_Set2_Item4','Weak+MM_Set3_Item1','Weak+MM_Set3_Item2','Weak+MM_Set3_Item3','Weak+MM_Set3_Item4','Weak+MM_Set4_Item1','Weak+MM_Set4_Item2','Weak+MM_Set4_Item3','Weak+MM_Set4_Item4','Weak+MM_Set5_Item1','Weak+MM_Set5_Item2','Weak+MM_Set5_Item3','Weak+MM_Set5_Item4','Weak+MM_Set6_Item1','Weak+MM_Set6_Item2','Weak+MM_Set6_Item3','Weak+MM_Set6_Item4','Weak+MM_Set7_Item1','Weak+MM_Set7_Item2','Weak+MM_Set7_Item3','Weak+MM_Set7_Item4','Weak+MM_Set8_Item1','Weak+MM_Set8_Item2','Weak+MM_Set8_Item3','Weak+MM_Set8_Item4','Strong+Noise_Set1_Item1','Strong+Noise_Set1_Item2','Strong+Noise_Set1_Item3','Strong+Noise_Set1_Item4','Strong+Noise_Set2_Item1','Strong+Noise_Set2_Item2','Strong+Noise_Set2_Item3','Strong+Noise_Set2_Item4','Strong+Noise_Set3_Item1','Strong+Noise_Set3_Item2','Strong+Noise_Set3_Item3','Strong+Noise_Set3_Item4','Strong+Noise_Set4_Item1','Strong+Noise_Set4_Item2','Strong+Noise_Set4_Item3','Strong+Noise_Set4_Item4','Strong+Noise_Set5_Item1','Strong+Noise_Set5_Item2','Strong+Noise_Set5_Item3','Strong+Noise_Set5_Item4','Strong+Noise_Set6_Item1','Strong+Noise_Set6_Item2','Strong+Noise_Set6_Item3','Strong+Noise_Set6_Item4','Strong+Noise_Set7_Item1','Strong+Noise_Set7_Item2','Strong+Noise_Set7_Item3','Strong+Noise_Set7_Item4','Strong+Noise_Set8_Item1','Strong+Noise_Set8_Item2','Strong+Noise_Set8_Item3','Strong+Noise_Set8_Item4','Weak+Noise_Set1_Item1','Weak+Noise_Set1_Item2','Weak+Noise_Set1_Item3','Weak+Noise_Set1_Item4','Weak+Noise_Set2_Item1','Weak+Noise_Set2_Item2','Weak+Noise_Set2_Item3','Weak+Noise_Set2_Item4','Weak+Noise_Set3_Item1','Weak+Noise_Set3_Item2','Weak+Noise_Set3_Item3','Weak+Noise_Set3_Item4','Weak+Noise_Set4_Item1','Weak+Noise_Set4_Item2','Weak+Noise_Set4_Item3','Weak+Noise_Set4_Item4','Weak+Noise_Set5_Item1','Weak+Noise_Set5_Item2','Weak+Noise_Set5_Item3','Weak+Noise_Set5_Item4','Weak+Noise_Set6_Item1','Weak+Noise_Set6_Item2','Weak+Noise_Set6_Item3','Weak+Noise_Set6_Item4','Weak+Noise_Set7_Item1','Weak+Noise_Set7_Item2','Weak+Noise_Set7_Item3','Weak+Noise_Set7_Item4','Weak+Noise_Set8_Item1','Weak+Noise_Set8_Item2','Weak+Noise_Set8_Item3','Weak+Noise_Set8_Item4','Noise+Speech_Set1_Item1','Noise+Speech_Set1_Item2','Noise+Speech_Set1_Item3','Noise+Speech_Set1_Item4','Noise+Speech_Set2_Item1','Noise+Speech_Set2_Item2','Noise+Speech_Set2_Item3','Noise+Speech_Set2_Item4','Noise+Speech_Set3_Item1','Noise+Speech_Set3_Item2','Noise+Speech_Set3_Item3','Noise+Speech_Set3_Item4','Noise+Speech_Set4_Item1','Noise+Speech_Set4_Item2','Noise+Speech_Set4_Item3','Noise+Speech_Set4_Item4','Noise+Speech_Set5_Item1','Noise+Speech_Set5_Item2','Noise+Speech_Set5_Item3','Noise+Speech_Set5_Item4','Noise+Speech_Set6_Item1','Noise+Speech_Set6_Item2','Noise+Speech_Set6_Item3','Noise+Speech_Set6_Item4','Noise+Speech_Set7_Item1','Noise+Speech_Set7_Item2','Noise+Speech_Set7_Item3','Noise+Speech_Set7_Item4','Noise+Speech_Set8_Item1','Noise+Speech_Set8_Item2','Noise+Speech_Set8_Item3','Noise+Speech_Set8_Item4'};
% temp = load([GLMDir filesep 'SPM.mat']);
% labelnames = {};
% for i = 1:length(temp.SPM.Sess(1).U)
%     if ~strncmp(temp.SPM.Sess(1).U(i).name,{'Match','Mismatch','Written'},5)
%         continue
%     else
%         labelnames(end+1) = temp.SPM.Sess(1).U(i).name;
%     end
% end
% labels = 1:length(labelnames);
% 
% modeltemplate = NaN(size(results.other_average.output{notempty_data(1)}));
% 
% labelnames_denumbered = {};
% for i = 1:length(labelnames)
%     labelnames_denumbered{i} = labelnames{i}(isletter(labelnames{i})|isspace(labelnames{i}));
% end
% modelNames = unique(labelnames_denumbered,'stable');
% 
% for j = 1:length(basemodelNames)
%     for i = 1:length(modelNames)
%         this_model = ((j-1)*length(modelNames))+i;
%         models{this_model} = modeltemplate;
%         models{this_model}(strcmp(modelNames{i},labelnames_denumbered),strcmp(modelNames{i},labelnames_denumbered))=basemodels.(basemodelNames{j});
%         this_model_name{this_model} = [modelNames{i} ' ' basemodelNames{j}];
%         %Optional check - view matrix
%         %         imagesc(models{this_model},'AlphaData',~isnan(models{this_model}))
%         %         title(this_model_name{this_model})
%         %         pause
%     end
% end
% 
% MisMatch_Cross_decode_base = zeros(16,16);
% MisMatch_Cross_decode_base(9:17:end/2) = 1;
% MisMatch_Cross_decode_base(end/2+1:17:end) = 1;
% MisMatch_Cross_decode_base = 1-MisMatch_Cross_decode_base;
% 
% cross_decode_label_pairs = {
%     'Match Unclear', 'Mismatch Unclear';
%     'Match Clear', 'Mismatch Unclear';
%     'Match Unclear', 'Mismatch Clear';
%     'Match Clear', 'Mismatch Clear'};
% 
% for i = 1:size(cross_decode_label_pairs,1)
%     models{end+1} = modeltemplate;
%     models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = MisMatch_Cross_decode_base;
%     models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = MisMatch_Cross_decode_base;
%     this_model_name{end+1} = [cross_decode_label_pairs{i,1} ' to ' cross_decode_label_pairs{i,2} ' Cross-decode'];
%     %Optional check - view matrix
%     %             imagesc(models{end},'AlphaData',~isnan(models{end}))
%     %             title(this_model_name{end})
%     %             pause
% end
% 
% %Now attempt cross-condition shared segments RSA without cross decoding, recognising that the MisMatch cue
% %was consistently 8 elements after/before the auditory word
% basemodels.shared_segments_cross = circshift(basemodels.shared_segments,[8 0]); %I think this is correct, but need to 100% check the off-diagonals
% for i = 1:size(cross_decode_label_pairs,1)
%     models{end+1} = modeltemplate;
%     models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.shared_segments_cross;
%     models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.shared_segments_cross';
%     this_model_name{end+1} = [cross_decode_label_pairs{i,1} ' to ' cross_decode_label_pairs{i,2} ' Shared Segments - cross'];
%     %Optional check - view matrix
%     %                     imagesc(models{end},'AlphaData',~isnan(models{end}))
%     %                     title(this_model_name{end})
%     %                     pause
% end
% 
% 
% %Now attempt cross-condition shared segments RSA without cross decoding, recognising that the MisMatch cue
% %was consistently 8 elements after/before the auditory word
% basemodels.shared_segments_cross_noself = basemodels.shared_segments;
% basemodels.shared_segments_cross_noself(1:17:end) = NaN;
% basemodels.shared_segments_cross_noself = circshift(basemodels.shared_segments_cross_noself,[8 0]);
% for i = 1:size(cross_decode_label_pairs,1)
%     models{end+1} = modeltemplate;
%     models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.shared_segments_cross_noself;
%     models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.shared_segments_cross_noself';
%     this_model_name{end+1} = [cross_decode_label_pairs{i,1} ' to ' cross_decode_label_pairs{i,2} ' Shared Segments - no self'];
%     %Optional check - view matrix
%     %             imagesc(models{end},'AlphaData',~isnan(models{end}))
%     %             title(this_model_name{end})
%     %             pause
% end