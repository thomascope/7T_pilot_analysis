function module_make_effect_maps(GLMDir,downsamp_ratio)
%For taking already calculated crossnobis distances and doing RSA

redo_maps = 0; %If you want to calculate them again for some reason.

if ~exist('downsamp_ratio','var')
    downsamp_ratio = 1;
end

addpath('/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/RSA_scripts/es_scripts_fMRI')
addpath('/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/RSA_scripts/decoding_toolbox_v3.999')
addpath(genpath('/group/language/data/ediz.sohoglu/matlab/rsatoolbox'));

%Define input data location
if downsamp_ratio == 1
    %cfg.results.dir = fullfile(GLMDir,'TDTcrossnobis');
    cfg.results.dir = fullfile(GLMDir,'TDTcrossnobis_parallel');
else
    cfg.results.dir = fullfile(GLMDir,['TDTcrossnobis_downsamp_' num2str(downsamp_ratio)]);
end

% Set the label names to the regressor names which you want to use for
% your similarity analysis, e.g.
%labelnames = {'Strong+M_Set1_Item1','Strong+M_Set1_Item2','Strong+M_Set1_Item3','Strong+M_Set1_Item4','Strong+M_Set2_Item1','Strong+M_Set2_Item2','Strong+M_Set2_Item3','Strong+M_Set2_Item4','Strong+M_Set3_Item1','Strong+M_Set3_Item2','Strong+M_Set3_Item3','Strong+M_Set3_Item4','Strong+M_Set4_Item1','Strong+M_Set4_Item2','Strong+M_Set4_Item3','Strong+M_Set4_Item4','Strong+M_Set5_Item1','Strong+M_Set5_Item2','Strong+M_Set5_Item3','Strong+M_Set5_Item4','Strong+M_Set6_Item1','Strong+M_Set6_Item2','Strong+M_Set6_Item3','Strong+M_Set6_Item4','Strong+M_Set7_Item1','Strong+M_Set7_Item2','Strong+M_Set7_Item3','Strong+M_Set7_Item4','Strong+M_Set8_Item1','Strong+M_Set8_Item2','Strong+M_Set8_Item3','Strong+M_Set8_Item4','Weak+M_Set1_Item1','Weak+M_Set1_Item2','Weak+M_Set1_Item3','Weak+M_Set1_Item4','Weak+M_Set2_Item1','Weak+M_Set2_Item2','Weak+M_Set2_Item3','Weak+M_Set2_Item4','Weak+M_Set3_Item1','Weak+M_Set3_Item2','Weak+M_Set3_Item3','Weak+M_Set3_Item4','Weak+M_Set4_Item1','Weak+M_Set4_Item2','Weak+M_Set4_Item3','Weak+M_Set4_Item4','Weak+M_Set5_Item1','Weak+M_Set5_Item2','Weak+M_Set5_Item3','Weak+M_Set5_Item4','Weak+M_Set6_Item1','Weak+M_Set6_Item2','Weak+M_Set6_Item3','Weak+M_Set6_Item4','Weak+M_Set7_Item1','Weak+M_Set7_Item2','Weak+M_Set7_Item3','Weak+M_Set7_Item4','Weak+M_Set8_Item1','Weak+M_Set8_Item2','Weak+M_Set8_Item3','Weak+M_Set8_Item4','Strong+MM_Set1_Item1','Strong+MM_Set1_Item2','Strong+MM_Set1_Item3','Strong+MM_Set1_Item4','Strong+MM_Set2_Item1','Strong+MM_Set2_Item2','Strong+MM_Set2_Item3','Strong+MM_Set2_Item4','Strong+MM_Set3_Item1','Strong+MM_Set3_Item2','Strong+MM_Set3_Item3','Strong+MM_Set3_Item4','Strong+MM_Set4_Item1','Strong+MM_Set4_Item2','Strong+MM_Set4_Item3','Strong+MM_Set4_Item4','Strong+MM_Set5_Item1','Strong+MM_Set5_Item2','Strong+MM_Set5_Item3','Strong+MM_Set5_Item4','Strong+MM_Set6_Item1','Strong+MM_Set6_Item2','Strong+MM_Set6_Item3','Strong+MM_Set6_Item4','Strong+MM_Set7_Item1','Strong+MM_Set7_Item2','Strong+MM_Set7_Item3','Strong+MM_Set7_Item4','Strong+MM_Set8_Item1','Strong+MM_Set8_Item2','Strong+MM_Set8_Item3','Strong+MM_Set8_Item4','Weak+MM_Set1_Item1','Weak+MM_Set1_Item2','Weak+MM_Set1_Item3','Weak+MM_Set1_Item4','Weak+MM_Set2_Item1','Weak+MM_Set2_Item2','Weak+MM_Set2_Item3','Weak+MM_Set2_Item4','Weak+MM_Set3_Item1','Weak+MM_Set3_Item2','Weak+MM_Set3_Item3','Weak+MM_Set3_Item4','Weak+MM_Set4_Item1','Weak+MM_Set4_Item2','Weak+MM_Set4_Item3','Weak+MM_Set4_Item4','Weak+MM_Set5_Item1','Weak+MM_Set5_Item2','Weak+MM_Set5_Item3','Weak+MM_Set5_Item4','Weak+MM_Set6_Item1','Weak+MM_Set6_Item2','Weak+MM_Set6_Item3','Weak+MM_Set6_Item4','Weak+MM_Set7_Item1','Weak+MM_Set7_Item2','Weak+MM_Set7_Item3','Weak+MM_Set7_Item4','Weak+MM_Set8_Item1','Weak+MM_Set8_Item2','Weak+MM_Set8_Item3','Weak+MM_Set8_Item4','Strong+Noise_Set1_Item1','Strong+Noise_Set1_Item2','Strong+Noise_Set1_Item3','Strong+Noise_Set1_Item4','Strong+Noise_Set2_Item1','Strong+Noise_Set2_Item2','Strong+Noise_Set2_Item3','Strong+Noise_Set2_Item4','Strong+Noise_Set3_Item1','Strong+Noise_Set3_Item2','Strong+Noise_Set3_Item3','Strong+Noise_Set3_Item4','Strong+Noise_Set4_Item1','Strong+Noise_Set4_Item2','Strong+Noise_Set4_Item3','Strong+Noise_Set4_Item4','Strong+Noise_Set5_Item1','Strong+Noise_Set5_Item2','Strong+Noise_Set5_Item3','Strong+Noise_Set5_Item4','Strong+Noise_Set6_Item1','Strong+Noise_Set6_Item2','Strong+Noise_Set6_Item3','Strong+Noise_Set6_Item4','Strong+Noise_Set7_Item1','Strong+Noise_Set7_Item2','Strong+Noise_Set7_Item3','Strong+Noise_Set7_Item4','Strong+Noise_Set8_Item1','Strong+Noise_Set8_Item2','Strong+Noise_Set8_Item3','Strong+Noise_Set8_Item4','Weak+Noise_Set1_Item1','Weak+Noise_Set1_Item2','Weak+Noise_Set1_Item3','Weak+Noise_Set1_Item4','Weak+Noise_Set2_Item1','Weak+Noise_Set2_Item2','Weak+Noise_Set2_Item3','Weak+Noise_Set2_Item4','Weak+Noise_Set3_Item1','Weak+Noise_Set3_Item2','Weak+Noise_Set3_Item3','Weak+Noise_Set3_Item4','Weak+Noise_Set4_Item1','Weak+Noise_Set4_Item2','Weak+Noise_Set4_Item3','Weak+Noise_Set4_Item4','Weak+Noise_Set5_Item1','Weak+Noise_Set5_Item2','Weak+Noise_Set5_Item3','Weak+Noise_Set5_Item4','Weak+Noise_Set6_Item1','Weak+Noise_Set6_Item2','Weak+Noise_Set6_Item3','Weak+Noise_Set6_Item4','Weak+Noise_Set7_Item1','Weak+Noise_Set7_Item2','Weak+Noise_Set7_Item3','Weak+Noise_Set7_Item4','Weak+Noise_Set8_Item1','Weak+Noise_Set8_Item2','Weak+Noise_Set8_Item3','Weak+Noise_Set8_Item4','Noise+Speech_Set1_Item1','Noise+Speech_Set1_Item2','Noise+Speech_Set1_Item3','Noise+Speech_Set1_Item4','Noise+Speech_Set2_Item1','Noise+Speech_Set2_Item2','Noise+Speech_Set2_Item3','Noise+Speech_Set2_Item4','Noise+Speech_Set3_Item1','Noise+Speech_Set3_Item2','Noise+Speech_Set3_Item3','Noise+Speech_Set3_Item4','Noise+Speech_Set4_Item1','Noise+Speech_Set4_Item2','Noise+Speech_Set4_Item3','Noise+Speech_Set4_Item4','Noise+Speech_Set5_Item1','Noise+Speech_Set5_Item2','Noise+Speech_Set5_Item3','Noise+Speech_Set5_Item4','Noise+Speech_Set6_Item1','Noise+Speech_Set6_Item2','Noise+Speech_Set6_Item3','Noise+Speech_Set6_Item4','Noise+Speech_Set7_Item1','Noise+Speech_Set7_Item2','Noise+Speech_Set7_Item3','Noise+Speech_Set7_Item4','Noise+Speech_Set8_Item1','Noise+Speech_Set8_Item2','Noise+Speech_Set8_Item3','Noise+Speech_Set8_Item4'};
temp = load([GLMDir filesep 'SPM.mat']);
labelnames = {};
for i = 1:length(temp.SPM.Sess(1).U)
    if ~strncmp(temp.SPM.Sess(1).U(i).name,{'Match','Mismatch','Written'},5)
        continue
    else
        labelnames(end+1) = temp.SPM.Sess(1).U(i).name;
    end
end
labels = 1:length(labelnames);

%% Make effect-maps (by correlating neural RDMs to model RDMs)

version = 'spearman'; % how to assess accuracy of model RDMs (pearson, spearman, weighted average)

if downsamp_ratio == 1
    outputDir = fullfile(GLMDir,'TDTcrossnobis',version);
else
    outputDir = fullfile(GLMDir,['TDTcrossnobis_downsamp_' num2str(downsamp_ratio)],version);
end
if ~exist(outputDir,'dir') mkdir(outputDir); end

clear models
%Squares based on all shared features

basemodels.shared_segments = zeros(16,16);
basemodels.shared_segments(1:17:end) = 1;
basemodels.shared_segments(2:68:end) = 2/3;
basemodels.shared_segments(3:68:end) = 2/3;
basemodels.shared_segments(4:68:end) = 1/3;
basemodels.shared_segments(17:68:end) = 2/3;
basemodels.shared_segments(19:68:end) = 1/3;
basemodels.shared_segments(20:68:end) = 2/3;
basemodels.shared_segments(33:68:end) = 2/3;
basemodels.shared_segments(34:68:end) = 1/3;
basemodels.shared_segments(36:68:end) = 2/3;
basemodels.shared_segments(49:68:end) = 1/3;
basemodels.shared_segments(50:68:end) = 2/3;
basemodels.shared_segments(51:68:end) = 2/3;

basemodels.shared_segments(1,16) = 1/3;
basemodels.shared_segments(1,14) = 1/3;
basemodels.shared_segments(16,1) = 1/3;
basemodels.shared_segments(14,1) = 1/3;
basemodels.shared_segments(3,16) = 1/3;
basemodels.shared_segments(3,14) = 1/3;
basemodels.shared_segments(16,3) = 1/3;
basemodels.shared_segments(14,3) = 1/3;

basemodels.shared_segments(5,9) = 1/3;
basemodels.shared_segments(7,9) = 1/3;
basemodels.shared_segments(9,5) = 1/3;
basemodels.shared_segments(9,7) = 1/3;
basemodels.shared_segments(5,11) = 1/3;
basemodels.shared_segments(7,11) = 1/3;
basemodels.shared_segments(11,5) = 1/3;
basemodels.shared_segments(11,7) = 1/3;

basemodels.shared_segments(6,10) = 1/3;
basemodels.shared_segments(8,10) = 1/3;
basemodels.shared_segments(10,6) = 1/3;
basemodels.shared_segments(10,8) = 1/3;
basemodels.shared_segments(6,12) = 1/3;
basemodels.shared_segments(8,12) = 1/3;
basemodels.shared_segments(12,6) = 1/3;
basemodels.shared_segments(12,8) = 1/3;

basemodels.shared_segments(15,3) = 1/3;
basemodels.shared_segments(15,4) = 1/3;
basemodels.shared_segments(16,4) = 1/3;
basemodels.shared_segments(16,3) = 1/3;
basemodels.shared_segments(3,16) = 1/3;
basemodels.shared_segments(4,16) = 1/3;
basemodels.shared_segments(4,15) = 1/3;
basemodels.shared_segments(3,15) = 1/3;

basemodels.shared_segments = 1-basemodels.shared_segments;

basemodels.shared_segments_cross = circshift(basemodels.shared_segments,[8 0]); %I think this is correct, but need to 100% check the off-diagonals

basemodelNames = {'shared_segments'};

load(fullfile(cfg.results.dir,'res_other_average.mat'));
data = results.other_average.output;
notempty_data = find(~cellfun(@isempty,results.other_average.output));
modeltemplate = NaN(size(results.other_average.output{notempty_data(1)}));

labelnames_denumbered = {};
for i = 1:length(labelnames)
    labelnames_denumbered{i} = labelnames{i}(isletter(labelnames{i})|isspace(labelnames{i}));
end
modelNames = unique(labelnames_denumbered,'stable');

for j = 1:length(basemodelNames)
    for i = 1:length(modelNames)
        this_model = ((j-1)*length(modelNames))+i;
        models{this_model} = modeltemplate;
        models{this_model}(strcmp(modelNames{i},labelnames_denumbered),strcmp(modelNames{i},labelnames_denumbered))=basemodels.(basemodelNames{j});
        this_model_name{this_model} = [modelNames{i} ' ' basemodelNames{j}];
        %Optional check - view matrix
%                 imagesc(models{this_model},'AlphaData',~isnan(models{this_model}))
%                 title(this_model_name{this_model})
%                 pause
    end
end

%Now create combined Match and MisMatch within-condition RSA
% Now add combined conditions
combine_label_sets = {
    'Match Unclear', 'Mismatch Unclear';
    'Match Clear', 'Mismatch Clear';
    'Match Unclear', 'Match Clear';
    'Mismatch Unclear', 'Mismatch Clear';
    };

basemodelNames = {'shared_segments'};

for i = 1:size(combine_label_sets,1)
    for j = 1:length(basemodelNames)
        models{end+1} = modeltemplate;
        this_model_name{end+1} = [combine_label_sets{i,1} ' and ' combine_label_sets{i,2} ' ' basemodelNames{j}];
        models{end}(strcmp(combine_label_sets{i,1},labelnames_denumbered),strcmp(combine_label_sets{i,1},labelnames_denumbered)) = basemodels.(basemodelNames{j});
        models{end}(strcmp(combine_label_sets{i,2},labelnames_denumbered),strcmp(combine_label_sets{i,2},labelnames_denumbered)) = basemodels.(basemodelNames{j});
    end
end

combine_label_sets = {
    'Match Unclear', 'Mismatch Unclear', 'Match Clear', 'Mismatch Clear';
    };
for j = 1:length(basemodelNames)
    models{end+1} = modeltemplate;
    this_model_name{end+1} = ['All ' basemodelNames{j}];
    for i = 1:size(combine_label_sets,2)
        models{end}(strcmp(combine_label_sets{1,i},labelnames_denumbered),strcmp(combine_label_sets{1,i},labelnames_denumbered)) = basemodels.(basemodelNames{j});
    end
end

MisMatch_Cross_decode_base = zeros(16,16);
MisMatch_Cross_decode_base(9:17:end/2) = 1;
MisMatch_Cross_decode_base(end/2+1:17:end) = 1;
MisMatch_Cross_decode_base = 1-MisMatch_Cross_decode_base;

Match_Cross_decode_base = 1-eye(16);

cross_decode_label_pairs = {
    'Match Unclear', 'Mismatch Unclear';
    'Match Clear', 'Mismatch Unclear';
    'Match Unclear', 'Mismatch Clear';
    'Match Clear', 'Mismatch Clear';
    'Match Unclear', 'Match Clear';
    'Mismatch Unclear', 'Mismatch Clear';
    'Match Unclear', 'Written';
    'Match Clear', 'Written';
    'Mismatch Unclear', 'Written';
    'Mismatch Clear', 'Written'};

% for i = 1:size(cross_decode_label_pairs,1)
%     models{end+1} = modeltemplate;
%     models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = MisMatch_Cross_decode_base;
%     models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = MisMatch_Cross_decode_base';
%     this_model_name{end+1} = [cross_decode_label_pairs{i,1} ' to ' cross_decode_label_pairs{i,2} ' Cross-decode'];
%     %Optional check - view matrix
%     %             imagesc(models{end},'AlphaData',~isnan(models{end}))
%     %             title(this_model_name{end})
%     %             pause
% end

%Now attempt cross-condition shared segments RSA without cross decoding, recognising that the MisMatch cue
%was consistently 8 elements after/before the auditory word
for i = 1:size(cross_decode_label_pairs,1)
    models{end+1} = modeltemplate;
    models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.shared_segments_cross;
    models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.shared_segments_cross';
    this_model_name{end+1} = [cross_decode_label_pairs{i,1} ' to ' cross_decode_label_pairs{i,2} ' Shared Segments - cross'];
    %Optional check - view matrix
    %                     imagesc(models{end},'AlphaData',~isnan(models{end}))
    %                     title(this_model_name{end})
    %                     pause
end


%Now attempt cross-condition shared segments RSA without cross decoding, recognising that the MisMatch cue
%was consistently 8 elements after/before the auditory word
basemodels.shared_segments_cross_noself = basemodels.shared_segments;
basemodels.shared_segments_cross_noself(1:17:end) = NaN;
basemodels.shared_segments_cross_noself = circshift(basemodels.shared_segments_cross_noself,[8 0]);
for i = 1:size(cross_decode_label_pairs,1)
    models{end+1} = modeltemplate;
    models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.shared_segments_cross_noself;
    models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.shared_segments_cross_noself';
    this_model_name{end+1} = [cross_decode_label_pairs{i,1} ' to ' cross_decode_label_pairs{i,2} ' Shared Segments - no self'];
    %Optional check - view matrix
    %                 imagesc(models{end},'AlphaData',~isnan(models{end}))
    %                 title(this_model_name{end})
    %                 colorbar
    %                 pause
end

cross_decode_label_pairs = {
    'Match Unclear', 'Mismatch Unclear';
    'Match Clear', 'Mismatch Unclear';
    'Match Unclear', 'Mismatch Clear';
    'Match Clear', 'Mismatch Clear'
    'Match Unclear', 'Match Clear';
    'Mismatch Unclear', 'Mismatch Clear';
    'Match Unclear', 'Written';
    'Match Clear', 'Written';
    'Mismatch Unclear', 'Written';
    'Mismatch Clear', 'Written'
    'Match Unclear', 'Written';
    'Match Clear', 'Written';
    };

% for i = 1:size(cross_decode_label_pairs,1)
%     models{end+1} = modeltemplate;
%     models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = Match_Cross_decode_base;
%     models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = Match_Cross_decode_base';
%     this_model_name{end+1} = [cross_decode_label_pairs{i,1} ' to ' cross_decode_label_pairs{i,2} ' Cross-decode_Match'];
%     %Optional check - view matrix
%     %             imagesc(models{end},'AlphaData',~isnan(models{end}))
%     %             title(this_model_name{end})
%     %             pause
% end

% %Now attempt cross-condition shared segments RSA without cross decoding, recognising that the MisMatch cue
% %was consistently 8 elements after/before the auditory word
% for i = 1:size(cross_decode_label_pairs,1)
%     models{end+1} = modeltemplate;
%     models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.shared_segments;
%     models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.shared_segments';
%     this_model_name{end+1} = [cross_decode_label_pairs{i,1} ' to ' cross_decode_label_pairs{i,2} ' SS_Match'];
%     %Optional check - view matrix
%     %                     imagesc(models{end},'AlphaData',~isnan(models{end}))
%     %                     title(this_model_name{end})
%     %                     pause
% end


%Now attempt cross-condition shared segments RSA without cross decoding, recognising that the MisMatch cue
%was consistently 8 elements after/before the auditory word
basemodels.shared_segments(1:17:end) = NaN;
for i = 1:size(cross_decode_label_pairs,1)
    models{end+1} = modeltemplate;
    models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.shared_segments;
    models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.shared_segments';
    this_model_name{end+1} = [cross_decode_label_pairs{i,1} ' to ' cross_decode_label_pairs{i,2} ' SS_Match - no self'];
    %Optional check - view matrix
    %                     imagesc(models{end},'AlphaData',~isnan(models{end}))
    %                     title(this_model_name{end})
    %                     pause
end
basemodels.shared_segments(1:17:end) = 1;

% Now add combined conditions
cross_decode_label_pairs = {
    'Match Unclear', 'Mismatch Unclear';
    'Match Clear', 'Mismatch Unclear';
    'Match Unclear', 'Mismatch Clear';
    'Match Clear', 'Mismatch Clear';
    'Match Unclear', 'Match Clear';
    'Mismatch Unclear', 'Mismatch Clear';};

models{end+1} = modeltemplate;
this_model_name{end+1} = ['All spoken Cross-decode_Match'];
for i = 1:size(cross_decode_label_pairs,1)
    models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = Match_Cross_decode_base;
    models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = Match_Cross_decode_base';
end

models{end+1} = modeltemplate;
this_model_name{end+1} = ['All spoken SS_Match'];
for i = 1:size(cross_decode_label_pairs,1)
    models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.shared_segments;
    models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.shared_segments';
end
models{end+1} = modeltemplate;
this_model_name{end+1} = ['All spoken SS_Match - no self'];
basemodels.shared_segments(1:17:end) = NaN;
for i = 1:size(cross_decode_label_pairs,1)
    models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.shared_segments;
    models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.shared_segments';
end
basemodels.shared_segments(1:17:end) = 0;

cross_decode_label_pairs = {
    'Match Unclear', 'Written';
    'Match Clear', 'Written';
    'Mismatch Unclear', 'Written';
    'Mismatch Clear', 'Written'
    };
% models{end+1} = modeltemplate;
% this_model_name{end+1} = ['Spoken to Written Cross-decode_Match'];
% for i = 1:size(cross_decode_label_pairs,1)
%     models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = Match_Cross_decode_base;
%     models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = Match_Cross_decode_base';
% end
% models{end+1} = modeltemplate;
% this_model_name{end+1} = ['Spoken to Written SS_Match - no self'];
% basemodels.shared_segments(1:17:end) = NaN;
% for i = 1:size(cross_decode_label_pairs,1)
%     models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.shared_segments;
%     models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.shared_segments';
% end
% basemodels.shared_segments(1:17:end) = 0;

% Now look at Match-Mismatch written word cross decoding
cross_decode_label_pairs = {
    'Match Unclear', 'Mismatch Unclear';
    'Match Clear', 'Mismatch Unclear';
    'Match Unclear', 'Mismatch Clear';
    'Match Clear', 'Mismatch Clear';
    };

models{end+1} = modeltemplate;
this_model_name{end+1} = ['Match to Mismatch Shared Segments - no self'];
for i = 1:size(cross_decode_label_pairs,1)
    models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.shared_segments_cross_noself;
    models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.shared_segments_cross_noself';
end

models{end+1} = modeltemplate;
this_model_name{end+1} = ['Match to Mismatch SS_Match - no self'];
basemodels.shared_segments(1:17:end) = NaN;
for i = 1:size(cross_decode_label_pairs,1)
    models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.shared_segments;
    models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.shared_segments';
end
basemodels.shared_segments(1:17:end) = 0;

basemodels.shared_segments(1:17:end) = NaN;
basemodels.combined_SS = basemodels.shared_segments-basemodels.shared_segments_cross_noself;
basemodels.shared_segments(1:17:end) = 0;

basemodels.combined_SS = (basemodels.combined_SS +1)/2; %Scale zero to 1
% models{end+1} = modeltemplate;
% this_model_name{end+1} = ['Match to Mismatch combined_SS - no self - rescaled'];
% for i = 1:size(cross_decode_label_pairs,1)
%     models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.combined_SS;
%     models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.combined_SS';
% end

basemodels.only_cross = basemodels.combined_SS;
basemodels.only_cross(basemodels.shared_segments~=1) = NaN;
models{end+1} = modeltemplate;
this_model_name{end+1} = ['Match to Mismatch only cross'];
for i = 1:size(cross_decode_label_pairs,1)
    models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.only_cross;
    models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.only_cross';
end

basemodels.only_not_cross = basemodels.combined_SS;
basemodels.only_not_cross(basemodels.shared_segments_cross_noself~=1) = NaN;
models{end+1} = modeltemplate;
this_model_name{end+1} = ['Match to Mismatch only not cross'];
for i = 1:size(cross_decode_label_pairs,1)
    models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.only_not_cross;
    models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.only_not_cross';
end

% Now look at these three models Clear-Unclear
% then in every individual combination
cross_decode_label_pairs = {
    'Match Unclear', 'Match Clear';
    'Mismatch Unclear', 'Mismatch Clear';
    'Match Unclear', 'Mismatch Clear';
    'Match Clear', 'Mismatch Unclear';
    };
models{end+1} = modeltemplate;
this_model_name{end+1} = ['Clear to Unclear combined_SS - no self - rescaled'];
for i = 1:size(cross_decode_label_pairs,1)
    models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.combined_SS;
    models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.combined_SS';
end

models{end+1} = modeltemplate;
this_model_name{end+1} = ['Clear to Unclear only cross'];
for i = 1:size(cross_decode_label_pairs,1)
    models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.only_cross;
    models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.only_cross';
end

models{end+1} = modeltemplate;
this_model_name{end+1} = ['Clear to Unclear only not cross'];
for i = 1:size(cross_decode_label_pairs,1)
    models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.only_not_cross;
    models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.only_not_cross';
end

cross_decode_label_pairs = {
    'Match Unclear', 'Mismatch Unclear';
    'Match Clear', 'Mismatch Unclear';
    'Match Unclear', 'Mismatch Clear';
    'Match Clear', 'Mismatch Clear'
    'Match Unclear', 'Match Clear';
    'Mismatch Unclear', 'Mismatch Clear';
    };

for i = 1:size(cross_decode_label_pairs,1)
    models{end+1} = modeltemplate;
    models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.combined_SS;
    models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.combined_SS';
    this_model_name{end+1} = [cross_decode_label_pairs{i,1} ' to ' cross_decode_label_pairs{i,2} ' combined_SS - no self'];
    %Optional check - view matrix
    %             imagesc(models{end},'AlphaData',~isnan(models{end}))
    %             title(this_model_name{end})
    %             pause
end

for i = 1:size(cross_decode_label_pairs,1)
    models{end+1} = modeltemplate;
    models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.only_cross;
    models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.only_cross';
    this_model_name{end+1} = [cross_decode_label_pairs{i,1} ' to ' cross_decode_label_pairs{i,2} ' only cross'];
    %Optional check - view matrix
    %             imagesc(models{end},'AlphaData',~isnan(models{end}))
    %             title(this_model_name{end})
    %             pause
end

for i = 1:size(cross_decode_label_pairs,1)
    models{end+1} = modeltemplate;
    models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.only_not_cross;
    models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.only_not_cross';
    this_model_name{end+1} = [cross_decode_label_pairs{i,1} ' to ' cross_decode_label_pairs{i,2} ' only not cross'];
    %Optional check - view matrix
    %             imagesc(models{end},'AlphaData',~isnan(models{end}))
    %             title(this_model_name{end})
    %             pause
end

% % models = modelRDMs; close all
% %modelNames = fieldnames(models);
%
% %modelNames = {'ProbM' 'ProbMM' 'EntropyM' 'EntropyMM'};


V = spm_vol(fullfile(GLMDir,'mask.nii'));
mask = spm_read_vols(V);
mask_index = results.mask_index;
downsamped_V.mat = V.mat;
downsamped_V.mat(1:3,1:3)=downsamped_V.mat(1:3,1:3)*downsamp_ratio;

clear results % to free memory

for m=1:length(this_model_name) %Parallelising here impossible due to out of memory on serialisation unless data downsampled
    fprintf('\nComputing effect-map for model %s\n',this_model_name{m});
    if ~exist(fullfile(outputDir,['effect-map_' this_model_name{m} '.nii'])) || redo_maps == 1
        modelRDM = vectorizeRDMs(models{m})';
        effectMap = NaN(size(mask));
        for vx=1:numel(data)
            neuralRDM = vectorizeRDMs(data{vx})';
            if isempty(neuralRDM)
                continue
            end
            notempty = vx;
            if ~isempty(strfind(version,'pearson'))
                effectMap(mask_index(vx)) = fisherTransform(corr(modelRDM,neuralRDM,'type','Pearson','Rows','pairwise'));
            elseif ~isempty(strfind(version,'spearman'))
                effectMap(mask_index(vx)) = fisherTransform(corr(modelRDM,neuralRDM,'type','Spearman','Rows','pairwise'));
            elseif ~isempty(strfind(version,'average'))
                %effectMap(mask_index(vx)) = mean(neuralRDM(find(~isnan(modelRDM)),:),1);
                effectMap(mask_index(vx)) = mean(neuralRDM(find(modelRDM==1),:),1);
            end
            if ~mod(vx,100)
                disp(['Processing voxel ' num2str(vx) ' of ' num2str(numel(data))])
            end
        end
        dims = size(effectMap);
        downsamped_effectMap = effectMap(1:downsamp_ratio:dims(1),1:downsamp_ratio:dims(2),1:downsamp_ratio:dims(3));
        
        saveMRImage(downsamped_effectMap,fullfile(outputDir,['effect-map_' this_model_name{m} '.nii']),downsamped_V.mat);
    else
        disp('Already exists - moving on - set redo_maps to 1 if you want to re-make')
    end
end