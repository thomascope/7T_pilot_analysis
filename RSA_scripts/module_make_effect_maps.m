function module_make_effect_maps(GLMDir,downsamp_ratio)
%For taking already calculated crossnobis distances and doing RSA

redo_maps = 0; %If you want to calculate them again for some reason.
save_design_matrices = 0; %If you want to output the design matrices for visualisation

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

basemodels.vowels = zeros(16,16);
basemodels.vowels(1:17:end) = 1;
basemodels.vowels(2:68:end) = 1/3;
basemodels.vowels(3:68:end) = 1/3;
basemodels.vowels(4:68:end) = 1/3;
basemodels.vowels(17:68:end) = 1/3;
basemodels.vowels(19:68:end) = 1/3;
basemodels.vowels(20:68:end) = 1/3;
basemodels.vowels(33:68:end) = 1/3;
basemodels.vowels(34:68:end) = 1/3;
basemodels.vowels(36:68:end) = 1/3;
basemodels.vowels(49:68:end) = 1/3;
basemodels.vowels(50:68:end) = 1/3;
basemodels.vowels(51:68:end) = 1/3;
basemodels.vowels = 1-basemodels.vowels;

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
basemodels.shared_segments_mismatch = circshift(basemodels.shared_segments,[8 8]); %Consistent relationship for mismatch - 8 items later

basemodels.shared_segments_cross = circshift(basemodels.shared_segments,[8 0]); %I think this is correct, but need to 100% check the off-diagonals


basemodels.shared_segments_both = (basemodels.shared_segments_mismatch + basemodels.shared_segments)/2;

basemodels.shared_features = [ % Based on Helen Blank's feature model. Similar but not identical to shared segments.
    1.0000    0.7533    0.7533    0.5067    0.5067    0.3833    0.3833    0.2600    0.5067    0.3833    0.3833    0.2600    0.2600    0.3833    0.2600    0.3833
    0.7533    1.0000    0.5067    0.7533    0.2600    0.2600    0.1367    0.1367    0.2600    0.2600    0.1367    0.1367    0.1367    0.1367    0.1367    0.1367
    0.7533    0.5067    1.0000    0.7533    0.2600    0.1367    0.2600    0.1367    0.5067    0.3833    0.3833    0.2600    0.3833    0.5067    0.5067    0.6300
    0.5067    0.7533    0.7533    1.0000    0.0133    0.0133    0.0133    0.0133    0.2600    0.2600    0.1367    0.1367    0.2600    0.2600    0.3833    0.3833
    0.5067    0.2600    0.2600    0.0133    1.0000    0.8767    0.8767    0.7533    0.6300    0.5067    0.5067    0.3833   -0.1100    0.0133   -0.1100    0.0133
    0.3833    0.2600    0.1367    0.0133    0.8767    1.0000    0.7533    0.8767    0.5067    0.6300    0.3833    0.5067   -0.1100   -0.1100   -0.1100   -0.1100
    0.3833    0.1367    0.2600    0.0133    0.8767    0.7533    1.0000    0.8767    0.6300    0.5067    0.5067    0.3833   -0.1100    0.0133   -0.1100    0.0133
    0.2600    0.1367    0.1367    0.0133    0.7533    0.8767    0.8767    1.0000    0.5067    0.6300    0.3833    0.5067   -0.1100   -0.1100   -0.1100   -0.1100
    0.5067    0.2600    0.5067    0.2600    0.6300    0.5067    0.6300    0.5067    1.0000    0.8767    0.8767    0.7533    0.1367    0.2600    0.2600    0.3833
    0.3833    0.2600    0.3833    0.2600    0.5067    0.6300    0.5067    0.6300    0.8767    1.0000    0.7533    0.8767    0.1367    0.1367    0.2600    0.2600
    0.3833    0.1367    0.3833    0.1367    0.5067    0.3833    0.5067    0.3833    0.8767    0.7533    1.0000    0.8767    0.0133    0.1367    0.1367    0.2600
    0.2600    0.1367    0.2600    0.1367    0.3833    0.5067    0.3833    0.5067    0.7533    0.8767    0.8767    1.0000    0.0133    0.0133    0.1367    0.1367
    0.2600    0.1367    0.3833    0.2600   -0.1100   -0.1100   -0.1100   -0.1100    0.1367    0.1367    0.0133    0.0133    1.0000    0.8767    0.8767    0.7533
    0.3833    0.1367    0.5067    0.2600    0.0133   -0.1100    0.0133   -0.1100    0.2600    0.1367    0.1367    0.0133    0.8767    1.0000    0.7533    0.8767
    0.2600    0.1367    0.5067    0.3833   -0.1100   -0.1100   -0.1100   -0.1100    0.2600    0.2600    0.1367    0.1367    0.8767    0.7533    1.0000    0.8767
    0.3833    0.1367    0.6300    0.3833    0.0133   -0.1100    0.0133   -0.1100    0.3833    0.2600    0.2600    0.1367    0.7533    0.8767    0.8767    1.0000
    ];
basemodels.shared_features = 1-basemodels.shared_features;

basemodelNames = {'vowels','shared_segments','shared_segments_mismatch','shared_segments_both','shared_features'};


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

basemodelNames = {'vowels','shared_segments','shared_segments_mismatch','shared_features'};

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
    %Optional check - view matrix
    %                     imagesc(models{end},'AlphaData',~isnan(models{end}))
    %                     title(this_model_name{end},'Interpreter','none')
    %                     colorbar
    %                     pause
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

for i = 1:size(cross_decode_label_pairs,1)
    models{end+1} = modeltemplate;
    models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = MisMatch_Cross_decode_base;
    models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = MisMatch_Cross_decode_base';
    this_model_name{end+1} = [cross_decode_label_pairs{i,1} ' to ' cross_decode_label_pairs{i,2} ' Cross-decode'];
    %Optional check - view matrix
    %             imagesc(models{end},'AlphaData',~isnan(models{end}))
    %             title(this_model_name{end})
    %             pause
end

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

for i = 1:size(cross_decode_label_pairs,1)
    models{end+1} = modeltemplate;
    models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = Match_Cross_decode_base;
    models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = Match_Cross_decode_base';
    this_model_name{end+1} = [cross_decode_label_pairs{i,1} ' to ' cross_decode_label_pairs{i,2} ' Cross-decode_Match'];
    %Optional check - view matrix
    %             imagesc(models{end},'AlphaData',~isnan(models{end}))
    %             title(this_model_name{end})
    %             pause
end

%Now attempt cross-condition shared segments RSA without cross decoding, recognising that the MisMatch cue
%was consistently 8 elements after/before the auditory word
for i = 1:size(cross_decode_label_pairs,1)
    models{end+1} = modeltemplate;
    models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.shared_segments;
    models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.shared_segments';
    this_model_name{end+1} = [cross_decode_label_pairs{i,1} ' to ' cross_decode_label_pairs{i,2} ' SS_Match'];
    %Optional check - view matrix
    %                     imagesc(models{end},'AlphaData',~isnan(models{end}))
    %                     title(this_model_name{end})
    %                     pause
end


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
models{end+1} = modeltemplate;
this_model_name{end+1} = ['Spoken to Written Cross-decode_Match'];
for i = 1:size(cross_decode_label_pairs,1)
    models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = Match_Cross_decode_base;
    models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = Match_Cross_decode_base';
end
models{end+1} = modeltemplate;
this_model_name{end+1} = ['Spoken to Written SS_Match - no self'];
basemodels.shared_segments(1:17:end) = NaN;
for i = 1:size(cross_decode_label_pairs,1)
    models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.shared_segments;
    models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.shared_segments';
end
basemodels.shared_segments(1:17:end) = 0;

models{end+1} = modeltemplate;
this_model_name{end+1} = ['Spoken to Written Cross-decode_written'];
for i = 1:2
    models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = Match_Cross_decode_base;
    models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = Match_Cross_decode_base';
end
for i = 3:4
    models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = MisMatch_Cross_decode_base;
    models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = MisMatch_Cross_decode_base';
end
models{end+1} = modeltemplate;
this_model_name{end+1} = ['Spoken to Written Cross-decode_written-lowpe'];
for i = 2
    models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = Match_Cross_decode_base;
    models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = Match_Cross_decode_base';
end
for i = 3
    models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = MisMatch_Cross_decode_base;
    models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = MisMatch_Cross_decode_base';
end
models{end+1} = modeltemplate;
this_model_name{end+1} = ['Spoken to Written Cross-decode_written-highpe'];
for i = 1
    models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = Match_Cross_decode_base;
    models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = Match_Cross_decode_base';
end
for i = 4
    models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = MisMatch_Cross_decode_base;
    models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = MisMatch_Cross_decode_base';
end
models{end+1} = modeltemplate;
this_model_name{end+1} = ['Spoken to Written SS_written - no self'];
basemodels.shared_segments(1:17:end) = NaN;
basemodels.shared_segments_cross(9:17:end/2) = NaN;
basemodels.shared_segments_cross(end/2+1:17:end) = NaN;
for i = 1:2
    models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.shared_segments;
    models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.shared_segments';
end
for i = 3:4
    models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.shared_segments_cross;
    models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.shared_segments_cross';
end
basemodels.shared_segments_cross(9:17:end/2) = 0;
basemodels.shared_segments_cross(end/2+1:17:end) = 0;
basemodels.shared_segments(1:17:end) = 0;

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
models{end+1} = modeltemplate;
this_model_name{end+1} = ['Match to Mismatch combined_SS - no self - rescaled'];
for i = 1:size(cross_decode_label_pairs,1)
    models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.combined_SS;
    models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.combined_SS';
end

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

models{end+1} = modeltemplate;
this_model_name{end+1} = ['Match to Mismatch cross-decode written'];
for i = 1:size(cross_decode_label_pairs,1)
    models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = MisMatch_Cross_decode_base;
    models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = MisMatch_Cross_decode_base';
end

models{end+1} = modeltemplate;
this_model_name{end+1} = ['Match to Mismatch cross-decode spoken'];
for i = 1:size(cross_decode_label_pairs,1)
    models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = Match_Cross_decode_base;
    models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = Match_Cross_decode_base';
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

models{end+1} = [ % Acoustic cross-spectral coherence of the input words
    1	0.171993798459698	0.164625519206801	0.231267842236098	0.231350166174658	0.246711030825474	0.213208295988683	0.229833338423719	0.318378408079613	0.351497039664561	0.144094407612773	0.226996463762615	0.345749031558995	0.140564522824956	0.186537734790775	0.165667859725031	1	0.171993798459698	0.164625519206801	0.231267842236098	0.231350166174658	0.246711030825474	0.213208295988683	0.229833338423719	0.318378408079613	0.351497039664561	0.144094407612773	0.226996463762615	0.345749031558995	0.140564522824956	0.186537734790775	0.165667859725031	0.134735352555133	0.227214056458885	0.153367199073079	0.160072862647111	0.298078652654840	0.185048818439264	0.238202610510690	0.244657386305671	0.262800581807974	0.276177229421184	0.336676352839481	0.192605147348659	0.210038103397686	0.220212132020989	0.176620548805626	0.202230550437096	0.134735352555133	0.227214056458885	0.153367199073079	0.160072862647111	0.298078652654840	0.185048818439264	0.238202610510690	0.244657386305671	0.262800581807974	0.276177229421184	0.336676352839481	0.192605147348659	0.210038103397686	0.220212132020989	0.176620548805626	0.202230550437096
    0.171993798459698	1	0.125212088733127	0.196111972609968	0.276239308925127	0.241445932956491	0.267070777504442	0.256807104104544	0.327256498082194	0.256642250686977	0.320405241914020	0.166228755761452	0.124095521169373	0.125589307679648	0.222587021298084	0.161343613293296	0.171993798459698	1	0.125212088733127	0.196111972609968	0.276239308925127	0.241445932956491	0.267070777504442	0.256807104104544	0.327256498082194	0.256642250686977	0.320405241914020	0.166228755761452	0.124095521169373	0.125589307679648	0.222587021298084	0.161343613293296	0.173670491285841	0.288216912365522	0.158001480883866	0.331153070144339	0.194075382869973	0.258008929738211	0.256601897716442	0.290654300814771	0.319924976813934	0.236311880485008	0.193399058787936	0.174715944086205	0.251205002561756	0.155368943598672	0.125854532171030	0.172439908547135	0.173670491285841	0.288216912365522	0.158001480883866	0.331153070144339	0.194075382869973	0.258008929738211	0.256601897716442	0.290654300814771	0.319924976813934	0.236311880485008	0.193399058787936	0.174715944086205	0.251205002561756	0.155368943598672	0.125854532171030	0.172439908547135
    0.164625519206801	0.125212088733127	1	0.122074197351473	0.195143524820889	0.203307795388911	0.202601217198439	0.223895055783473	0.104686206593351	0.0897727602674931	0.234319363120630	0.261843417004457	0.249410865733282	0.232155430216566	0.211296806321270	0.181946481255141	0.164625519206801	0.125212088733127	1	0.122074197351473	0.195143524820889	0.203307795388911	0.202601217198439	0.223895055783473	0.104686206593351	0.0897727602674931	0.234319363120630	0.261843417004457	0.249410865733282	0.232155430216566	0.211296806321270	0.181946481255141	0.355438791892943	0.105011488754895	0.156555090607610	0.125346993220384	0.127546712821935	0.184111341471906	0.174741938753492	0.161328541375099	0.132588234744717	0.150571222907838	0.155217995450231	0.204002216638815	0.226516656921773	0.177749281176850	0.177912875604621	0.221049465091916	0.355438791892943	0.105011488754895	0.156555090607610	0.125346993220384	0.127546712821935	0.184111341471906	0.174741938753492	0.161328541375099	0.132588234744717	0.150571222907838	0.155217995450231	0.204002216638815	0.226516656921773	0.177749281176850	0.177912875604621	0.221049465091916
    0.231267842236098	0.196111972609968	0.122074197351473	1	0.278288418516294	0.259530520960673	0.222656482938236	0.309178647326195	0.218167575763155	0.127462626971511	0.233861564972150	0.189140223832126	0.252557763978371	0.189453445988760	0.199529730989422	0.297543855401901	0.231267842236098	0.196111972609968	0.122074197351473	1	0.278288418516294	0.259530520960673	0.222656482938236	0.309178647326195	0.218167575763155	0.127462626971511	0.233861564972150	0.189140223832126	0.252557763978371	0.189453445988760	0.199529730989422	0.297543855401901	0.138436483264120	0.179991998384704	0.196398065234584	0.281801135230052	0.342431532552267	0.197020666103735	0.284260054814097	0.306282197159308	0.267137279136037	0.180893477778343	0.362942511121388	0.189051164175225	0.245304855183173	0.335138555543597	0.157826223919072	0.205458969300860	0.138436483264120	0.179991998384704	0.196398065234584	0.281801135230052	0.342431532552267	0.197020666103735	0.284260054814097	0.306282197159308	0.267137279136037	0.180893477778343	0.362942511121388	0.189051164175225	0.245304855183173	0.335138555543597	0.157826223919072	0.205458969300860
    0.231350166174658	0.276239308925127	0.195143524820889	0.278288418516294	1	0.514404571410664	0.350631309245543	0.540904747159626	0.429512817819356	0.268421862728620	0.445841051867186	0.214840881926063	0.228870176601458	0.196661761399125	0.293908654529095	0.252775414080312	0.231350166174658	0.276239308925127	0.195143524820889	0.278288418516294	1	0.514404571410664	0.350631309245543	0.540904747159626	0.429512817819356	0.268421862728620	0.445841051867186	0.214840881926063	0.228870176601458	0.196661761399125	0.293908654529095	0.252775414080312	0.165966452020543	0.256167797179826	0.237889094103677	0.384173275699749	0.368798215585550	0.419237674218629	0.488725641082833	0.494811701913341	0.414268088297245	0.301894085772881	0.368203777790144	0.358403653653469	0.463927096573370	0.192235561370051	0.148883560819950	0.255283901846791	0.165966452020543	0.256167797179826	0.237889094103677	0.384173275699749	0.368798215585550	0.419237674218629	0.488725641082833	0.494811701913341	0.414268088297245	0.301894085772881	0.368203777790144	0.358403653653469	0.463927096573370	0.192235561370051	0.148883560819950	0.255283901846791
    0.246711030825474	0.241445932956491	0.203307795388911	0.259530520960673	0.514404571410664	1	0.391592493830467	0.544910193311925	0.394221616243283	0.315985379112176	0.407242678540514	0.238934742576366	0.285327317210685	0.221210508105597	0.283421560120358	0.248411578301046	0.246711030825474	0.241445932956491	0.203307795388911	0.259530520960673	0.514404571410664	1	0.391592493830467	0.544910193311925	0.394221616243283	0.315985379112176	0.407242678540514	0.238934742576366	0.285327317210685	0.221210508105597	0.283421560120358	0.248411578301046	0.160598464390537	0.246367399139255	0.247601458225063	0.299630847277878	0.353290735799105	0.329583131194398	0.406912056368750	0.422318930549578	0.330528228726873	0.343883123574280	0.322818832106714	0.380010451159630	0.375553378164681	0.158512710264046	0.178713286156517	0.206601766747013	0.160598464390537	0.246367399139255	0.247601458225063	0.299630847277878	0.353290735799105	0.329583131194398	0.406912056368750	0.422318930549578	0.330528228726873	0.343883123574280	0.322818832106714	0.380010451159630	0.375553378164681	0.158512710264046	0.178713286156517	0.206601766747013
    0.213208295988683	0.267070777504442	0.202601217198439	0.222656482938236	0.350631309245543	0.391592493830467	1	0.389363699618271	0.288066335326684	0.278354253430404	0.369391580664048	0.276767013488525	0.273201019595578	0.209021812469469	0.314025712029172	0.203477291051575	0.213208295988683	0.267070777504442	0.202601217198439	0.222656482938236	0.350631309245543	0.391592493830467	1	0.389363699618271	0.288066335326684	0.278354253430404	0.369391580664048	0.276767013488525	0.273201019595578	0.209021812469469	0.314025712029172	0.203477291051575	0.187873073925526	0.245990326829341	0.184067479775407	0.242997536756078	0.233877443204227	0.258738603216210	0.237779012053277	0.259812235085185	0.297618588058598	0.376420241926500	0.194651469874191	0.235946694846246	0.212040428446836	0.180430521800859	0.216952897766093	0.160941735562725	0.187873073925526	0.245990326829341	0.184067479775407	0.242997536756078	0.233877443204227	0.258738603216210	0.237779012053277	0.259812235085185	0.297618588058598	0.376420241926500	0.194651469874191	0.235946694846246	0.212040428446836	0.180430521800859	0.216952897766093	0.160941735562725
    0.229833338423719	0.256807104104544	0.223895055783473	0.309178647326195	0.540904747159626	0.544910193311925	0.389363699618271	1	0.380555376691428	0.256396883310910	0.483693150841733	0.297393710933921	0.305329865786777	0.265914244590708	0.343516731042126	0.309709719448287	0.229833338423719	0.256807104104544	0.223895055783473	0.309178647326195	0.540904747159626	0.544910193311925	0.389363699618271	1	0.380555376691428	0.256396883310910	0.483693150841733	0.297393710933921	0.305329865786777	0.265914244590708	0.343516731042126	0.309709719448287	0.161237062608672	0.242600708434698	0.282127225606543	0.382027794502961	0.392846364702732	0.428719565972320	0.489612201661955	0.515545813640042	0.355080589179632	0.304338042723125	0.391411058109650	0.440688572535383	0.474463527483964	0.194559122028008	0.196668230175339	0.248086166274221	0.161237062608672	0.242600708434698	0.282127225606543	0.382027794502961	0.392846364702732	0.428719565972320	0.489612201661955	0.515545813640042	0.355080589179632	0.304338042723125	0.391411058109650	0.440688572535383	0.474463527483964	0.194559122028008	0.196668230175339	0.248086166274221
    0.318378408079613	0.327256498082194	0.104686206593351	0.218167575763155	0.429512817819356	0.394221616243283	0.288066335326684	0.380555376691428	1	0.456389173427060	0.302695229942900	0.126358246947378	0.156110200461239	0.0994905099299154	0.196953822162161	0.151575379870957	0.318378408079613	0.327256498082194	0.104686206593351	0.218167575763155	0.429512817819356	0.394221616243283	0.288066335326684	0.380555376691428	1	0.456389173427060	0.302695229942900	0.126358246947378	0.156110200461239	0.0994905099299154	0.196953822162161	0.151575379870957	0.142171564257541	0.348498229969896	0.166330902174613	0.344051373683414	0.315111126481233	0.346697694484368	0.397214414895051	0.426482935033363	0.449983936871318	0.331400876549848	0.301921963123240	0.245949176015459	0.364937088905025	0.137260384502674	0.0984872355958472	0.198903649734215	0.142171564257541	0.348498229969896	0.166330902174613	0.344051373683414	0.315111126481233	0.346697694484368	0.397214414895051	0.426482935033363	0.449983936871318	0.331400876549848	0.301921963123240	0.245949176015459	0.364937088905025	0.137260384502674	0.0984872355958472	0.198903649734215
    0.351497039664561	0.256642250686977	0.0897727602674931	0.127462626971511	0.268421862728620	0.315985379112176	0.278354253430404	0.256396883310910	0.456389173427060	1	0.209840381981601	0.211967986974051	0.238591953291001	0.103234120703231	0.145408926765413	0.0891055998525839	0.351497039664561	0.256642250686977	0.0897727602674931	0.127462626971511	0.268421862728620	0.315985379112176	0.278354253430404	0.256396883310910	0.456389173427060	1	0.209840381981601	0.211967986974051	0.238591953291001	0.103234120703231	0.145408926765413	0.0891055998525839	0.144450005975650	0.328899815080986	0.172561189924276	0.179401236672933	0.237993322355195	0.244090267580614	0.225779144247316	0.259439168982143	0.334655469346739	0.472056460404923	0.197605829777781	0.262489421725674	0.204053926746631	0.0868603132242725	0.113088587100550	0.126335122236854	0.144450005975650	0.328899815080986	0.172561189924276	0.179401236672933	0.237993322355195	0.244090267580614	0.225779144247316	0.259439168982143	0.334655469346739	0.472056460404923	0.197605829777781	0.262489421725674	0.204053926746631	0.0868603132242725	0.113088587100550	0.126335122236854
    0.144094407612773	0.320405241914020	0.234319363120630	0.233861564972150	0.445841051867186	0.407242678540514	0.369391580664048	0.483693150841733	0.302695229942900	0.209840381981601	1	0.367014357855516	0.259621327924974	0.268137420167537	0.401965629092308	0.280603239781068	0.144094407612773	0.320405241914020	0.234319363120630	0.233861564972150	0.445841051867186	0.407242678540514	0.369391580664048	0.483693150841733	0.302695229942900	0.209840381981601	1	0.367014357855516	0.259621327924974	0.268137420167537	0.401965629092308	0.280603239781068	0.185045963637839	0.230621147338107	0.223620976155853	0.389218444585022	0.222100397821173	0.391549984585214	0.371597503966509	0.405186841441442	0.365728800397777	0.325338632802839	0.234284764849756	0.356423707451147	0.430567652576515	0.201202640961916	0.188398126858605	0.214103184190751	0.185045963637839	0.230621147338107	0.223620976155853	0.389218444585022	0.222100397821173	0.391549984585214	0.371597503966509	0.405186841441442	0.365728800397777	0.325338632802839	0.234284764849756	0.356423707451147	0.430567652576515	0.201202640961916	0.188398126858605	0.214103184190751
    0.226996463762615	0.166228755761452	0.261843417004457	0.189140223832126	0.214840881926063	0.238934742576366	0.276767013488525	0.297393710933921	0.126358246947378	0.211967986974051	0.367014357855516	1	0.380166616052408	0.301284435228734	0.355713134972479	0.249782437611580	0.226996463762615	0.166228755761452	0.261843417004457	0.189140223832126	0.214840881926063	0.238934742576366	0.276767013488525	0.297393710933921	0.126358246947378	0.211967986974051	0.367014357855516	1	0.380166616052408	0.301284435228734	0.355713134972479	0.249782437611580	0.200141851055337	0.140964588095181	0.222178343670103	0.214474683217675	0.164399475282110	0.250255414156282	0.211219450709919	0.236420964685782	0.204757998601009	0.287411091053269	0.201923713117651	0.337905005742051	0.281344602087222	0.241598343196672	0.225351602112520	0.181966122544000	0.200141851055337	0.140964588095181	0.222178343670103	0.214474683217675	0.164399475282110	0.250255414156282	0.211219450709919	0.236420964685782	0.204757998601009	0.287411091053269	0.201923713117651	0.337905005742051	0.281344602087222	0.241598343196672	0.225351602112520	0.181966122544000
    0.345749031558995	0.124095521169373	0.249410865733282	0.252557763978371	0.228870176601458	0.285327317210685	0.273201019595578	0.305329865786777	0.156110200461239	0.238591953291001	0.259621327924974	0.380166616052408	1	0.342408869116889	0.252558136830167	0.221694819653684	0.345749031558995	0.124095521169373	0.249410865733282	0.252557763978371	0.228870176601458	0.285327317210685	0.273201019595578	0.305329865786777	0.156110200461239	0.238591953291001	0.259621327924974	0.380166616052408	1	0.342408869116889	0.252558136830167	0.221694819653684	0.160464678105875	0.154103221330077	0.168945275101099	0.140221432295813	0.329975873220386	0.158490458520637	0.191569207123187	0.206721277320307	0.189705473665356	0.305240889445556	0.349823875877798	0.251851533104586	0.189959042834344	0.171994402676464	0.246868236409744	0.223233013220737	0.160464678105875	0.154103221330077	0.168945275101099	0.140221432295813	0.329975873220386	0.158490458520637	0.191569207123187	0.206721277320307	0.189705473665356	0.305240889445556	0.349823875877798	0.251851533104586	0.189959042834344	0.171994402676464	0.246868236409744	0.223233013220737
    0.140564522824956	0.125589307679648	0.232155430216566	0.189453445988760	0.196661761399125	0.221210508105597	0.209021812469469	0.265914244590708	0.0994905099299154	0.103234120703231	0.268137420167537	0.301284435228734	0.342408869116889	1	0.282478416733625	0.313647391314683	0.140564522824956	0.125589307679648	0.232155430216566	0.189453445988760	0.196661761399125	0.221210508105597	0.209021812469469	0.265914244590708	0.0994905099299154	0.103234120703231	0.268137420167537	0.301284435228734	0.342408869116889	1	0.282478416733625	0.313647391314683	0.266779427191495	0.0869095597410673	0.271251740423712	0.181456465610374	0.184490367591133	0.194274117433449	0.199568957473710	0.215356962511248	0.146755406974320	0.184024888821981	0.203844847782750	0.251985640799158	0.250932620433989	0.306173680748464	0.244212520546245	0.194598448235030	0.266779427191495	0.0869095597410673	0.271251740423712	0.181456465610374	0.184490367591133	0.194274117433449	0.199568957473710	0.215356962511248	0.146755406974320	0.184024888821981	0.203844847782750	0.251985640799158	0.250932620433989	0.306173680748464	0.244212520546245	0.194598448235030
    0.186537734790775	0.222587021298084	0.211296806321270	0.199529730989422	0.293908654529095	0.283421560120358	0.314025712029172	0.343516731042126	0.196953822162161	0.145408926765413	0.401965629092308	0.355713134972479	0.252558136830167	0.282478416733625	1	0.254035307869095	0.186537734790775	0.222587021298084	0.211296806321270	0.199529730989422	0.293908654529095	0.283421560120358	0.314025712029172	0.343516731042126	0.196953822162161	0.145408926765413	0.401965629092308	0.355713134972479	0.252558136830167	0.282478416733625	1	0.254035307869095	0.184790520389430	0.195748174989036	0.143985932779031	0.312044394341670	0.186278817895424	0.277696610296234	0.263292195639382	0.295776267565491	0.263482123875376	0.246628517932146	0.201473876493201	0.276081114681235	0.290965230513544	0.288494413544832	0.324120487388837	0.214320506259475	0.184790520389430	0.195748174989036	0.143985932779031	0.312044394341670	0.186278817895424	0.277696610296234	0.263292195639382	0.295776267565491	0.263482123875376	0.246628517932146	0.201473876493201	0.276081114681235	0.290965230513544	0.288494413544832	0.324120487388837	0.214320506259475
    0.165667859725031	0.161343613293296	0.181946481255141	0.297543855401901	0.252775414080312	0.248411578301046	0.203477291051575	0.309709719448287	0.151575379870957	0.0891055998525839	0.280603239781068	0.249782437611580	0.221694819653684	0.313647391314683	0.254035307869095	1	0.165667859725031	0.161343613293296	0.181946481255141	0.297543855401901	0.252775414080312	0.248411578301046	0.203477291051575	0.309709719448287	0.151575379870957	0.0891055998525839	0.280603239781068	0.249782437611580	0.221694819653684	0.313647391314683	0.254035307869095	1	0.188908568954180	0.130666084046134	0.231153257653851	0.215809119317150	0.263484082424131	0.237582011612920	0.286781069416484	0.309485703130578	0.205849622984209	0.139276537033501	0.284903071676509	0.268899941359785	0.304286170642738	0.350843129007248	0.162236432601847	0.336763198993654	0.188908568954180	0.130666084046134	0.231153257653851	0.215809119317150	0.263484082424131	0.237582011612920	0.286781069416484	0.309485703130578	0.205849622984209	0.139276537033501	0.284903071676509	0.268899941359785	0.304286170642738	0.350843129007248	0.162236432601847	0.336763198993654
    1	0.171993798459698	0.164625519206801	0.231267842236098	0.231350166174658	0.246711030825474	0.213208295988683	0.229833338423719	0.318378408079613	0.351497039664561	0.144094407612773	0.226996463762615	0.345749031558995	0.140564522824956	0.186537734790775	0.165667859725031	1	0.171993798459698	0.164625519206801	0.231267842236098	0.231350166174658	0.246711030825474	0.213208295988683	0.229833338423719	0.318378408079613	0.351497039664561	0.144094407612773	0.226996463762615	0.345749031558995	0.140564522824956	0.186537734790775	0.165667859725031	0.134735352555133	0.227214056458885	0.153367199073079	0.160072862647111	0.298078652654840	0.185048818439264	0.238202610510690	0.244657386305671	0.262800581807974	0.276177229421184	0.336676352839481	0.192605147348659	0.210038103397686	0.220212132020989	0.176620548805626	0.202230550437096	0.134735352555133	0.227214056458885	0.153367199073079	0.160072862647111	0.298078652654840	0.185048818439264	0.238202610510690	0.244657386305671	0.262800581807974	0.276177229421184	0.336676352839481	0.192605147348659	0.210038103397686	0.220212132020989	0.176620548805626	0.202230550437096
    0.171993798459698	1	0.125212088733127	0.196111972609968	0.276239308925127	0.241445932956491	0.267070777504442	0.256807104104544	0.327256498082194	0.256642250686977	0.320405241914020	0.166228755761452	0.124095521169373	0.125589307679648	0.222587021298084	0.161343613293296	0.171993798459698	1	0.125212088733127	0.196111972609968	0.276239308925127	0.241445932956491	0.267070777504442	0.256807104104544	0.327256498082194	0.256642250686977	0.320405241914020	0.166228755761452	0.124095521169373	0.125589307679648	0.222587021298084	0.161343613293296	0.173670491285841	0.288216912365522	0.158001480883866	0.331153070144339	0.194075382869973	0.258008929738211	0.256601897716442	0.290654300814771	0.319924976813934	0.236311880485008	0.193399058787936	0.174715944086205	0.251205002561756	0.155368943598672	0.125854532171030	0.172439908547135	0.173670491285841	0.288216912365522	0.158001480883866	0.331153070144339	0.194075382869973	0.258008929738211	0.256601897716442	0.290654300814771	0.319924976813934	0.236311880485008	0.193399058787936	0.174715944086205	0.251205002561756	0.155368943598672	0.125854532171030	0.172439908547135
    0.164625519206801	0.125212088733127	1	0.122074197351473	0.195143524820889	0.203307795388911	0.202601217198439	0.223895055783473	0.104686206593351	0.0897727602674931	0.234319363120630	0.261843417004457	0.249410865733282	0.232155430216566	0.211296806321270	0.181946481255141	0.164625519206801	0.125212088733127	1	0.122074197351473	0.195143524820889	0.203307795388911	0.202601217198439	0.223895055783473	0.104686206593351	0.0897727602674931	0.234319363120630	0.261843417004457	0.249410865733282	0.232155430216566	0.211296806321270	0.181946481255141	0.355438791892943	0.105011488754895	0.156555090607610	0.125346993220384	0.127546712821935	0.184111341471906	0.174741938753492	0.161328541375099	0.132588234744717	0.150571222907838	0.155217995450231	0.204002216638815	0.226516656921773	0.177749281176850	0.177912875604621	0.221049465091916	0.355438791892943	0.105011488754895	0.156555090607610	0.125346993220384	0.127546712821935	0.184111341471906	0.174741938753492	0.161328541375099	0.132588234744717	0.150571222907838	0.155217995450231	0.204002216638815	0.226516656921773	0.177749281176850	0.177912875604621	0.221049465091916
    0.231267842236098	0.196111972609968	0.122074197351473	1	0.278288418516294	0.259530520960673	0.222656482938236	0.309178647326195	0.218167575763155	0.127462626971511	0.233861564972150	0.189140223832126	0.252557763978371	0.189453445988760	0.199529730989422	0.297543855401901	0.231267842236098	0.196111972609968	0.122074197351473	1	0.278288418516294	0.259530520960673	0.222656482938236	0.309178647326195	0.218167575763155	0.127462626971511	0.233861564972150	0.189140223832126	0.252557763978371	0.189453445988760	0.199529730989422	0.297543855401901	0.138436483264120	0.179991998384704	0.196398065234584	0.281801135230052	0.342431532552267	0.197020666103735	0.284260054814097	0.306282197159308	0.267137279136037	0.180893477778343	0.362942511121388	0.189051164175225	0.245304855183173	0.335138555543597	0.157826223919072	0.205458969300860	0.138436483264120	0.179991998384704	0.196398065234584	0.281801135230052	0.342431532552267	0.197020666103735	0.284260054814097	0.306282197159308	0.267137279136037	0.180893477778343	0.362942511121388	0.189051164175225	0.245304855183173	0.335138555543597	0.157826223919072	0.205458969300860
    0.231350166174658	0.276239308925127	0.195143524820889	0.278288418516294	1	0.514404571410664	0.350631309245543	0.540904747159626	0.429512817819356	0.268421862728620	0.445841051867186	0.214840881926063	0.228870176601458	0.196661761399125	0.293908654529095	0.252775414080312	0.231350166174658	0.276239308925127	0.195143524820889	0.278288418516294	1	0.514404571410664	0.350631309245543	0.540904747159626	0.429512817819356	0.268421862728620	0.445841051867186	0.214840881926063	0.228870176601458	0.196661761399125	0.293908654529095	0.252775414080312	0.165966452020543	0.256167797179826	0.237889094103677	0.384173275699749	0.368798215585550	0.419237674218629	0.488725641082833	0.494811701913341	0.414268088297245	0.301894085772881	0.368203777790144	0.358403653653469	0.463927096573370	0.192235561370051	0.148883560819950	0.255283901846791	0.165966452020543	0.256167797179826	0.237889094103677	0.384173275699749	0.368798215585550	0.419237674218629	0.488725641082833	0.494811701913341	0.414268088297245	0.301894085772881	0.368203777790144	0.358403653653469	0.463927096573370	0.192235561370051	0.148883560819950	0.255283901846791
    0.246711030825474	0.241445932956491	0.203307795388911	0.259530520960673	0.514404571410664	1	0.391592493830467	0.544910193311925	0.394221616243283	0.315985379112176	0.407242678540514	0.238934742576366	0.285327317210685	0.221210508105597	0.283421560120358	0.248411578301046	0.246711030825474	0.241445932956491	0.203307795388911	0.259530520960673	0.514404571410664	1	0.391592493830467	0.544910193311925	0.394221616243283	0.315985379112176	0.407242678540514	0.238934742576366	0.285327317210685	0.221210508105597	0.283421560120358	0.248411578301046	0.160598464390537	0.246367399139255	0.247601458225063	0.299630847277878	0.353290735799105	0.329583131194398	0.406912056368750	0.422318930549578	0.330528228726873	0.343883123574280	0.322818832106714	0.380010451159630	0.375553378164681	0.158512710264046	0.178713286156517	0.206601766747013	0.160598464390537	0.246367399139255	0.247601458225063	0.299630847277878	0.353290735799105	0.329583131194398	0.406912056368750	0.422318930549578	0.330528228726873	0.343883123574280	0.322818832106714	0.380010451159630	0.375553378164681	0.158512710264046	0.178713286156517	0.206601766747013
    0.213208295988683	0.267070777504442	0.202601217198439	0.222656482938236	0.350631309245543	0.391592493830467	1	0.389363699618271	0.288066335326684	0.278354253430404	0.369391580664048	0.276767013488525	0.273201019595578	0.209021812469469	0.314025712029172	0.203477291051575	0.213208295988683	0.267070777504442	0.202601217198439	0.222656482938236	0.350631309245543	0.391592493830467	1	0.389363699618271	0.288066335326684	0.278354253430404	0.369391580664048	0.276767013488525	0.273201019595578	0.209021812469469	0.314025712029172	0.203477291051575	0.187873073925526	0.245990326829341	0.184067479775407	0.242997536756078	0.233877443204227	0.258738603216210	0.237779012053277	0.259812235085185	0.297618588058598	0.376420241926500	0.194651469874191	0.235946694846246	0.212040428446836	0.180430521800859	0.216952897766093	0.160941735562725	0.187873073925526	0.245990326829341	0.184067479775407	0.242997536756078	0.233877443204227	0.258738603216210	0.237779012053277	0.259812235085185	0.297618588058598	0.376420241926500	0.194651469874191	0.235946694846246	0.212040428446836	0.180430521800859	0.216952897766093	0.160941735562725
    0.229833338423719	0.256807104104544	0.223895055783473	0.309178647326195	0.540904747159626	0.544910193311925	0.389363699618271	1	0.380555376691428	0.256396883310910	0.483693150841733	0.297393710933921	0.305329865786777	0.265914244590708	0.343516731042126	0.309709719448287	0.229833338423719	0.256807104104544	0.223895055783473	0.309178647326195	0.540904747159626	0.544910193311925	0.389363699618271	1	0.380555376691428	0.256396883310910	0.483693150841733	0.297393710933921	0.305329865786777	0.265914244590708	0.343516731042126	0.309709719448287	0.161237062608672	0.242600708434698	0.282127225606543	0.382027794502961	0.392846364702732	0.428719565972320	0.489612201661955	0.515545813640042	0.355080589179632	0.304338042723125	0.391411058109650	0.440688572535383	0.474463527483964	0.194559122028008	0.196668230175339	0.248086166274221	0.161237062608672	0.242600708434698	0.282127225606543	0.382027794502961	0.392846364702732	0.428719565972320	0.489612201661955	0.515545813640042	0.355080589179632	0.304338042723125	0.391411058109650	0.440688572535383	0.474463527483964	0.194559122028008	0.196668230175339	0.248086166274221
    0.318378408079613	0.327256498082194	0.104686206593351	0.218167575763155	0.429512817819356	0.394221616243283	0.288066335326684	0.380555376691428	1	0.456389173427060	0.302695229942900	0.126358246947378	0.156110200461239	0.0994905099299154	0.196953822162161	0.151575379870957	0.318378408079613	0.327256498082194	0.104686206593351	0.218167575763155	0.429512817819356	0.394221616243283	0.288066335326684	0.380555376691428	1	0.456389173427060	0.302695229942900	0.126358246947378	0.156110200461239	0.0994905099299154	0.196953822162161	0.151575379870957	0.142171564257541	0.348498229969896	0.166330902174613	0.344051373683414	0.315111126481233	0.346697694484368	0.397214414895051	0.426482935033363	0.449983936871318	0.331400876549848	0.301921963123240	0.245949176015459	0.364937088905025	0.137260384502674	0.0984872355958472	0.198903649734215	0.142171564257541	0.348498229969896	0.166330902174613	0.344051373683414	0.315111126481233	0.346697694484368	0.397214414895051	0.426482935033363	0.449983936871318	0.331400876549848	0.301921963123240	0.245949176015459	0.364937088905025	0.137260384502674	0.0984872355958472	0.198903649734215
    0.351497039664561	0.256642250686977	0.0897727602674931	0.127462626971511	0.268421862728620	0.315985379112176	0.278354253430404	0.256396883310910	0.456389173427060	1	0.209840381981601	0.211967986974051	0.238591953291001	0.103234120703231	0.145408926765413	0.0891055998525839	0.351497039664561	0.256642250686977	0.0897727602674931	0.127462626971511	0.268421862728620	0.315985379112176	0.278354253430404	0.256396883310910	0.456389173427060	1	0.209840381981601	0.211967986974051	0.238591953291001	0.103234120703231	0.145408926765413	0.0891055998525839	0.144450005975650	0.328899815080986	0.172561189924276	0.179401236672933	0.237993322355195	0.244090267580614	0.225779144247316	0.259439168982143	0.334655469346739	0.472056460404923	0.197605829777781	0.262489421725674	0.204053926746631	0.0868603132242725	0.113088587100550	0.126335122236854	0.144450005975650	0.328899815080986	0.172561189924276	0.179401236672933	0.237993322355195	0.244090267580614	0.225779144247316	0.259439168982143	0.334655469346739	0.472056460404923	0.197605829777781	0.262489421725674	0.204053926746631	0.0868603132242725	0.113088587100550	0.126335122236854
    0.144094407612773	0.320405241914020	0.234319363120630	0.233861564972150	0.445841051867186	0.407242678540514	0.369391580664048	0.483693150841733	0.302695229942900	0.209840381981601	1	0.367014357855516	0.259621327924974	0.268137420167537	0.401965629092308	0.280603239781068	0.144094407612773	0.320405241914020	0.234319363120630	0.233861564972150	0.445841051867186	0.407242678540514	0.369391580664048	0.483693150841733	0.302695229942900	0.209840381981601	1	0.367014357855516	0.259621327924974	0.268137420167537	0.401965629092308	0.280603239781068	0.185045963637839	0.230621147338107	0.223620976155853	0.389218444585022	0.222100397821173	0.391549984585214	0.371597503966509	0.405186841441442	0.365728800397777	0.325338632802839	0.234284764849756	0.356423707451147	0.430567652576515	0.201202640961916	0.188398126858605	0.214103184190751	0.185045963637839	0.230621147338107	0.223620976155853	0.389218444585022	0.222100397821173	0.391549984585214	0.371597503966509	0.405186841441442	0.365728800397777	0.325338632802839	0.234284764849756	0.356423707451147	0.430567652576515	0.201202640961916	0.188398126858605	0.214103184190751
    0.226996463762615	0.166228755761452	0.261843417004457	0.189140223832126	0.214840881926063	0.238934742576366	0.276767013488525	0.297393710933921	0.126358246947378	0.211967986974051	0.367014357855516	1	0.380166616052408	0.301284435228734	0.355713134972479	0.249782437611580	0.226996463762615	0.166228755761452	0.261843417004457	0.189140223832126	0.214840881926063	0.238934742576366	0.276767013488525	0.297393710933921	0.126358246947378	0.211967986974051	0.367014357855516	1	0.380166616052408	0.301284435228734	0.355713134972479	0.249782437611580	0.200141851055337	0.140964588095181	0.222178343670103	0.214474683217675	0.164399475282110	0.250255414156282	0.211219450709919	0.236420964685782	0.204757998601009	0.287411091053269	0.201923713117651	0.337905005742051	0.281344602087222	0.241598343196672	0.225351602112520	0.181966122544000	0.200141851055337	0.140964588095181	0.222178343670103	0.214474683217675	0.164399475282110	0.250255414156282	0.211219450709919	0.236420964685782	0.204757998601009	0.287411091053269	0.201923713117651	0.337905005742051	0.281344602087222	0.241598343196672	0.225351602112520	0.181966122544000
    0.345749031558995	0.124095521169373	0.249410865733282	0.252557763978371	0.228870176601458	0.285327317210685	0.273201019595578	0.305329865786777	0.156110200461239	0.238591953291001	0.259621327924974	0.380166616052408	1	0.342408869116889	0.252558136830167	0.221694819653684	0.345749031558995	0.124095521169373	0.249410865733282	0.252557763978371	0.228870176601458	0.285327317210685	0.273201019595578	0.305329865786777	0.156110200461239	0.238591953291001	0.259621327924974	0.380166616052408	1	0.342408869116889	0.252558136830167	0.221694819653684	0.160464678105875	0.154103221330077	0.168945275101099	0.140221432295813	0.329975873220386	0.158490458520637	0.191569207123187	0.206721277320307	0.189705473665356	0.305240889445556	0.349823875877798	0.251851533104586	0.189959042834344	0.171994402676464	0.246868236409744	0.223233013220737	0.160464678105875	0.154103221330077	0.168945275101099	0.140221432295813	0.329975873220386	0.158490458520637	0.191569207123187	0.206721277320307	0.189705473665356	0.305240889445556	0.349823875877798	0.251851533104586	0.189959042834344	0.171994402676464	0.246868236409744	0.223233013220737
    0.140564522824956	0.125589307679648	0.232155430216566	0.189453445988760	0.196661761399125	0.221210508105597	0.209021812469469	0.265914244590708	0.0994905099299154	0.103234120703231	0.268137420167537	0.301284435228734	0.342408869116889	1	0.282478416733625	0.313647391314683	0.140564522824956	0.125589307679648	0.232155430216566	0.189453445988760	0.196661761399125	0.221210508105597	0.209021812469469	0.265914244590708	0.0994905099299154	0.103234120703231	0.268137420167537	0.301284435228734	0.342408869116889	1	0.282478416733625	0.313647391314683	0.266779427191495	0.0869095597410673	0.271251740423712	0.181456465610374	0.184490367591133	0.194274117433449	0.199568957473710	0.215356962511248	0.146755406974320	0.184024888821981	0.203844847782750	0.251985640799158	0.250932620433989	0.306173680748464	0.244212520546245	0.194598448235030	0.266779427191495	0.0869095597410673	0.271251740423712	0.181456465610374	0.184490367591133	0.194274117433449	0.199568957473710	0.215356962511248	0.146755406974320	0.184024888821981	0.203844847782750	0.251985640799158	0.250932620433989	0.306173680748464	0.244212520546245	0.194598448235030
    0.186537734790775	0.222587021298084	0.211296806321270	0.199529730989422	0.293908654529095	0.283421560120358	0.314025712029172	0.343516731042126	0.196953822162161	0.145408926765413	0.401965629092308	0.355713134972479	0.252558136830167	0.282478416733625	1	0.254035307869095	0.186537734790775	0.222587021298084	0.211296806321270	0.199529730989422	0.293908654529095	0.283421560120358	0.314025712029172	0.343516731042126	0.196953822162161	0.145408926765413	0.401965629092308	0.355713134972479	0.252558136830167	0.282478416733625	1	0.254035307869095	0.184790520389430	0.195748174989036	0.143985932779031	0.312044394341670	0.186278817895424	0.277696610296234	0.263292195639382	0.295776267565491	0.263482123875376	0.246628517932146	0.201473876493201	0.276081114681235	0.290965230513544	0.288494413544832	0.324120487388837	0.214320506259475	0.184790520389430	0.195748174989036	0.143985932779031	0.312044394341670	0.186278817895424	0.277696610296234	0.263292195639382	0.295776267565491	0.263482123875376	0.246628517932146	0.201473876493201	0.276081114681235	0.290965230513544	0.288494413544832	0.324120487388837	0.214320506259475
    0.165667859725031	0.161343613293296	0.181946481255141	0.297543855401901	0.252775414080312	0.248411578301046	0.203477291051575	0.309709719448287	0.151575379870957	0.0891055998525839	0.280603239781068	0.249782437611580	0.221694819653684	0.313647391314683	0.254035307869095	1	0.165667859725031	0.161343613293296	0.181946481255141	0.297543855401901	0.252775414080312	0.248411578301046	0.203477291051575	0.309709719448287	0.151575379870957	0.0891055998525839	0.280603239781068	0.249782437611580	0.221694819653684	0.313647391314683	0.254035307869095	1	0.188908568954180	0.130666084046134	0.231153257653851	0.215809119317150	0.263484082424131	0.237582011612920	0.286781069416484	0.309485703130578	0.205849622984209	0.139276537033501	0.284903071676509	0.268899941359785	0.304286170642738	0.350843129007248	0.162236432601847	0.336763198993654	0.188908568954180	0.130666084046134	0.231153257653851	0.215809119317150	0.263484082424131	0.237582011612920	0.286781069416484	0.309485703130578	0.205849622984209	0.139276537033501	0.284903071676509	0.268899941359785	0.304286170642738	0.350843129007248	0.162236432601847	0.336763198993654
    0.134735352555133	0.173670491285841	0.355438791892943	0.138436483264120	0.165966452020543	0.160598464390537	0.187873073925526	0.161237062608672	0.142171564257541	0.144450005975650	0.185045963637839	0.200141851055337	0.160464678105875	0.266779427191495	0.184790520389430	0.188908568954180	0.134735352555133	0.173670491285841	0.355438791892943	0.138436483264120	0.165966452020543	0.160598464390537	0.187873073925526	0.161237062608672	0.142171564257541	0.144450005975650	0.185045963637839	0.200141851055337	0.160464678105875	0.266779427191495	0.184790520389430	0.188908568954180	1	0.119563571207939	0.134517649418353	0.145645342006557	0.124878680641713	0.212129259364463	0.168453869184571	0.165063486104721	0.168622540386407	0.156279017501782	0.0931167304667288	0.205201036957932	0.181695384085607	0.331381888723394	0.0694982205575988	0.0937803258353237	1	0.119563571207939	0.134517649418353	0.145645342006557	0.124878680641713	0.212129259364463	0.168453869184571	0.165063486104721	0.168622540386407	0.156279017501782	0.0931167304667288	0.205201036957932	0.181695384085607	0.331381888723394	0.0694982205575988	0.0937803258353237
    0.227214056458885	0.288216912365522	0.105011488754895	0.179991998384704	0.256167797179826	0.246367399139255	0.245990326829341	0.242600708434698	0.348498229969896	0.328899815080986	0.230621147338107	0.140964588095181	0.154103221330077	0.0869095597410673	0.195748174989036	0.130666084046134	0.227214056458885	0.288216912365522	0.105011488754895	0.179991998384704	0.256167797179826	0.246367399139255	0.245990326829341	0.242600708434698	0.348498229969896	0.328899815080986	0.230621147338107	0.140964588095181	0.154103221330077	0.0869095597410673	0.195748174989036	0.130666084046134	0.119563571207939	1	0.130300072561244	0.311823212613157	0.225668991451732	0.222888736595654	0.250347801641747	0.278881204611943	0.345335936696517	0.264580041902274	0.244921364491539	0.159471642239823	0.223885938932743	0.140128521425098	0.109383174686410	0.185481661708238	0.119563571207939	1	0.130300072561244	0.311823212613157	0.225668991451732	0.222888736595654	0.250347801641747	0.278881204611943	0.345335936696517	0.264580041902274	0.244921364491539	0.159471642239823	0.223885938932743	0.140128521425098	0.109383174686410	0.185481661708238
    0.153367199073079	0.158001480883866	0.156555090607610	0.196398065234584	0.237889094103677	0.247601458225063	0.184067479775407	0.282127225606543	0.166330902174613	0.172561189924276	0.223620976155853	0.222178343670103	0.168945275101099	0.271251740423712	0.143985932779031	0.231153257653851	0.153367199073079	0.158001480883866	0.156555090607610	0.196398065234584	0.237889094103677	0.247601458225063	0.184067479775407	0.282127225606543	0.166330902174613	0.172561189924276	0.223620976155853	0.222178343670103	0.168945275101099	0.271251740423712	0.143985932779031	0.231153257653851	0.134517649418353	0.130300072561244	1	0.293871915677903	0.252438920240734	0.318898807571854	0.323727296097083	0.349661378890827	0.0952671517244940	0.0871863300749362	0.267043889089147	0.440808634948773	0.365709285771055	0.317878714930271	0.107976903089386	0.147715667980565	0.134517649418353	0.130300072561244	1	0.293871915677903	0.252438920240734	0.318898807571854	0.323727296097083	0.349661378890827	0.0952671517244940	0.0871863300749362	0.267043889089147	0.440808634948773	0.365709285771055	0.317878714930271	0.107976903089386	0.147715667980565
    0.160072862647111	0.331153070144339	0.125346993220384	0.281801135230052	0.384173275699749	0.299630847277878	0.242997536756078	0.382027794502961	0.344051373683414	0.179401236672933	0.389218444585022	0.214474683217675	0.140221432295813	0.181456465610374	0.312044394341670	0.215809119317150	0.160072862647111	0.331153070144339	0.125346993220384	0.281801135230052	0.384173275699749	0.299630847277878	0.242997536756078	0.382027794502961	0.344051373683414	0.179401236672933	0.389218444585022	0.214474683217675	0.140221432295813	0.181456465610374	0.312044394341670	0.215809119317150	0.145645342006557	0.311823212613157	0.293871915677903	1	0.364173992674786	0.503458147248227	0.578195376186041	0.645342396741406	0.464914860750329	0.151752731806341	0.402303235121321	0.371449199942161	0.583114246861958	0.226871281996533	0.199785846755072	0.175053651101674	0.145645342006557	0.311823212613157	0.293871915677903	1	0.364173992674786	0.503458147248227	0.578195376186041	0.645342396741406	0.464914860750329	0.151752731806341	0.402303235121321	0.371449199942161	0.583114246861958	0.226871281996533	0.199785846755072	0.175053651101674
    0.298078652654840	0.194075382869973	0.127546712821935	0.342431532552267	0.368798215585550	0.353290735799105	0.233877443204227	0.392846364702732	0.315111126481233	0.237993322355195	0.222100397821173	0.164399475282110	0.329975873220386	0.184490367591133	0.186278817895424	0.263484082424131	0.298078652654840	0.194075382869973	0.127546712821935	0.342431532552267	0.368798215585550	0.353290735799105	0.233877443204227	0.392846364702732	0.315111126481233	0.237993322355195	0.222100397821173	0.164399475282110	0.329975873220386	0.184490367591133	0.186278817895424	0.263484082424131	0.124878680641713	0.225668991451732	0.252438920240734	0.364173992674786	1	0.367316113266388	0.550454917488211	0.557751258890851	0.359265234826973	0.169866450660821	0.620217762992566	0.387660620005432	0.384377072928048	0.260248742258431	0.228935748587690	0.273146684676011	0.124878680641713	0.225668991451732	0.252438920240734	0.364173992674786	1	0.367316113266388	0.550454917488211	0.557751258890851	0.359265234826973	0.169866450660821	0.620217762992566	0.387660620005432	0.384377072928048	0.260248742258431	0.228935748587690	0.273146684676011
    0.185048818439264	0.258008929738211	0.184111341471906	0.197020666103735	0.419237674218629	0.329583131194398	0.258738603216210	0.428719565972320	0.346697694484368	0.244090267580614	0.391549984585214	0.250255414156282	0.158490458520637	0.194274117433449	0.277696610296234	0.237582011612920	0.185048818439264	0.258008929738211	0.184111341471906	0.197020666103735	0.419237674218629	0.329583131194398	0.258738603216210	0.428719565972320	0.346697694484368	0.244090267580614	0.391549984585214	0.250255414156282	0.158490458520637	0.194274117433449	0.277696610296234	0.237582011612920	0.212129259364463	0.222888736595654	0.318898807571854	0.503458147248227	0.367316113266388	1	0.625482044460074	0.648679878691443	0.380768251978448	0.159245633927180	0.340769057534646	0.527966278310987	0.622308582378193	0.159278171243285	0.165307212791028	0.228766141810542	0.212129259364463	0.222888736595654	0.318898807571854	0.503458147248227	0.367316113266388	1	0.625482044460074	0.648679878691443	0.380768251978448	0.159245633927180	0.340769057534646	0.527966278310987	0.622308582378193	0.159278171243285	0.165307212791028	0.228766141810542
    0.238202610510690	0.256601897716442	0.174741938753492	0.284260054814097	0.488725641082833	0.406912056368750	0.237779012053277	0.489612201661955	0.397214414895051	0.225779144247316	0.371597503966509	0.211219450709919	0.191569207123187	0.199568957473710	0.263292195639382	0.286781069416484	0.238202610510690	0.256601897716442	0.174741938753492	0.284260054814097	0.488725641082833	0.406912056368750	0.237779012053277	0.489612201661955	0.397214414895051	0.225779144247316	0.371597503966509	0.211219450709919	0.191569207123187	0.199568957473710	0.263292195639382	0.286781069416484	0.168453869184571	0.250347801641747	0.323727296097083	0.578195376186041	0.550454917488211	0.625482044460074	1	0.787427738904120	0.488849534035235	0.133131753882770	0.574543025512610	0.544452393239199	0.723960474824533	0.227125840323131	0.147911624893899	0.348981055037568	0.168453869184571	0.250347801641747	0.323727296097083	0.578195376186041	0.550454917488211	0.625482044460074	1	0.787427738904120	0.488849534035235	0.133131753882770	0.574543025512610	0.544452393239199	0.723960474824533	0.227125840323131	0.147911624893899	0.348981055037568
    0.244657386305671	0.290654300814771	0.161328541375099	0.306282197159308	0.494811701913341	0.422318930549578	0.259812235085185	0.515545813640042	0.426482935033363	0.259439168982143	0.405186841441442	0.236420964685782	0.206721277320307	0.215356962511248	0.295776267565491	0.309485703130578	0.244657386305671	0.290654300814771	0.161328541375099	0.306282197159308	0.494811701913341	0.422318930549578	0.259812235085185	0.515545813640042	0.426482935033363	0.259439168982143	0.405186841441442	0.236420964685782	0.206721277320307	0.215356962511248	0.295776267565491	0.309485703130578	0.165063486104721	0.278881204611943	0.349661378890827	0.645342396741406	0.557751258890851	0.648679878691443	0.787427738904120	1	0.473766372691625	0.123271187683695	0.584899865748324	0.581411325805605	0.746588047839827	0.235302037429423	0.197305953135593	0.343325448828626	0.165063486104721	0.278881204611943	0.349661378890827	0.645342396741406	0.557751258890851	0.648679878691443	0.787427738904120	1	0.473766372691625	0.123271187683695	0.584899865748324	0.581411325805605	0.746588047839827	0.235302037429423	0.197305953135593	0.343325448828626
    0.262800581807974	0.319924976813934	0.132588234744717	0.267137279136037	0.414268088297245	0.330528228726873	0.297618588058598	0.355080589179632	0.449983936871318	0.334655469346739	0.365728800397777	0.204757998601009	0.189705473665356	0.146755406974320	0.263482123875376	0.205849622984209	0.262800581807974	0.319924976813934	0.132588234744717	0.267137279136037	0.414268088297245	0.330528228726873	0.297618588058598	0.355080589179632	0.449983936871318	0.334655469346739	0.365728800397777	0.204757998601009	0.189705473665356	0.146755406974320	0.263482123875376	0.205849622984209	0.168622540386407	0.345335936696517	0.0952671517244940	0.464914860750329	0.359265234826973	0.380768251978448	0.488849534035235	0.473766372691625	1	0.468147862419377	0.377863390845746	0.158254032823682	0.415123513075344	0.242973018072230	0.0927973466199387	0.335052662095399	0.168622540386407	0.345335936696517	0.0952671517244940	0.464914860750329	0.359265234826973	0.380768251978448	0.488849534035235	0.473766372691625	1	0.468147862419377	0.377863390845746	0.158254032823682	0.415123513075344	0.242973018072230	0.0927973466199387	0.335052662095399
    0.276177229421184	0.236311880485008	0.150571222907838	0.180893477778343	0.301894085772881	0.343883123574280	0.376420241926500	0.304338042723125	0.331400876549848	0.472056460404923	0.325338632802839	0.287411091053269	0.305240889445556	0.184024888821981	0.246628517932146	0.139276537033501	0.276177229421184	0.236311880485008	0.150571222907838	0.180893477778343	0.301894085772881	0.343883123574280	0.376420241926500	0.304338042723125	0.331400876549848	0.472056460404923	0.325338632802839	0.287411091053269	0.305240889445556	0.184024888821981	0.246628517932146	0.139276537033501	0.156279017501782	0.264580041902274	0.0871863300749362	0.151752731806341	0.169866450660821	0.159245633927180	0.133131753882770	0.123271187683695	0.468147862419377	1	0.141019735100310	0.124003277948560	0.104261779936692	0.152726428303252	0.131272218167682	0.163819925333801	0.156279017501782	0.264580041902274	0.0871863300749362	0.151752731806341	0.169866450660821	0.159245633927180	0.133131753882770	0.123271187683695	0.468147862419377	1	0.141019735100310	0.124003277948560	0.104261779936692	0.152726428303252	0.131272218167682	0.163819925333801
    0.336676352839481	0.193399058787936	0.155217995450231	0.362942511121388	0.368203777790144	0.322818832106714	0.194651469874191	0.391411058109650	0.301921963123240	0.197605829777781	0.234284764849756	0.201923713117651	0.349823875877798	0.203844847782750	0.201473876493201	0.284903071676509	0.336676352839481	0.193399058787936	0.155217995450231	0.362942511121388	0.368203777790144	0.322818832106714	0.194651469874191	0.391411058109650	0.301921963123240	0.197605829777781	0.234284764849756	0.201923713117651	0.349823875877798	0.203844847782750	0.201473876493201	0.284903071676509	0.0931167304667288	0.244921364491539	0.267043889089147	0.402303235121321	0.620217762992566	0.340769057534646	0.574543025512610	0.584899865748324	0.377863390845746	0.141019735100310	1	0.403702468364406	0.477293590020078	0.316068880836508	0.199059245516642	0.360291224650277	0.0931167304667288	0.244921364491539	0.267043889089147	0.402303235121321	0.620217762992566	0.340769057534646	0.574543025512610	0.584899865748324	0.377863390845746	0.141019735100310	1	0.403702468364406	0.477293590020078	0.316068880836508	0.199059245516642	0.360291224650277
    0.192605147348659	0.174715944086205	0.204002216638815	0.189051164175225	0.358403653653469	0.380010451159630	0.235946694846246	0.440688572535383	0.245949176015459	0.262489421725674	0.356423707451147	0.337905005742051	0.251851533104586	0.251985640799158	0.276081114681235	0.268899941359785	0.192605147348659	0.174715944086205	0.204002216638815	0.189051164175225	0.358403653653469	0.380010451159630	0.235946694846246	0.440688572535383	0.245949176015459	0.262489421725674	0.356423707451147	0.337905005742051	0.251851533104586	0.251985640799158	0.276081114681235	0.268899941359785	0.205201036957932	0.159471642239823	0.440808634948773	0.371449199942161	0.387660620005432	0.527966278310987	0.544452393239199	0.581411325805605	0.158254032823682	0.124003277948560	0.403702468364406	1	0.572196449299261	0.117528655794355	0.243072917849270	0.197024622410144	0.205201036957932	0.159471642239823	0.440808634948773	0.371449199942161	0.387660620005432	0.527966278310987	0.544452393239199	0.581411325805605	0.158254032823682	0.124003277948560	0.403702468364406	1	0.572196449299261	0.117528655794355	0.243072917849270	0.197024622410144
    0.210038103397686	0.251205002561756	0.226516656921773	0.245304855183173	0.463927096573370	0.375553378164681	0.212040428446836	0.474463527483964	0.364937088905025	0.204053926746631	0.430567652576515	0.281344602087222	0.189959042834344	0.250932620433989	0.290965230513544	0.304286170642738	0.210038103397686	0.251205002561756	0.226516656921773	0.245304855183173	0.463927096573370	0.375553378164681	0.212040428446836	0.474463527483964	0.364937088905025	0.204053926746631	0.430567652576515	0.281344602087222	0.189959042834344	0.250932620433989	0.290965230513544	0.304286170642738	0.181695384085607	0.223885938932743	0.365709285771055	0.583114246861958	0.384377072928048	0.622308582378193	0.723960474824533	0.746588047839827	0.415123513075344	0.104261779936692	0.477293590020078	0.572196449299261	1	0.239750890072918	0.136152337622859	0.340516329798941	0.181695384085607	0.223885938932743	0.365709285771055	0.583114246861958	0.384377072928048	0.622308582378193	0.723960474824533	0.746588047839827	0.415123513075344	0.104261779936692	0.477293590020078	0.572196449299261	1	0.239750890072918	0.136152337622859	0.340516329798941
    0.220212132020989	0.155368943598672	0.177749281176850	0.335138555543597	0.192235561370051	0.158512710264046	0.180430521800859	0.194559122028008	0.137260384502674	0.0868603132242725	0.201202640961916	0.241598343196672	0.171994402676464	0.306173680748464	0.288494413544832	0.350843129007248	0.220212132020989	0.155368943598672	0.177749281176850	0.335138555543597	0.192235561370051	0.158512710264046	0.180430521800859	0.194559122028008	0.137260384502674	0.0868603132242725	0.201202640961916	0.241598343196672	0.171994402676464	0.306173680748464	0.288494413544832	0.350843129007248	0.331381888723394	0.140128521425098	0.317878714930271	0.226871281996533	0.260248742258431	0.159278171243285	0.227125840323131	0.235302037429423	0.242973018072230	0.152726428303252	0.316068880836508	0.117528655794355	0.239750890072918	1	0.146455638641751	0.211402010776324	0.331381888723394	0.140128521425098	0.317878714930271	0.226871281996533	0.260248742258431	0.159278171243285	0.227125840323131	0.235302037429423	0.242973018072230	0.152726428303252	0.316068880836508	0.117528655794355	0.239750890072918	1	0.146455638641751	0.211402010776324
    0.176620548805626	0.125854532171030	0.177912875604621	0.157826223919072	0.148883560819950	0.178713286156517	0.216952897766093	0.196668230175339	0.0984872355958472	0.113088587100550	0.188398126858605	0.225351602112520	0.246868236409744	0.244212520546245	0.324120487388837	0.162236432601847	0.176620548805626	0.125854532171030	0.177912875604621	0.157826223919072	0.148883560819950	0.178713286156517	0.216952897766093	0.196668230175339	0.0984872355958472	0.113088587100550	0.188398126858605	0.225351602112520	0.246868236409744	0.244212520546245	0.324120487388837	0.162236432601847	0.0694982205575988	0.109383174686410	0.107976903089386	0.199785846755072	0.228935748587690	0.165307212791028	0.147911624893899	0.197305953135593	0.0927973466199387	0.131272218167682	0.199059245516642	0.243072917849270	0.136152337622859	0.146455638641751	1	0.137183258846915	0.0694982205575988	0.109383174686410	0.107976903089386	0.199785846755072	0.228935748587690	0.165307212791028	0.147911624893899	0.197305953135593	0.0927973466199387	0.131272218167682	0.199059245516642	0.243072917849270	0.136152337622859	0.146455638641751	1	0.137183258846915
    0.202230550437096	0.172439908547135	0.221049465091916	0.205458969300860	0.255283901846791	0.206601766747013	0.160941735562725	0.248086166274221	0.198903649734215	0.126335122236854	0.214103184190751	0.181966122544000	0.223233013220737	0.194598448235030	0.214320506259475	0.336763198993654	0.202230550437096	0.172439908547135	0.221049465091916	0.205458969300860	0.255283901846791	0.206601766747013	0.160941735562725	0.248086166274221	0.198903649734215	0.126335122236854	0.214103184190751	0.181966122544000	0.223233013220737	0.194598448235030	0.214320506259475	0.336763198993654	0.0937803258353237	0.185481661708238	0.147715667980565	0.175053651101674	0.273146684676011	0.228766141810542	0.348981055037568	0.343325448828626	0.335052662095399	0.163819925333801	0.360291224650277	0.197024622410144	0.340516329798941	0.211402010776324	0.137183258846915	1	0.0937803258353237	0.185481661708238	0.147715667980565	0.175053651101674	0.273146684676011	0.228766141810542	0.348981055037568	0.343325448828626	0.335052662095399	0.163819925333801	0.360291224650277	0.197024622410144	0.340516329798941	0.211402010776324	0.137183258846915	1
    0.134735352555133	0.173670491285841	0.355438791892943	0.138436483264120	0.165966452020543	0.160598464390537	0.187873073925526	0.161237062608672	0.142171564257541	0.144450005975650	0.185045963637839	0.200141851055337	0.160464678105875	0.266779427191495	0.184790520389430	0.188908568954180	0.134735352555133	0.173670491285841	0.355438791892943	0.138436483264120	0.165966452020543	0.160598464390537	0.187873073925526	0.161237062608672	0.142171564257541	0.144450005975650	0.185045963637839	0.200141851055337	0.160464678105875	0.266779427191495	0.184790520389430	0.188908568954180	1	0.119563571207939	0.134517649418353	0.145645342006557	0.124878680641713	0.212129259364463	0.168453869184571	0.165063486104721	0.168622540386407	0.156279017501782	0.0931167304667288	0.205201036957932	0.181695384085607	0.331381888723394	0.0694982205575988	0.0937803258353237	1	0.119563571207939	0.134517649418353	0.145645342006557	0.124878680641713	0.212129259364463	0.168453869184571	0.165063486104721	0.168622540386407	0.156279017501782	0.0931167304667288	0.205201036957932	0.181695384085607	0.331381888723394	0.0694982205575988	0.0937803258353237
    0.227214056458885	0.288216912365522	0.105011488754895	0.179991998384704	0.256167797179826	0.246367399139255	0.245990326829341	0.242600708434698	0.348498229969896	0.328899815080986	0.230621147338107	0.140964588095181	0.154103221330077	0.0869095597410673	0.195748174989036	0.130666084046134	0.227214056458885	0.288216912365522	0.105011488754895	0.179991998384704	0.256167797179826	0.246367399139255	0.245990326829341	0.242600708434698	0.348498229969896	0.328899815080986	0.230621147338107	0.140964588095181	0.154103221330077	0.0869095597410673	0.195748174989036	0.130666084046134	0.119563571207939	1	0.130300072561244	0.311823212613157	0.225668991451732	0.222888736595654	0.250347801641747	0.278881204611943	0.345335936696517	0.264580041902274	0.244921364491539	0.159471642239823	0.223885938932743	0.140128521425098	0.109383174686410	0.185481661708238	0.119563571207939	1	0.130300072561244	0.311823212613157	0.225668991451732	0.222888736595654	0.250347801641747	0.278881204611943	0.345335936696517	0.264580041902274	0.244921364491539	0.159471642239823	0.223885938932743	0.140128521425098	0.109383174686410	0.185481661708238
    0.153367199073079	0.158001480883866	0.156555090607610	0.196398065234584	0.237889094103677	0.247601458225063	0.184067479775407	0.282127225606543	0.166330902174613	0.172561189924276	0.223620976155853	0.222178343670103	0.168945275101099	0.271251740423712	0.143985932779031	0.231153257653851	0.153367199073079	0.158001480883866	0.156555090607610	0.196398065234584	0.237889094103677	0.247601458225063	0.184067479775407	0.282127225606543	0.166330902174613	0.172561189924276	0.223620976155853	0.222178343670103	0.168945275101099	0.271251740423712	0.143985932779031	0.231153257653851	0.134517649418353	0.130300072561244	1	0.293871915677903	0.252438920240734	0.318898807571854	0.323727296097083	0.349661378890827	0.0952671517244940	0.0871863300749362	0.267043889089147	0.440808634948773	0.365709285771055	0.317878714930271	0.107976903089386	0.147715667980565	0.134517649418353	0.130300072561244	1	0.293871915677903	0.252438920240734	0.318898807571854	0.323727296097083	0.349661378890827	0.0952671517244940	0.0871863300749362	0.267043889089147	0.440808634948773	0.365709285771055	0.317878714930271	0.107976903089386	0.147715667980565
    0.160072862647111	0.331153070144339	0.125346993220384	0.281801135230052	0.384173275699749	0.299630847277878	0.242997536756078	0.382027794502961	0.344051373683414	0.179401236672933	0.389218444585022	0.214474683217675	0.140221432295813	0.181456465610374	0.312044394341670	0.215809119317150	0.160072862647111	0.331153070144339	0.125346993220384	0.281801135230052	0.384173275699749	0.299630847277878	0.242997536756078	0.382027794502961	0.344051373683414	0.179401236672933	0.389218444585022	0.214474683217675	0.140221432295813	0.181456465610374	0.312044394341670	0.215809119317150	0.145645342006557	0.311823212613157	0.293871915677903	1	0.364173992674786	0.503458147248227	0.578195376186041	0.645342396741406	0.464914860750329	0.151752731806341	0.402303235121321	0.371449199942161	0.583114246861958	0.226871281996533	0.199785846755072	0.175053651101674	0.145645342006557	0.311823212613157	0.293871915677903	1	0.364173992674786	0.503458147248227	0.578195376186041	0.645342396741406	0.464914860750329	0.151752731806341	0.402303235121321	0.371449199942161	0.583114246861958	0.226871281996533	0.199785846755072	0.175053651101674
    0.298078652654840	0.194075382869973	0.127546712821935	0.342431532552267	0.368798215585550	0.353290735799105	0.233877443204227	0.392846364702732	0.315111126481233	0.237993322355195	0.222100397821173	0.164399475282110	0.329975873220386	0.184490367591133	0.186278817895424	0.263484082424131	0.298078652654840	0.194075382869973	0.127546712821935	0.342431532552267	0.368798215585550	0.353290735799105	0.233877443204227	0.392846364702732	0.315111126481233	0.237993322355195	0.222100397821173	0.164399475282110	0.329975873220386	0.184490367591133	0.186278817895424	0.263484082424131	0.124878680641713	0.225668991451732	0.252438920240734	0.364173992674786	1	0.367316113266388	0.550454917488211	0.557751258890851	0.359265234826973	0.169866450660821	0.620217762992566	0.387660620005432	0.384377072928048	0.260248742258431	0.228935748587690	0.273146684676011	0.124878680641713	0.225668991451732	0.252438920240734	0.364173992674786	1	0.367316113266388	0.550454917488211	0.557751258890851	0.359265234826973	0.169866450660821	0.620217762992566	0.387660620005432	0.384377072928048	0.260248742258431	0.228935748587690	0.273146684676011
    0.185048818439264	0.258008929738211	0.184111341471906	0.197020666103735	0.419237674218629	0.329583131194398	0.258738603216210	0.428719565972320	0.346697694484368	0.244090267580614	0.391549984585214	0.250255414156282	0.158490458520637	0.194274117433449	0.277696610296234	0.237582011612920	0.185048818439264	0.258008929738211	0.184111341471906	0.197020666103735	0.419237674218629	0.329583131194398	0.258738603216210	0.428719565972320	0.346697694484368	0.244090267580614	0.391549984585214	0.250255414156282	0.158490458520637	0.194274117433449	0.277696610296234	0.237582011612920	0.212129259364463	0.222888736595654	0.318898807571854	0.503458147248227	0.367316113266388	1	0.625482044460074	0.648679878691443	0.380768251978448	0.159245633927180	0.340769057534646	0.527966278310987	0.622308582378193	0.159278171243285	0.165307212791028	0.228766141810542	0.212129259364463	0.222888736595654	0.318898807571854	0.503458147248227	0.367316113266388	1	0.625482044460074	0.648679878691443	0.380768251978448	0.159245633927180	0.340769057534646	0.527966278310987	0.622308582378193	0.159278171243285	0.165307212791028	0.228766141810542
    0.238202610510690	0.256601897716442	0.174741938753492	0.284260054814097	0.488725641082833	0.406912056368750	0.237779012053277	0.489612201661955	0.397214414895051	0.225779144247316	0.371597503966509	0.211219450709919	0.191569207123187	0.199568957473710	0.263292195639382	0.286781069416484	0.238202610510690	0.256601897716442	0.174741938753492	0.284260054814097	0.488725641082833	0.406912056368750	0.237779012053277	0.489612201661955	0.397214414895051	0.225779144247316	0.371597503966509	0.211219450709919	0.191569207123187	0.199568957473710	0.263292195639382	0.286781069416484	0.168453869184571	0.250347801641747	0.323727296097083	0.578195376186041	0.550454917488211	0.625482044460074	1	0.787427738904120	0.488849534035235	0.133131753882770	0.574543025512610	0.544452393239199	0.723960474824533	0.227125840323131	0.147911624893899	0.348981055037568	0.168453869184571	0.250347801641747	0.323727296097083	0.578195376186041	0.550454917488211	0.625482044460074	1	0.787427738904120	0.488849534035235	0.133131753882770	0.574543025512610	0.544452393239199	0.723960474824533	0.227125840323131	0.147911624893899	0.348981055037568
    0.244657386305671	0.290654300814771	0.161328541375099	0.306282197159308	0.494811701913341	0.422318930549578	0.259812235085185	0.515545813640042	0.426482935033363	0.259439168982143	0.405186841441442	0.236420964685782	0.206721277320307	0.215356962511248	0.295776267565491	0.309485703130578	0.244657386305671	0.290654300814771	0.161328541375099	0.306282197159308	0.494811701913341	0.422318930549578	0.259812235085185	0.515545813640042	0.426482935033363	0.259439168982143	0.405186841441442	0.236420964685782	0.206721277320307	0.215356962511248	0.295776267565491	0.309485703130578	0.165063486104721	0.278881204611943	0.349661378890827	0.645342396741406	0.557751258890851	0.648679878691443	0.787427738904120	1	0.473766372691625	0.123271187683695	0.584899865748324	0.581411325805605	0.746588047839827	0.235302037429423	0.197305953135593	0.343325448828626	0.165063486104721	0.278881204611943	0.349661378890827	0.645342396741406	0.557751258890851	0.648679878691443	0.787427738904120	1	0.473766372691625	0.123271187683695	0.584899865748324	0.581411325805605	0.746588047839827	0.235302037429423	0.197305953135593	0.343325448828626
    0.262800581807974	0.319924976813934	0.132588234744717	0.267137279136037	0.414268088297245	0.330528228726873	0.297618588058598	0.355080589179632	0.449983936871318	0.334655469346739	0.365728800397777	0.204757998601009	0.189705473665356	0.146755406974320	0.263482123875376	0.205849622984209	0.262800581807974	0.319924976813934	0.132588234744717	0.267137279136037	0.414268088297245	0.330528228726873	0.297618588058598	0.355080589179632	0.449983936871318	0.334655469346739	0.365728800397777	0.204757998601009	0.189705473665356	0.146755406974320	0.263482123875376	0.205849622984209	0.168622540386407	0.345335936696517	0.0952671517244940	0.464914860750329	0.359265234826973	0.380768251978448	0.488849534035235	0.473766372691625	1	0.468147862419377	0.377863390845746	0.158254032823682	0.415123513075344	0.242973018072230	0.0927973466199387	0.335052662095399	0.168622540386407	0.345335936696517	0.0952671517244940	0.464914860750329	0.359265234826973	0.380768251978448	0.488849534035235	0.473766372691625	1	0.468147862419377	0.377863390845746	0.158254032823682	0.415123513075344	0.242973018072230	0.0927973466199387	0.335052662095399
    0.276177229421184	0.236311880485008	0.150571222907838	0.180893477778343	0.301894085772881	0.343883123574280	0.376420241926500	0.304338042723125	0.331400876549848	0.472056460404923	0.325338632802839	0.287411091053269	0.305240889445556	0.184024888821981	0.246628517932146	0.139276537033501	0.276177229421184	0.236311880485008	0.150571222907838	0.180893477778343	0.301894085772881	0.343883123574280	0.376420241926500	0.304338042723125	0.331400876549848	0.472056460404923	0.325338632802839	0.287411091053269	0.305240889445556	0.184024888821981	0.246628517932146	0.139276537033501	0.156279017501782	0.264580041902274	0.0871863300749362	0.151752731806341	0.169866450660821	0.159245633927180	0.133131753882770	0.123271187683695	0.468147862419377	1	0.141019735100310	0.124003277948560	0.104261779936692	0.152726428303252	0.131272218167682	0.163819925333801	0.156279017501782	0.264580041902274	0.0871863300749362	0.151752731806341	0.169866450660821	0.159245633927180	0.133131753882770	0.123271187683695	0.468147862419377	1	0.141019735100310	0.124003277948560	0.104261779936692	0.152726428303252	0.131272218167682	0.163819925333801
    0.336676352839481	0.193399058787936	0.155217995450231	0.362942511121388	0.368203777790144	0.322818832106714	0.194651469874191	0.391411058109650	0.301921963123240	0.197605829777781	0.234284764849756	0.201923713117651	0.349823875877798	0.203844847782750	0.201473876493201	0.284903071676509	0.336676352839481	0.193399058787936	0.155217995450231	0.362942511121388	0.368203777790144	0.322818832106714	0.194651469874191	0.391411058109650	0.301921963123240	0.197605829777781	0.234284764849756	0.201923713117651	0.349823875877798	0.203844847782750	0.201473876493201	0.284903071676509	0.0931167304667288	0.244921364491539	0.267043889089147	0.402303235121321	0.620217762992566	0.340769057534646	0.574543025512610	0.584899865748324	0.377863390845746	0.141019735100310	1	0.403702468364406	0.477293590020078	0.316068880836508	0.199059245516642	0.360291224650277	0.0931167304667288	0.244921364491539	0.267043889089147	0.402303235121321	0.620217762992566	0.340769057534646	0.574543025512610	0.584899865748324	0.377863390845746	0.141019735100310	1	0.403702468364406	0.477293590020078	0.316068880836508	0.199059245516642	0.360291224650277
    0.192605147348659	0.174715944086205	0.204002216638815	0.189051164175225	0.358403653653469	0.380010451159630	0.235946694846246	0.440688572535383	0.245949176015459	0.262489421725674	0.356423707451147	0.337905005742051	0.251851533104586	0.251985640799158	0.276081114681235	0.268899941359785	0.192605147348659	0.174715944086205	0.204002216638815	0.189051164175225	0.358403653653469	0.380010451159630	0.235946694846246	0.440688572535383	0.245949176015459	0.262489421725674	0.356423707451147	0.337905005742051	0.251851533104586	0.251985640799158	0.276081114681235	0.268899941359785	0.205201036957932	0.159471642239823	0.440808634948773	0.371449199942161	0.387660620005432	0.527966278310987	0.544452393239199	0.581411325805605	0.158254032823682	0.124003277948560	0.403702468364406	1	0.572196449299261	0.117528655794355	0.243072917849270	0.197024622410144	0.205201036957932	0.159471642239823	0.440808634948773	0.371449199942161	0.387660620005432	0.527966278310987	0.544452393239199	0.581411325805605	0.158254032823682	0.124003277948560	0.403702468364406	1	0.572196449299261	0.117528655794355	0.243072917849270	0.197024622410144
    0.210038103397686	0.251205002561756	0.226516656921773	0.245304855183173	0.463927096573370	0.375553378164681	0.212040428446836	0.474463527483964	0.364937088905025	0.204053926746631	0.430567652576515	0.281344602087222	0.189959042834344	0.250932620433989	0.290965230513544	0.304286170642738	0.210038103397686	0.251205002561756	0.226516656921773	0.245304855183173	0.463927096573370	0.375553378164681	0.212040428446836	0.474463527483964	0.364937088905025	0.204053926746631	0.430567652576515	0.281344602087222	0.189959042834344	0.250932620433989	0.290965230513544	0.304286170642738	0.181695384085607	0.223885938932743	0.365709285771055	0.583114246861958	0.384377072928048	0.622308582378193	0.723960474824533	0.746588047839827	0.415123513075344	0.104261779936692	0.477293590020078	0.572196449299261	1	0.239750890072918	0.136152337622859	0.340516329798941	0.181695384085607	0.223885938932743	0.365709285771055	0.583114246861958	0.384377072928048	0.622308582378193	0.723960474824533	0.746588047839827	0.415123513075344	0.104261779936692	0.477293590020078	0.572196449299261	1	0.239750890072918	0.136152337622859	0.340516329798941
    0.220212132020989	0.155368943598672	0.177749281176850	0.335138555543597	0.192235561370051	0.158512710264046	0.180430521800859	0.194559122028008	0.137260384502674	0.0868603132242725	0.201202640961916	0.241598343196672	0.171994402676464	0.306173680748464	0.288494413544832	0.350843129007248	0.220212132020989	0.155368943598672	0.177749281176850	0.335138555543597	0.192235561370051	0.158512710264046	0.180430521800859	0.194559122028008	0.137260384502674	0.0868603132242725	0.201202640961916	0.241598343196672	0.171994402676464	0.306173680748464	0.288494413544832	0.350843129007248	0.331381888723394	0.140128521425098	0.317878714930271	0.226871281996533	0.260248742258431	0.159278171243285	0.227125840323131	0.235302037429423	0.242973018072230	0.152726428303252	0.316068880836508	0.117528655794355	0.239750890072918	1	0.146455638641751	0.211402010776324	0.331381888723394	0.140128521425098	0.317878714930271	0.226871281996533	0.260248742258431	0.159278171243285	0.227125840323131	0.235302037429423	0.242973018072230	0.152726428303252	0.316068880836508	0.117528655794355	0.239750890072918	1	0.146455638641751	0.211402010776324
    0.176620548805626	0.125854532171030	0.177912875604621	0.157826223919072	0.148883560819950	0.178713286156517	0.216952897766093	0.196668230175339	0.0984872355958472	0.113088587100550	0.188398126858605	0.225351602112520	0.246868236409744	0.244212520546245	0.324120487388837	0.162236432601847	0.176620548805626	0.125854532171030	0.177912875604621	0.157826223919072	0.148883560819950	0.178713286156517	0.216952897766093	0.196668230175339	0.0984872355958472	0.113088587100550	0.188398126858605	0.225351602112520	0.246868236409744	0.244212520546245	0.324120487388837	0.162236432601847	0.0694982205575988	0.109383174686410	0.107976903089386	0.199785846755072	0.228935748587690	0.165307212791028	0.147911624893899	0.197305953135593	0.0927973466199387	0.131272218167682	0.199059245516642	0.243072917849270	0.136152337622859	0.146455638641751	1	0.137183258846915	0.0694982205575988	0.109383174686410	0.107976903089386	0.199785846755072	0.228935748587690	0.165307212791028	0.147911624893899	0.197305953135593	0.0927973466199387	0.131272218167682	0.199059245516642	0.243072917849270	0.136152337622859	0.146455638641751	1	0.137183258846915
    0.202230550437096	0.172439908547135	0.221049465091916	0.205458969300860	0.255283901846791	0.206601766747013	0.160941735562725	0.248086166274221	0.198903649734215	0.126335122236854	0.214103184190751	0.181966122544000	0.223233013220737	0.194598448235030	0.214320506259475	0.336763198993654	0.202230550437096	0.172439908547135	0.221049465091916	0.205458969300860	0.255283901846791	0.206601766747013	0.160941735562725	0.248086166274221	0.198903649734215	0.126335122236854	0.214103184190751	0.181966122544000	0.223233013220737	0.194598448235030	0.214320506259475	0.336763198993654	0.0937803258353237	0.185481661708238	0.147715667980565	0.175053651101674	0.273146684676011	0.228766141810542	0.348981055037568	0.343325448828626	0.335052662095399	0.163819925333801	0.360291224650277	0.197024622410144	0.340516329798941	0.211402010776324	0.137183258846915	1	0.0937803258353237	0.185481661708238	0.147715667980565	0.175053651101674	0.273146684676011	0.228766141810542	0.348981055037568	0.343325448828626	0.335052662095399	0.163819925333801	0.360291224650277	0.197024622410144	0.340516329798941	0.211402010776324	0.137183258846915	1
    ];
models{end}(65:80,:) = NaN; %No acoustic similarity in written only case
models{end}(:,65:80) = NaN;
models{end}(models{end}==1) = NaN; %Exclude self
models{end} = 1-models{end}; %Dissimilarity
this_model_name{end+1} = ['Acoustic multispectral coherence - noself'];


V = spm_vol(fullfile(GLMDir,'mask.nii'));
mask = spm_read_vols(V);
mask_index = results.mask_index;
downsamped_V.mat = V.mat;
downsamped_V.mat(1:3,1:3)=downsamped_V.mat(1:3,1:3)*downsamp_ratio;

clear results % to free memory

%% Save design matrices for visualisation
if save_design_matrices
    figure
    set(gcf,'Position',[100 100 1600 800]);
    set(gcf, 'PaperPositionMode', 'auto');
    set(gcf,'color','w');
    multivariate_matrix_dir = './multivariate_design matrices/';
    if ~exist(multivariate_matrix_dir,'dir')
        mkdir(multivariate_matrix_dir)
    end
    for m=1:length(this_model_name)
        if ~all(all(isnan(models{m}(1:64,1:64))))
            b = imagesc(models{m}(1:64,1:64),[floor(min(min(models{m}(1:64,1:64)))) ceil(max(max(models{m}(1:64,1:64))))]);
            set(b,'AlphaData',~isnan(models{m}(1:64,1:64)))
            axis square
            text(8.5,-1,'Match 3', 'HorizontalAlignment', 'center', 'fontweight', 'bold' )
            text(16+8.5,-1,'Match 15', 'HorizontalAlignment', 'center', 'fontweight', 'bold' )
            text(32+8.5,-1,'Mismatch 3', 'HorizontalAlignment', 'center', 'fontweight', 'bold' )
            text(48+8.5,-1,'Mismatch 15', 'HorizontalAlignment', 'center', 'fontweight', 'bold' )
            set(gca,'xtick',[0:16:64],'ytick',[0:16:64])
            drawnow
            saveas(gcf,[multivariate_matrix_dir this_model_name{m} '.pdf'])
            saveas(gcf,[multivariate_matrix_dir this_model_name{m} '.png'])
            colorbar
            drawnow
            saveas(gcf,[multivariate_matrix_dir this_model_name{m} '_withcolorbar.pdf'])
            saveas(gcf,[multivariate_matrix_dir this_model_name{m} '_withcolorbar.png'])
        end
    end
end

%% Now make effect maps
%parfor m=1:length(this_model_name) %Parallelising here impossible due to out of memory on serialisation unless data downsampled
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