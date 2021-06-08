function module_make_partial_roi_RSA(GLMDir,mask_names)
%For taking already calculated crossnobis distances and doing RSA

redo_maps = 0; %To re-caclulate maps

addpath('/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/RSA_scripts/es_scripts_fMRI')
addpath('/group/language/data/thomascope/7T_full_paradigm_pilot_analysis_scripts/RSA_scripts/decoding_toolbox_v3.999')
addpath(genpath('/group/language/data/ediz.sohoglu/matlab/rsatoolbox'));

for i = 1:length(mask_names)
    cfg.results.dir = fullfile(GLMDir,'TDTcrossnobis_ROI',mask_names{i});
    if ~exist(fullfile(cfg.results.dir,'res_other_average.mat'),'file')
        cfg.results.dir = fullfile(GLMDir,['TDTcrossnobis_ROI',mask_names{i}]); % Stupid coding error earlier in analysis led to misnamed directories
        if ~exist(fullfile(cfg.results.dir,'res_other_average.mat'),'file')
            disp([cfg.results.dir '/res_other_average.mat does not exist, moving on'])
            continue
        end
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
    
    outputDir = fullfile(cfg.results.dir,'RSA',version);
    
    if redo_maps
        if exist(outputDir,'dir'); rmdir(outputDir,'s'); mkdir(outputDir); else; mkdir(outputDir); end
    else
        if ~exist(outputDir,'dir'); mkdir(outputDir); end
    end
    
    clear models this_model_name
    
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
    
    basemodels.shared_segments_cross = circshift(basemodels.shared_segments,[8 0]); %I think this is correct, but need to 100% check the off-diagonals
    
    basemodelNames = {'vowels','shared_segments','shared_segments_cross'};
    
    load(fullfile(cfg.results.dir,'res_other_average.mat'));
    data = results.other_average.output;
    notempty_data = find(~cellfun(@isempty,results.other_average.output));
    modeltemplate = NaN(size(results.other_average.output{notempty_data(1)}));
    
    labelnames_denumbered = {};
    for i = 1:length(labelnames)
        labelnames_denumbered{i} = labelnames{i}(isletter(labelnames{i})|isspace(labelnames{i}));
    end
    modelNames = unique(labelnames_denumbered,'stable');
    
    MisMatch_Cross_decode_base = zeros(16,16);
    MisMatch_Cross_decode_base(9:17:end/2) = 1;
    MisMatch_Cross_decode_base(end/2+1:17:end) = 1;
    MisMatch_Cross_decode_base = 1-MisMatch_Cross_decode_base;
    
    Match_Cross_decode_base = 1-eye(16);
    
    % Now arrange models in cells each containing a pair to be partialled one
    % against the other (see Carlin ... Rowe 2011 Current Biology)
    
    basemodels.shared_segments_cross_noself = basemodels.shared_segments;
    basemodels.shared_segments_cross_noself(1:17:end) = NaN;
    basemodels.shared_segments_cross_noself = circshift(basemodels.shared_segments_cross_noself,[8 0]);
    
    basemodels.shared_segments(1:17:end) = NaN;
    basemodels.combined_SS = basemodels.shared_segments-basemodels.shared_segments_cross_noself;
    basemodels.shared_segments(1:17:end) = 0;
    basemodels.combined_SS = (basemodels.combined_SS +1)/2; %Scale zero to 1 - This actually makes no difference to correlations but helps visualisation
    
    basemodels.only_cross = basemodels.combined_SS;
    basemodels.only_cross(basemodels.shared_segments~=1) = NaN;
    
    basemodels.only_not_cross = basemodels.combined_SS;
    basemodels.only_not_cross(basemodels.shared_segments_cross_noself~=1) = NaN;
    
    cross_decode_label_pairs = {
        'Match Unclear', 'Mismatch Unclear';
        'Match Clear', 'Mismatch Unclear';
        'Match Unclear', 'Mismatch Clear';
        'Match Clear', 'Mismatch Clear';
        };
    
    % Model 1
    models{1}{1} = modeltemplate;
    this_model_name{1}{1} = ['M to MM combined_SS - no self'];
    for i = 1:size(cross_decode_label_pairs,1)
        models{1}{1}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.combined_SS;
        models{1}{1}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.combined_SS';
    end
    
    models{1}{2} = modeltemplate;
    this_model_name{1}{2} = ['M to MM only cross'];
    for i = 1:size(cross_decode_label_pairs,1)
        models{1}{2}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.only_cross;
        models{1}{2}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.only_cross';
    end
    
    % Next model
    models{end+1}{1} = modeltemplate;
    this_model_name{end+1}{1} = ['M to MM combined_SS - no self'];
    for i = 1:size(cross_decode_label_pairs,1)
        models{end}{1}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.combined_SS;
        models{end}{1}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.combined_SS';
    end
    models{end}{2} = modeltemplate;
    this_model_name{end}{2} = ['M to MM only not cross'];
    for i = 1:size(cross_decode_label_pairs,1)
        models{end}{2}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.only_not_cross;
        models{end}{2}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.only_not_cross';
    end
    
    % Next model
    models{end+1}{1} = modeltemplate;
    this_model_name{end+1}{1} = ['M to MM only cross'];
    for i = 1:size(cross_decode_label_pairs,1)
        models{end}{1}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.only_cross;
        models{end}{1}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.only_cross';
    end
    
    models{end}{2} = modeltemplate;
    this_model_name{end}{2} = ['M to MM only not cross'];
    for i = 1:size(cross_decode_label_pairs,1)
        models{end}{2}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.only_not_cross;
        models{end}{2}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.only_not_cross';
    end
    
    % Next model
    models{end+1}{1} = modeltemplate;
    this_model_name{end+1}{1} = ['M to MM Shared Segments Cross'];
    for i = 1:size(cross_decode_label_pairs,1)
        models{end}{1}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.shared_segments_cross_noself;
        models{end}{1}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.shared_segments_cross_noself';
    end
    
    models{end}{2} = modeltemplate;
    this_model_name{end}{2} = ['M to MM Shared Segments'];
    for i = 1:size(cross_decode_label_pairs,1)
        models{end}{2}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.shared_segments;
        models{end}{2}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.shared_segments';
    end
    
    % Next model
    models{end+1}{1} = modeltemplate;
    this_model_name{end+1}{1} = ['M to MM Shared Segments Cross'];
    for i = 1:size(cross_decode_label_pairs,1)
        models{end}{1}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.shared_segments_cross_noself;
        models{end}{1}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.shared_segments_cross_noself';
    end
    
    basemodels.shared_segments(1:17:end) = NaN;
    models{end}{2} = modeltemplate;
    this_model_name{end}{2} = ['M to MM Shared Segments - no self'];
    for i = 1:size(cross_decode_label_pairs,1)
        models{end}{2}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.shared_segments;
        models{end}{2}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.shared_segments';
    end
    basemodels.shared_segments(1:17:end) = 0;
    
    % Next model
    models{end+1}{1} = modeltemplate;
    this_model_name{end+1}{1} = ['M to MM Shared Segments Cross Negative'];
    for i = 1:size(cross_decode_label_pairs,1)
        models{end}{1}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = 1-basemodels.shared_segments_cross_noself;
        models{end}{1}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = 1-basemodels.shared_segments_cross_noself';
    end
    
    models{end}{2} = modeltemplate;
    this_model_name{end}{2} = ['M to MM Shared Segments'];
    for i = 1:size(cross_decode_label_pairs,1)
        models{end}{2}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.shared_segments;
        models{end}{2}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.shared_segments';
    end
    
    % Next model
    models{end+1}{1} = modeltemplate;
    this_model_name{end+1}{1} = ['M to MM Shared Segments Cross Negative'];
    for i = 1:size(cross_decode_label_pairs,1)
        models{end}{1}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = 1-basemodels.shared_segments_cross_noself;
        models{end}{1}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = 1-basemodels.shared_segments_cross_noself';
    end
    
    basemodels.shared_segments(1:17:end) = NaN;
    models{end}{2} = modeltemplate;
    this_model_name{end}{2} = ['M to MM Shared Segments - no self'];
    for i = 1:size(cross_decode_label_pairs,1)
        models{end}{2}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.shared_segments;
        models{end}{2}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.shared_segments';
    end
    basemodels.shared_segments(1:17:end) = 0;
    
    
    %
    %
    % cross_decode_label_pairs = {
    %     'Match Unclear', 'Mismatch Unclear';
    %     'Match Clear', 'Mismatch Unclear';
    %     'Match Unclear', 'Mismatch Clear';
    %     'Match Clear', 'Mismatch Clear';
    %     'Match Unclear', 'Match Clear';
    %     'Mismatch Unclear', 'Mismatch Clear';
    %     'Match Unclear', 'Written';
    %     'Match Clear', 'Written';
    %     'Mismatch Unclear', 'Written';
    %     'Mismatch Clear', 'Written'};
    %
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
    %
    % %Now attempt cross-condition shared segments RSA without cross decoding, recognising that the MisMatch cue
    % %was consistently 8 elements after/before the auditory word
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
    %     %                 imagesc(models{end},'AlphaData',~isnan(models{end}))
    %     %                 title(this_model_name{end})
    %     %                 colorbar
    %     %                 pause
    % end
    %
    % cross_decode_label_pairs = {
    %     'Match Unclear', 'Mismatch Unclear';
    %     'Match Clear', 'Mismatch Unclear';
    %     'Match Unclear', 'Mismatch Clear';
    %     'Match Clear', 'Mismatch Clear'
    %     'Match Unclear', 'Match Clear';
    %     'Mismatch Unclear', 'Mismatch Clear';
    %     'Match Unclear', 'Written';
    %     'Match Clear', 'Written';
    %     'Mismatch Unclear', 'Written';
    %     'Mismatch Clear', 'Written'
    %     'Match Unclear', 'Written';
    %     'Match Clear', 'Written';
    %     };
    %
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
    %
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
    %
    %
    % %Now attempt cross-condition shared segments RSA without cross decoding, recognising that the MisMatch cue
    % %was consistently 8 elements after/before the auditory word
    % basemodels.shared_segments(1:17:end) = NaN;
    % for i = 1:size(cross_decode_label_pairs,1)
    %     models{end+1} = modeltemplate;
    %     models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.shared_segments;
    %     models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.shared_segments';
    %     this_model_name{end+1} = [cross_decode_label_pairs{i,1} ' to ' cross_decode_label_pairs{i,2} ' SS_Match - no self'];
    %     %Optional check - view matrix
    %     %                     imagesc(models{end},'AlphaData',~isnan(models{end}))
    %     %                     title(this_model_name{end})
    %     %                     pause
    % end
    % basemodels.shared_segments(1:17:end) = 1;
    %
    % % Now add combined conditions
    % cross_decode_label_pairs = {
    %     'Match Unclear', 'Mismatch Unclear';
    %     'Match Clear', 'Mismatch Unclear';
    %     'Match Unclear', 'Mismatch Clear';
    %     'Match Clear', 'Mismatch Clear';
    %     'Match Unclear', 'Match Clear';
    %     'Mismatch Unclear', 'Mismatch Clear';};
    %
    % models{end+1} = modeltemplate;
    % this_model_name{end+1} = ['All spoken Cross-decode_Match'];
    % for i = 1:size(cross_decode_label_pairs,1)
    %     models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = Match_Cross_decode_base;
    %     models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = Match_Cross_decode_base';
    % end
    %
    % models{end+1} = modeltemplate;
    % this_model_name{end+1} = ['All spoken SS_Match'];
    % for i = 1:size(cross_decode_label_pairs,1)
    %     models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.shared_segments;
    %     models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.shared_segments';
    % end
    % models{end+1} = modeltemplate;
    % this_model_name{end+1} = ['All spoken SS_Match - no self'];
    % basemodels.shared_segments(1:17:end) = NaN;
    % for i = 1:size(cross_decode_label_pairs,1)
    %     models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.shared_segments;
    %     models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.shared_segments';
    % end
    % basemodels.shared_segments(1:17:end) = 0;
    %
    % cross_decode_label_pairs = {
    %     'Match Unclear', 'Written';
    %     'Match Clear', 'Written';
    %     'Mismatch Unclear', 'Written';
    %     'Mismatch Clear', 'Written'
    %     };
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
    %
    % % Now look at Match-Mismatch written word cross decoding
    % cross_decode_label_pairs = {
    %     'Match Unclear', 'Mismatch Unclear';
    %     'Match Clear', 'Mismatch Unclear';
    %     'Match Unclear', 'Mismatch Clear';
    %     'Match Clear', 'Mismatch Clear';
    %     };
    %
    % models{end+1} = modeltemplate;
    % this_model_name{end+1} = ['Match to Mismatch Shared Segments - no self'];
    % for i = 1:size(cross_decode_label_pairs,1)
    %     models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.shared_segments_cross_noself;
    %     models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.shared_segments_cross_noself';
    % end
    %
    % models{end+1} = modeltemplate;
    % this_model_name{end+1} = ['Match to Mismatch SS_Match - no self'];
    % basemodels.shared_segments(1:17:end) = NaN;
    % for i = 1:size(cross_decode_label_pairs,1)
    %     models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.shared_segments;
    %     models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.shared_segments';
    % end
    % basemodels.shared_segments(1:17:end) = 0;
    
    
    % % Now look at these three models Clear-Unclear
    % % then in every individual combination
    % cross_decode_label_pairs = {
    %     'Match Unclear', 'Match Clear';
    %     'Mismatch Unclear', 'Mismatch Clear';
    %     'Match Unclear', 'Mismatch Clear';
    %     'Match Clear', 'Mismatch Unclear';
    %     };
    % models{end+1} = modeltemplate;
    % this_model_name{end+1} = ['Clear to Unclear combined_SS - no self - rescaled'];
    % for i = 1:size(cross_decode_label_pairs,1)
    %     models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.combined_SS;
    %     models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.combined_SS';
    % end
    %
    % models{end+1} = modeltemplate;
    % this_model_name{end+1} = ['Clear to Unclear only cross'];
    % for i = 1:size(cross_decode_label_pairs,1)
    %     models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.only_cross;
    %     models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.only_cross';
    % end
    %
    % models{end+1} = modeltemplate;
    % this_model_name{end+1} = ['Clear to Unclear only not cross'];
    % for i = 1:size(cross_decode_label_pairs,1)
    %     models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.only_not_cross;
    %     models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.only_not_cross';
    % end
    %
    % cross_decode_label_pairs = {
    %     'Match Unclear', 'Mismatch Unclear';
    %     'Match Clear', 'Mismatch Unclear';
    %     'Match Unclear', 'Mismatch Clear';
    %     'Match Clear', 'Mismatch Clear'
    %     'Match Unclear', 'Match Clear';
    %     'Mismatch Unclear', 'Mismatch Clear';
    %     };
    %
    % for i = 1:size(cross_decode_label_pairs,1)
    %     models{end+1} = modeltemplate;
    %     models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.combined_SS;
    %     models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.combined_SS';
    %     this_model_name{end+1} = [cross_decode_label_pairs{i,1} ' to ' cross_decode_label_pairs{i,2} ' combined_SS - no self'];
    %     %Optional check - view matrix
    %     %             imagesc(models{end},'AlphaData',~isnan(models{end}))
    %     %             title(this_model_name{end})
    %     %             pause
    % end
    %
    % for i = 1:size(cross_decode_label_pairs,1)
    %     models{end+1} = modeltemplate;
    %     models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.only_cross;
    %     models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.only_cross';
    %     this_model_name{end+1} = [cross_decode_label_pairs{i,1} ' to ' cross_decode_label_pairs{i,2} ' only cross'];
    %     %Optional check - view matrix
    %     %             imagesc(models{end},'AlphaData',~isnan(models{end}))
    %     %             title(this_model_name{end})
    %     %             pause
    % end
    %
    % for i = 1:size(cross_decode_label_pairs,1)
    %     models{end+1} = modeltemplate;
    %     models{end}(strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered)) = basemodels.only_not_cross;
    %     models{end}(strcmp(cross_decode_label_pairs{i,2},labelnames_denumbered),strcmp(cross_decode_label_pairs{i,1},labelnames_denumbered)) = basemodels.only_not_cross';
    %     this_model_name{end+1} = [cross_decode_label_pairs{i,1} ' to ' cross_decode_label_pairs{i,2} ' only not cross'];
    %     %Optional check - view matrix
    %     %             imagesc(models{end},'AlphaData',~isnan(models{end}))
    %     %             title(this_model_name{end})
    %     %             pause
    % end
    
    % % models = modelRDMs; close all
    % %modelNames = fieldnames(models);
    %
    % %modelNames = {'ProbM' 'ProbMM' 'EntropyM' 'EntropyMM'};

    % % Optional - to visualise a model
    % b = imagesc(models{end}{1},[0 1]);
    % set(b,'AlphaData',~isnan(models{end}{1}))
    % title(this_model_name{end}{1},'Interpreter','none')
    % colorbar
    % pause
    % b = imagesc(models{end}{2},[0 1]);
    % set(b,'AlphaData',~isnan(models{end}{2}))
    % title(this_model_name{end}{2},'Interpreter','none')
    % colorbar
    
    
    roi_names = results.roi_names;
    for m=1:length(this_model_name)
        fprintf('\nComputing ROI effects for model %s partialling out model %s\n',this_model_name{m}{1},this_model_name{m}{2});
        shorter_string = min(numel(this_model_name{m}{1}),numel(this_model_name{m}{2}));
        common_letters = sum(cumprod(this_model_name{m}{1}(1:shorter_string)==this_model_name{m}{2}(1:shorter_string)));
        this_out_string{1} = [deblank(this_model_name{m}{1}(1:common_letters)) ': ' deblank(this_model_name{m}{1}(common_letters+1:end)) ' partialling ' deblank(this_model_name{m}{2}(common_letters+1:end))];
        this_out_string{2} = [deblank(this_model_name{m}{1}(1:common_letters)) ': ' deblank(this_model_name{m}{2}(common_letters+1:end)) ' partialling ' deblank(this_model_name{m}{1}(common_letters+1:end))];
        if ~exist(fullfile(outputDir,['roi_effects_' this_out_string{2} '.mat'])) || redo_maps == 1
            modelRDM = vectorizeRDMs(models{m}{1})';
            modeloutRDM = vectorizeRDMs(models{m}{2})';
            for vx=1:numel(data)
                neuralRDM = vectorizeRDMs(data{vx})';
                if ~isempty(strfind(version,'pearson'))
                    effect_1(vx) = fisherTransform(partialcorr(modelRDM,neuralRDM,modeloutRDM,'type','Pearson','Rows','pairwise'));
                    effect_2(vx) = fisherTransform(partialcorr(modeloutRDM,neuralRDM,modelRDM,'type','Pearson','Rows','pairwise'));
                elseif ~isempty(strfind(version,'spearman'))
                    effect_1(vx) = fisherTransform(partialcorr(modelRDM,neuralRDM,modeloutRDM,'type','Spearman','Rows','pairwise'));
                    effect_2(vx) = fisherTransform(partialcorr(modeloutRDM,neuralRDM,modelRDM,'type','Spearman','Rows','pairwise'));
                end
            end
            for direction = 1:2
                eval(['roi_effect = effect_' num2str(direction) ';'])
                save(fullfile(outputDir,['roi_effects_' this_out_string{direction} '.mat']),'roi_effect','roi_names')
            end
        else
            disp('Already exists - moving on - set redo_maps to 1 if you want to re-make')
        end
    end
end