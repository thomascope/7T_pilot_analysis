%% Analyse by condition and brain region
addpath(genpath('/imaging/mlr/users/tc02/toolboxes')); %Where is the RSA toolbox?

% masks = { %rw for re-sliced after warping into native space (I have also re-binarised). Underscores for spaces in the atlas search above.
%     'rwblank_mask'
%     'rwLeft_STG'
%     'rwLeft_PT'
%     'rwLeft_PrG'
%     'rwLeft_FO'
%     'rwLeft_TrIFG'
%     };

masks = {
    'rwBlank_2016_inflated'
    'rwL_STG_cross-segment_cluster'
    'rwLeft_Superior_Temporal_Gyrus'
    'rwLeft_Angular_Gyrus'
    'rwLeft_Precentral_Gyrus'
    'rwLeft_Frontal_Operculum'
    'rwLeft_Inferior_Frontal_Angular_Gyrus'
    'rwRight_Superior_Temporal_Gyrus'
    'rwRight_Angular_Gyrus'
    'rwRight_Precentral_Gyrus'
    'rwRight_Frontal_Operculum'
    'rwRight_Inferior_Frontal_Angular_Gyrus'
    'rwLeft_IFG_Written_Cluster'
    'rwLeft_Precentral_Written_Cluster'
    };

GLMDir = [preprocessedpathstem subjects{1} '/stats4_multi_3_noabsent']; %Template, first subject
temp = load([GLMDir filesep 'SPM.mat']);
labelnames = {};
for i = 1:length(temp.SPM.Sess(1).U)
    if ~strncmp(temp.SPM.Sess(1).U(i).name,{'Match','Mismatch','Written'},5)
        continue
    else
        labelnames(end+1) = temp.SPM.Sess(1).U(i).name;
    end
end        
labelnames_denumbered = {};
for i = 1:length(labelnames)
    labelnames_denumbered{i} = labelnames{i}(isletter(labelnames{i})|isspace(labelnames{i}));
end
conditionnames = unique(labelnames_denumbered,'stable');
clear temp labelnames_denumbered labelnames

nrun = size(subjects,2); % enter the number of runs here
mahalanobisroiworkedcorrectly = zeros(1,nrun);
parfor crun = 1:nrun
    addpath(genpath('./RSA_scripts'))
    GLMDir = [preprocessedpathstem subjects{crun} '/stats4_multi_3_noabsent']; %Where is the SPM model?
    mask_dir = [preprocessedpathstem subjects{crun}]; %Where are the native space ROI masks?
    try
        TDTCrossnobisAnalysis_roi(GLMDir,mask_dir,masks);
        mahalanobisroiworkedcorrectly(crun) = 1;
    catch
        mahalanobisroiworkedcorrectly(crun) = 0;
    end
end

%% Now do RSA on ROI data
nrun = size(subjects,2); % enter the number of runs here
RSAroiworkedcorrectly = zeros(1,nrun);
masks = {
    'rwBlank_2016_inflated'
    'rwL_STG_cross-segment_cluster'
    'rwLeft_Superior_Temporal_Gyrus'
    'rwLeft_Angular_Gyrus'
    'rwLeft_Precentral_Gyrus'
    'rwLeft_Frontal_Operculum'
    'rwLeft_Inferior_Frontal_Angular_Gyrus'
    'rwRight_Superior_Temporal_Gyrus'
    'rwRight_Angular_Gyrus'
    'rwRight_Precentral_Gyrus'
    'rwRight_Frontal_Operculum'
    'rwRight_Inferior_Frontal_Angular_Gyrus'
    'rwLeft_IFG_Written_Cluster'
    'rwLeft_Precentral_Written_Cluster'
    };

parfor crun = 1:nrun
    addpath(genpath('./RSA_scripts'))
    GLMDir = [preprocessedpathstem subjects{crun} '/stats4_multi_3_noabsent']; %Where is the SPM model?
    try
        module_roi_RSA(GLMDir,masks)
        RSAroiworkedcorrectly(crun) = 1;
    catch
        RSAroiworkedcorrectly(crun) = 0;
    end
end

%% Analyse by condition and brain region
addpath(genpath('/imaging/mlr/users/tc02/toolboxes')); %Where is the RSA toolbox?

% masks = { %rw for re-sliced after warping into native space (I have also re-binarised). Underscores for spaces in the atlas search above.
%     'rwblank_mask'
%     'rwLeft_STG'
%     'rwLeft_PT'
%     'rwLeft_PrG'
%     'rwLeft_FO'
%     'rwLeft_TrIFG'
%     };

masks = {
    'rwBlank_2016_inflated'
    'rwL_STG_cross-segment_cluster'
    'rwLeft_Superior_Temporal_Gyrus'
    'rwLeft_Angular_Gyrus'
    'rwLeft_Precentral_Gyrus'
    'rwLeft_Frontal_Operculum'
    'rwLeft_Inferior_Frontal_Angular_Gyrus'
    'rwRight_Superior_Temporal_Gyrus'
    'rwRight_Angular_Gyrus'
    'rwRight_Precentral_Gyrus'
    'rwRight_Frontal_Operculum'
    'rwRight_Inferior_Frontal_Angular_Gyrus'
    'rwLeft_IFG_Written_Cluster'
    'rwLeft_Precentral_Written_Cluster'
    };

GLMDir = [preprocessedpathstem subjects{1} '/stats4_multi_3']; %Template, first subject
temp = load([GLMDir filesep 'SPM.mat']);
labelnames = {};
for i = 1:length(temp.SPM.Sess(1).U)
    if ~strncmp(temp.SPM.Sess(1).U(i).name,{'Match','Mismatch','Written'},5)
        continue
    else
        labelnames(end+1) = temp.SPM.Sess(1).U(i).name;
    end
end        
labelnames_denumbered = {};
for i = 1:length(labelnames)
    labelnames_denumbered{i} = labelnames{i}(isletter(labelnames{i})|isspace(labelnames{i}));
end
conditionnames = unique(labelnames_denumbered,'stable');
clear temp labelnames_denumbered labelnames

nrun = size(subjects,2); % enter the number of runs here
mahalanobisroiworkedcorrectly = zeros(1,nrun);
parfor crun = 1:nrun
    addpath(genpath('./RSA_scripts'))
    GLMDir = [preprocessedpathstem subjects{crun} '/stats4_multi_3']; %Where is the SPM model?
    mask_dir = [preprocessedpathstem subjects{crun}]; %Where are the native space ROI masks?
    try
        TDTCrossnobisAnalysis_roi(GLMDir,mask_dir,masks);
        mahalanobisroiworkedcorrectly(crun) = 1;
    catch
        mahalanobisroiworkedcorrectly(crun) = 0;
    end
end

%% Now do RSA on ROI data
nrun = size(subjects,2); % enter the number of runs here
RSAroiworkedcorrectly = zeros(1,nrun);
masks = {
    'rwBlank_2016_inflated'
    'rwL_STG_cross-segment_cluster'
    'rwLeft_Superior_Temporal_Gyrus'
    'rwLeft_Angular_Gyrus'
    'rwLeft_Precentral_Gyrus'
    'rwLeft_Frontal_Operculum'
    'rwLeft_Inferior_Frontal_Angular_Gyrus'
    'rwRight_Superior_Temporal_Gyrus'
    'rwRight_Angular_Gyrus'
    'rwRight_Precentral_Gyrus'
    'rwRight_Frontal_Operculum'
    'rwRight_Inferior_Frontal_Angular_Gyrus'
    'rwLeft_IFG_Written_Cluster'
    'rwLeft_Precentral_Written_Cluster'
    };

parfor crun = 1:nrun
    addpath(genpath('./RSA_scripts'))
    GLMDir = [preprocessedpathstem subjects{crun} '/stats4_multi_3']; %Where is the SPM model?
    try
        module_roi_RSA(GLMDir,masks)
        RSAroiworkedcorrectly(crun) = 1;
    catch
        RSAroiworkedcorrectly(crun) = 0;
    end
end

%% Compare across conditions in STG as a sanity check then go on to do all ROIs
this_spm_model = 'stats4_multi_3_nowritten2';
%this_spm_model = 'stats4_multi_3';
%this_spm_model = 'stats4_multi_3_noabsent';

GLMDir = [preprocessedpathstem subjects{1} '/' this_spm_model]; %Template, first subject
temp = load([GLMDir filesep 'SPM.mat']);
labelnames = {};
for i = 1:length(temp.SPM.Sess(1).U)
    if ~strncmp(temp.SPM.Sess(1).U(i).name,{'Match','Mismatch','Written'},5)
        continue
    else
        labelnames(end+1) = temp.SPM.Sess(1).U(i).name;
    end
end        
labelnames_denumbered = {};
for i = 1:length(labelnames)
    labelnames_denumbered{i} = labelnames{i}(isletter(labelnames{i})|isspace(labelnames{i}));
end
conditionnames = unique(labelnames_denumbered,'stable');

clear this_model_name mask_names
this_model_name{1} = {
    'Match Unclear vowels'
    'Match Clear vowels'
    'Mismatch Unclear vowels'
    'Mismatch Clear vowels'
    'Match Unclear shared_segments'
    'Match Clear shared_segments'
    'Mismatch Unclear shared_segments'
    'Mismatch Clear shared_segments'
    'Written vowels'
    'Written shared_segments'
    };

this_model_name{2} = {
    'Match Unclear to Mismatch Unclear Cross-decode_Match'
    'Match Unclear to Mismatch Unclear SS_Match'
    'Match Unclear to Mismatch Unclear SS_Match - no self'
    'Match Unclear to Mismatch Clear Cross-decode_Match'
    'Match Unclear to Mismatch Clear SS_Match'
    'Match Unclear to Mismatch Clear SS_Match - no self'
    'Match Clear to Mismatch Unclear Cross-decode_Match'
    'Match Clear to Mismatch Unclear SS_Match'
    'Match Clear to Mismatch Unclear SS_Match - no self'
    'Match Clear to Mismatch Clear Cross-decode_Match'
    'Match Clear to Mismatch Clear SS_Match'
    'Match Clear to Mismatch Clear SS_Match - no self'
    'Match Unclear to Written Cross-decode_Match'
    'Match Clear to Written Cross-decode_Match'
    'Mismatch Unclear to Written Cross-decode_Match'
    'Mismatch Clear to Written Cross-decode_Match'
    };

this_model_name{3} = {
    'Match Unclear to Mismatch Unclear Cross-decode'
    'Match Unclear to Mismatch Unclear Shared Segments - cross'
    'Match Unclear to Mismatch Unclear Shared Segments - no self'
    'Match Unclear to Mismatch Clear Cross-decode'
    'Match Unclear to Mismatch Clear Shared Segments - cross'
    'Match Unclear to Mismatch Clear Shared Segments - no self'
    'Match Clear to Mismatch Unclear Cross-decode'
    'Match Clear to Mismatch Unclear Shared Segments - cross'
    'Match Clear to Mismatch Unclear Shared Segments - no self'
    'Match Clear to Mismatch Clear Cross-decode'
    'Match Clear to Mismatch Clear Shared Segments - cross'
    'Match Clear to Mismatch Clear Shared Segments - no self'
    'Match Unclear to Written Cross-decode'
    'Match Clear to Written Cross-decode'
    'Mismatch Unclear to Written Cross-decode'
    'Mismatch Clear to Written Cross-decode'
    };

nrun = size(subjects,2); % enter the number of runs here
% First load in the similarities
RSA_ROI_data_exist = zeros(1,nrun);
all_data = [];
mask_names{1} = {'rwLeft_Superior_Temporal_Gyrus';
    %'rwL_STG_cross-segment_cluster'
    'rwBlank_2016_inflated'
    'rwLeft_IFG_Written_Cluster'
    'rwLeft_Precentral_Written_Cluster'};

% mask_names{2} = {    
%     %'rwL_STG_cross-segment_cluster'
%     'rwLeft_Angular_Gyrus'
%     'rwLeft_Precentral_Gyrus'
%     'rwLeft_Frontal_Operculum'
%     'rwLeft_Inferior_Frontal_Angular_Gyrus'
%     'rwRight_Superior_Temporal_Gyrus'
%     'rwRight_Angular_Gyrus'
%     'rwRight_Precentral_Gyrus'
%     'rwRight_Frontal_Operculum'
%     'rwRight_Inferior_Frontal_Angular_Gyrus'
%     };

for j = 1:length(this_model_name)
    for k = 1:length(mask_names)
        for i = 1:length(mask_names{k})
            for crun = 1:nrun
                %ROI_RSA_dir = [preprocessedpathstem subjects{crun} '/stats4_multi_3_nowritten2/TDTcrossnobis_ROI/RSA/spearman']; %Where are the results>
                ROI_RSA_dir = [preprocessedpathstem subjects{crun} '/' this_spm_model '/TDTcrossnobis_ROI/' mask_names{k}{i} '/RSA/spearman'];
                if ~exist(fullfile(ROI_RSA_dir,['roi_effects_' this_model_name{j}{1} '.mat']),'file')
                    ROI_RSA_dir = [preprocessedpathstem subjects{crun} '/' this_spm_model '/TDTcrossnobis_ROI' mask_names{k}{i} '/RSA/spearman']; % Stupid coding error earlier in analysis led to misnamed directories
                end
                for m = 1:length(this_model_name{j})
                    try
                        temp_data = load(fullfile(ROI_RSA_dir,['roi_effects_' this_model_name{j}{m} '.mat']));
                        all_data(m,:,crun) = temp_data.roi_effect; %Create a matrix of condition by ROI by subject
                        RSA_ROI_data_exist(crun) = 1;
                    catch
                        warning(['No data for ' subjects{crun} ' probably because of SPM dropout, ignoring them'])
                        %error
                        RSA_ROI_data_exist(crun) = 0;
                        continue
                    end
                end
            end
            roi_names = temp_data.roi_names;
            clear temp_data
            disp(['Excluding subjects ' num2str(find(RSA_ROI_data_exist==0)) ' belonging to groups ' num2str(group(RSA_ROI_data_exist==0)) ' maybe check them'])
            all_data(:,:,RSA_ROI_data_exist==0) = NaN;
            
            %this_ROI = find(strcmp('rwLeft_Superior_Temporal_Gyrus',roi_names));
            this_ROI = find(strcmp(mask_names{k}{i},roi_names));
            figure
            set(gcf,'Position',[100 100 1600 800]);
            set(gcf, 'PaperPositionMode', 'auto');
            hold on
            errorbar([1:length(this_model_name{j})]-0.1,nanmean(squeeze(all_data(:,this_ROI,group==1&RSA_ROI_data_exist)),2),nanstd(squeeze(all_data(:,this_ROI,group==1&RSA_ROI_data_exist))')/sqrt(sum(group==1&RSA_ROI_data_exist)),'kx')
            errorbar([1:length(this_model_name{j})]+0.1,nanmean(squeeze(all_data(:,this_ROI,group==2&RSA_ROI_data_exist)),2),nanstd(squeeze(all_data(:,this_ROI,group==2&RSA_ROI_data_exist))')/sqrt(sum(group==2&RSA_ROI_data_exist)),'rx')
            xlim([0 length(this_model_name{j})+1])
            set(gca,'xtick',[1:length(this_model_name{j})],'xticklabels',this_model_name{j},'XTickLabelRotation',45,'TickLabelInterpreter','none')
            plot([0 length(this_model_name{j})+1],[0,0],'k--')
            title([mask_names{k}{i}(3:end) ' RSA'],'Interpreter','none')
            legend('Controls','Patients','location','southeast','AutoUpdate','off')
            [h,p] = ttest(squeeze(all_data(:,this_ROI,logical(RSA_ROI_data_exist)))');
            these_y_lims = ylim;
            if sum(h)~=0
                plot(find(h),these_y_lims(2)-diff(these_y_lims/10),'k*')
            end
            [h,p] = ttest2(squeeze(all_data(:,this_ROI,group==1&logical(RSA_ROI_data_exist)))',squeeze(all_data(:,this_ROI,group==2&logical(RSA_ROI_data_exist)))');
            if sum(h)~=0
                plot(find(h),these_y_lims(2)-diff(these_y_lims/20),'kx')
            end
            drawnow
        end
    end
end
