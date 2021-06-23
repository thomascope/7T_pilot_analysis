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
    'rwLeft_PrG_All_Shared_Segments'
    'rwLeft_PrG_SSMatchnoself_combined'
    'rwLeft_STG_Univariate8mm_15>3'
    'rwLeft_STG_Univariate3mm_15>3'
    'rwLeft_Precentral_Univariate_Interaction_combined'
    'rwLeft_Angular_Univariate_Interaction_combined'
    'rwLeft_PostSTG_Univariate_Interaction'
    'rwLeft_Precentral_Univariate_Interaction1'
    'rwLeft_Precentral_Univariate_Interaction2'
    'rwLeft_Precentral_Univariate_Interaction3'
    'rwLeft_Angular_Univariate_Interaction1'
    'rwLeft_Angular_Univariate_Interaction2'
    'rwLeft_Frontal_Univariate_MM>M'
    'rwLeft_Temporal_Univariate_MM>M'
    'rwLeft_IFG_cross_group_cluster'
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
partialRSAroiworkedcorrectly = zeros(1,nrun);
masks = {
    'rwLeft_PrG_All_Shared_Segments'
    'rwLeft_PrG_SSMatchnoself_combined'
    'rwLeft_STG_Univariate8mm_15>3'
    'rwLeft_STG_Univariate3mm_15>3'
    'rwLeft_Precentral_Univariate_Interaction_combined'
    'rwLeft_Angular_Univariate_Interaction_combined'
    'rwLeft_PostSTG_Univariate_Interaction'
    'rwLeft_Precentral_Univariate_Interaction1'
    'rwLeft_Precentral_Univariate_Interaction2'
    'rwLeft_Precentral_Univariate_Interaction3'
    'rwLeft_Angular_Univariate_Interaction1'
    'rwLeft_Angular_Univariate_Interaction2'
    'rwLeft_Frontal_Univariate_MM>M'
    'rwLeft_Temporal_Univariate_MM>M'
    'rwLeft_IFG_cross_group_cluster'
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
parfor crun = 1:nrun
    addpath(genpath('./RSA_scripts'))
    GLMDir = [preprocessedpathstem subjects{crun} '/stats4_multi_3']; %Where is the SPM model?
    try
        module_make_partial_roi_RSA(GLMDir,masks)
        partialRSAroiworkedcorrectly(crun) = 1;
    catch
        partialRSAroiworkedcorrectly(crun) = 0;
    end
end


%% Compare across conditions in STG as a sanity check then go on to do all ROIs
GLMDir = [preprocessedpathstem subjects{1} '/stats4_multi_3']; %Template, first subject
outdir = ['./ROI_figures/stats4_multi_3'];
mkdir(outdir)
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

% Add covariates of interest
model_run_date = '13-May-2021';
try
    load(['./modelparameters/modelparameters_' model_run_date '.mat'])
catch
    [all_sigma_pred,all_thresholds,controls_sigma_pred,controls_threshold,patients_sigma_pred,patients_threshold] = module_bayesian_behaviour(subjects,group,dates);
    save(['./modelparameters/modelparameters_' date '.mat'],'all_sigma_pred','all_thresholds','controls_sigma_pred','controls_threshold','patients_sigma_pred','patients_threshold');
end
this_age = [];
age_lookup = readtable('Pinfa_ages.csv');

extract_run_date = '22-Jun-2021';
try
    load(['./freesurfer_stats/roi_thicknesses_' extract_run_date '.mat'])
catch
    setenv('SUBJECTS_DIR',this_subjects_dir);
    all_roi_thicknesses = module_extract_freesurfer(Regions_of_interest,subjects,group);
    save(['./freesurfer_stats/roi_thicknesses_' date '.mat'],'all_roi_thicknesses');
end

for crun = 1:length(subjects)
    this_age(crun) = age_lookup.Age(strcmp(age_lookup.Study_ID,subjects{crun}));
end
covariates = [this_age',nanmean(all_sigma_pred)',all_roi_thicknesses{:,:}];
covariate_names = horzcat('Age','Prior_Precision',all_roi_thicknesses.Properties.VariableNames);

%Now build model space for testing

%Now build model space for testing
clear this_model_name mask_names
this_model_name{1} = {
    'Match Unclear shared_segments'
%     'Match Unclear shared_segments_mismatch'
%     'Match Unclear shared_segments_both'
    'Match Clear shared_segments'
%     'Match Clear shared_segments_mismatch'
%     'Match Clear shared_segments_both'
%     'Mismatch Unclear shared_segments'
    'Mismatch Unclear shared_segments_mismatch'
%     'Mismatch Unclear shared_segments_both'
%     'Mismatch Clear shared_segments'
    'Mismatch Clear shared_segments_mismatch'
%     'Mismatch Clear shared_segments_both'
    %     'Written vowels'
    'Written shared_segments'
    };

% this_model_name{2} = {
%     'Match Clear shared_segments: _mismatch partialling '
%     'Match Clear: shared_segments_mismatch partialling vowels'
%     'Match Clear shared_segments:  partialling _mismatch'
%     'Match Clear: shared_segments partialling vowels'
%     'Match Clear: vowels partialling shared_segments'
%     'Match Clear: vowels partialling shared_segments_mismatch'
%     'Match Unclear shared_segments: _mismatch partialling '
%     'Match Unclear: shared_segments_mismatch partialling vowels'
%     'Match Unclear shared_segments:  partialling _mismatch'
%     'Match Unclear: shared_segments partialling vowels'
%     'Match Unclear: vowels partialling shared_segments'
%     'Match Unclear: vowels partialling shared_segments_mismatch'
%     'Mismatch Clear shared_segments: _mismatch partialling '
%     'Mismatch Clear: shared_segments_mismatch partialling vowels'
%     'Mismatch Clear shared_segments:  partialling _mismatch'
%     'Mismatch Clear: shared_segments partialling vowels'
%     'Mismatch Clear: vowels partialling shared_segments'
%     'Mismatch Clear: vowels partialling shared_segments_mismatch'
%     'Mismatch Unclear shared_segments: _mismatch partialling '
%     'Mismatch Unclear: shared_segments_mismatch partialling vowels'
%     'Mismatch Unclear shared_segments:  partialling _mismatch'
%     'Mismatch Unclear: shared_segments partialling vowels'
%     'Mismatch Unclear: vowels partialling shared_segments'
%     'Mismatch Unclear: vowels partialling shared_segments_mismatch'
%     };

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
    'Mismatch Unclear to Written Cross-decode'
    'Mismatch Clear to Written Cross-decode'
    'Match Unclear to Written SS_Match'
    'Match Clear to Written SS_Match'
    'Mismatch Unclear to Written Shared Segments - no self'
    'Mismatch Clear to Written Shared Segments - no self'
%     'Match Unclear to Written SS_Match - no self'
%     'Match Clear to Written SS_Match - no self'
%     'Mismatch Unclear to Written SS_Match - no self'
%     'Mismatch Clear to Written SS_Match - no self'
    };

this_model_name{3} = {'All spoken Cross-decode_Match'
    'All spoken SS_Match'
    'All spoken SS_Match - no self'
    'Spoken to Written Cross-decode_Match'
    'Spoken to Written SS_Match - no self'
    'Spoken to Written Cross-decode_written'
    'Spoken to Written SS_written - no self'
    'Spoken to Written Cross-decode_written-lowpe'
    'Spoken to Written Cross-decode_written-highpe'
    'Match to Mismatch Shared Segments - no self'
    'Match to Mismatch SS_Match - no self'
    'Match to Mismatch combined_SS - no self - rescaled'
    'Match to Mismatch only cross'
    'Match to Mismatch only not cross'
    };

this_model_name{4} = {
    'Match Unclear to Match Clear Cross-decode_Match';
    'Mismatch Unclear to Mismatch Clear Cross-decode';
    'Match Unclear to Match Clear Cross-decode';
    'Mismatch Unclear to Mismatch Clear Cross-decode_Match';
    };

this_model_name{5} = {
    'M to MM Shared Segments:  Cross Negative partialling '
    'M to MM Shared Segments:  partialling  Cross Negative';
    };


% this_model_name{6} = {
%     'Match Unclear to Mismatch Unclear Cross-decode'
%     'Match Unclear to Mismatch Unclear Shared Segments - cross'
%     'Match Unclear to Mismatch Unclear Shared Segments - no self'
%     'Match Unclear to Mismatch Clear Cross-decode'
%     'Match Unclear to Mismatch Clear Shared Segments - cross'
%     'Match Unclear to Mismatch Clear Shared Segments - no self'
%     'Match Clear to Mismatch Unclear Cross-decode'
%     'Match Clear to Mismatch Unclear Shared Segments - cross'
%     'Match Clear to Mismatch Unclear Shared Segments - no self'
%     'Match Clear to Mismatch Clear Cross-decode'
%     'Match Clear to Mismatch Clear Shared Segments - cross'
%     'Match Clear to Mismatch Clear Shared Segments - no self'
%     'Match Unclear to Written Cross-decode'
%     'Match Clear to Written Cross-decode'
%     'Mismatch Unclear to Written Cross-decode'
%     'Mismatch Clear to Written Cross-decode'
%     };

nrun = size(subjects,2); % enter the number of runs here
% First load in the similarities
RSA_ROI_data_exist = zeros(1,nrun);
all_data = [];
mask_names{1} = {
    %     'rwLeft_IFG_cross_group_cluster'
    %     %     'rwLeft_Superior_Temporal_Gyrus';
    %     'rwL_STG_cross-segment_cluster'
    %     'rwBlank_2016_inflated'
    'rwLeft_Frontal_Univariate_MM>M'
    'rwLeft_Temporal_Univariate_MM>M'
    %     'rwLeft_PostSTG_Univariate_Interaction'
    %             'rwLeft_Precentral_Univariate_Interaction1'
    %             'rwLeft_Precentral_Univariate_Interaction2'
    %             'rwLeft_Precentral_Univariate_Interaction3'
    %     'rwLeft_Angular_Univariate_Interaction1'
    %     'rwLeft_Angular_Univariate_Interaction2'
    'rwLeft_Precentral_Univariate_Interaction_combined'
    %     'rwLeft_Angular_Univariate_Interaction_combined'
    %     'rwLeft_STG_Univariate8mm_15>3'
    'rwLeft_STG_Univariate3mm_15>3'
    'rwLeft_PrG_SSMatchnoself_combined'
    'rwLeft_PrG_All_Shared_Segments'
    };

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
all_rho = [];
all_corr_ps = [];
all_corrected_rho = [];
all_corrected_corr_ps = [];
for j = 1:length(this_model_name)
    for k = 1:length(mask_names)
        for i = 1:length(mask_names{k})
            all_data = [];
            all_corrected_data = [];
            for crun = 1:nrun
                %ROI_RSA_dir = [preprocessedpathstem subjects{crun} '/stats4_multi_3/TDTcrossnobis_ROI/RSA/spearman']; %Where are the results>
                ROI_RSA_dir = [preprocessedpathstem subjects{crun} '/stats4_multi_3/TDTcrossnobis_ROI/' mask_names{k}{i} '/RSA/spearman'];
                if ~exist(fullfile(ROI_RSA_dir,['roi_effects_' this_model_name{j}{1} '.mat']),'file')
                    ROI_RSA_dir = [preprocessedpathstem subjects{crun} '/stats4_multi_3/TDTcrossnobis_ROI' mask_names{k}{i} '/RSA/spearman']; % Stupid coding error earlier in analysis led to misnamed directories
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
            all_corrected_data(:,:,group==1) = es_removeBetween_rotated(all_data(:,:,group==1),[3,1,2]); %Subjects, conditions, measures columns = 3,1,2 here
            all_corrected_data(:,:,group==2) = es_removeBetween_rotated(all_data(:,:,group==2),[3,1,2]); %Subjects, conditions, measures columns = 3,1,2 here
            
            
            this_ROI = find(strcmp(mask_names{k}{i},roi_names));
            %Test covariates
            for m = 1:length(this_model_name{j})
                [all_rho(j,k,i,m,:),all_corr_ps(j,k,i,m,:)] = corr(covariates,squeeze(all_data(m,this_ROI,:)),'rows','pairwise');
                for this_corr = 1:size(all_corr_ps,5);
                if all_corr_ps(j,k,i,m,this_corr) < 0.05
                    disp(['Exploratory correlation in ' mask_names{k}{i}(3:end) ' ' this_model_name{j}{m} ' for ' covariate_names{this_corr}])
                end
                end
            end
            
            %this_ROI = find(strcmp('rwLeft_Superior_Temporal_Gyrus',roi_names));
            
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
            if verLessThan('matlab', '9.2')
                legend('Controls','Patients','location','southeast')
            else
                legend('Controls','Patients','location','southeast','AutoUpdate','off')
            end
            [h,p] = ttest(squeeze(all_data(:,this_ROI,logical(RSA_ROI_data_exist)))');
            these_y_lims = ylim;
            if sum(h)~=0
                plot(find(h),these_y_lims(2)-diff(these_y_lims/10),'g*')
            end
            [h,p] = ttest(squeeze(all_data(:,this_ROI,group==1&logical(RSA_ROI_data_exist)))');
            these_y_lims = ylim;
            if sum(h)~=0
                plot(find(h)-0.1,these_y_lims(2)-diff(these_y_lims/10),'k*')
            end
            [h,p] = ttest(squeeze(all_data(:,this_ROI,group==2&logical(RSA_ROI_data_exist)))');
            these_y_lims = ylim;
            if sum(h)~=0
                plot(find(h)+0.1,these_y_lims(2)-diff(these_y_lims/10),'r*')
            end
            
            [h,p] = ttest2(squeeze(all_data(:,this_ROI,group==1&logical(RSA_ROI_data_exist)))',squeeze(all_data(:,this_ROI,group==2&logical(RSA_ROI_data_exist)))');
            if sum(h)~=0
                plot(find(h),these_y_lims(2)-diff(these_y_lims/20),'gx')
            end
            for m = 1:length(this_model_name{j})
                for this_corr = 1:size(all_corr_ps,5);
                    if all_corr_ps(j,k,i,m,this_corr) < 0.05
                        text(m, these_y_lims(2)-(this_corr*diff(these_y_lims/100)),covariate_names{this_corr},'Interpreter','None')
                    end
                end
            end
            drawnow
           saveas(gcf,[outdir filesep mask_names{k}{i}(3:end) '_Model_set_' num2str(j) '.png'])
            
            for m = 1:length(this_model_name{j})
                [all_corrected_rho(j,k,i,m,:),all_corrected_corr_ps(j,k,i,m,:)] = corr(covariates,squeeze(all_data(m,this_ROI,:)),'rows','pairwise');
                for this_corr = 1:size(all_corr_ps,5);
                    if all_corr_ps(j,k,i,m,this_corr) < 0.05
                        disp(['Exploratory corrected correlation in ' mask_names{k}{i}(3:end) ' ' this_model_name{j}{m} ' for ' covariate_names{this_corr}])
                    end
                end
            end
            
            figure
            set(gcf,'Position',[100 100 1600 800]);
            set(gcf, 'PaperPositionMode', 'auto');
            hold on
            errorbar([1:length(this_model_name{j})]-0.1,nanmean(squeeze(all_corrected_data(:,this_ROI,group==1&RSA_ROI_data_exist)),2),nanstd(squeeze(all_corrected_data(:,this_ROI,group==1&RSA_ROI_data_exist))')/sqrt(sum(group==1&RSA_ROI_data_exist)),'kx')
            errorbar([1:length(this_model_name{j})]+0.1,nanmean(squeeze(all_corrected_data(:,this_ROI,group==2&RSA_ROI_data_exist)),2),nanstd(squeeze(all_corrected_data(:,this_ROI,group==2&RSA_ROI_data_exist))')/sqrt(sum(group==2&RSA_ROI_data_exist)),'rx')
            xlim([0 length(this_model_name{j})+1])
            set(gca,'xtick',[1:length(this_model_name{j})],'xticklabels',this_model_name{j},'XTickLabelRotation',45,'TickLabelInterpreter','none')
            plot([0 length(this_model_name{j})+1],[0,0],'k--')
            title(['Corrected ' mask_names{k}{i}(3:end) ' RSA'],'Interpreter','none')
            if verLessThan('matlab', '9.2')
                legend('Controls','Patients','location','southeast')
            else
                legend('Controls','Patients','location','southeast','AutoUpdate','off')
            end
            [h,p] = ttest(squeeze(all_corrected_data(:,this_ROI,logical(RSA_ROI_data_exist)))');
            these_y_lims = ylim;
            if sum(h)~=0
                plot(find(h),these_y_lims(2)-diff(these_y_lims/10),'g*')
            end
            [h,p] = ttest(squeeze(all_corrected_data(:,this_ROI,group==1&logical(RSA_ROI_data_exist)))');
            these_y_lims = ylim;
            if sum(h)~=0
                plot(find(h)-0.1,these_y_lims(2)-diff(these_y_lims/10),'k*')
            end
            [h,p] = ttest(squeeze(all_corrected_data(:,this_ROI,group==2&logical(RSA_ROI_data_exist)))');
            these_y_lims = ylim;
            if sum(h)~=0
                plot(find(h)+0.1,these_y_lims(2)-diff(these_y_lims/10),'r*')
            end
            
            [h,p] = ttest2(squeeze(all_corrected_data(:,this_ROI,group==1&logical(RSA_ROI_data_exist)))',squeeze(all_corrected_data(:,this_ROI,group==2&logical(RSA_ROI_data_exist)))');
            if sum(h)~=0
                plot(find(h),these_y_lims(2)-diff(these_y_lims/20),'gx')
            end
            for m = 1:length(this_model_name{j})
                for this_corr = 1:size(all_corrected_corr_ps,5);
                    if all_corrected_corr_ps(j,k,i,m,this_corr) < 0.05
                        text(m, these_y_lims(2)-(this_corr*diff(these_y_lims/100)),covariate_names{this_corr},'Interpreter','None')
                    end
                end
            end
            drawnow
            saveas(gcf,[outdir filesep 'Corrected_' mask_names{k}{i}(3:end) '_Model_set_' num2str(j) '.png'])
                       
        end
    end
end

