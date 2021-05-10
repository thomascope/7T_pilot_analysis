% %% Now do RSA on ROI data
% nrun = size(subjects,2); % enter the number of runs here
% RSAroiworkedcorrectly = zeros(1,nrun);
% masks = {
%     'rwBlank_2016_inflated'
%     'rwL_STG_cross-segment_cluster'
%     'rwLeft_Superior_Temporal_Gyrus'
%     'rwLeft_Angular_Gyrus'
%     'rwLeft_Precentral_Gyrus'
%     'rwLeft_Frontal_Operculum'
%     'rwLeft_Inferior_Frontal_Angular_Gyrus'
%     'rwRight_Superior_Temporal_Gyrus'
%     'rwRight_Angular_Gyrus'
%     'rwRight_Precentral_Gyrus'
%     'rwRight_Frontal_Operculum'
%     'rwRight_Inferior_Frontal_Angular_Gyrus'
%     'rwLeft_IFG_Written_Cluster'
%     'rwLeft_Precentral_Written_Cluster'
%     };
% parfor crun = 1:nrun
%     addpath(genpath('./RSA_scripts'))
%     GLMDir = [preprocessedpathstem subjects{crun} '/stats4_multi_3_nowritten2']; %Where is the SPM model?
%     try
%         module_roi_RSA(GLMDir,masks)
%         RSAroiworkedcorrectly(crun) = 1;
%     catch
%         RSAroiworkedcorrectly(crun) = 0;
%     end
% end


for j = 1:length(this_model_name)
    for k = 1:length(mask_names)
        for i = 1:length(mask_names{k})
            all_data = [];
            for crun = 1:nrun
                %ROI_RSA_dir = [preprocessedpathstem subjects{crun} '/stats4_multi_3_nowritten2/TDTcrossnobis_ROI/RSA/spearman']; %Where are the results>
                ROI_RSA_dir = [preprocessedpathstem subjects{crun} '/stats4_multi_3_nowritten2/TDTcrossnobis_ROI/' mask_names{k}{i} '/RSA/spearman'];
                if ~exist(fullfile(ROI_RSA_dir,['roi_effects_' this_model_name{j}{1} '.mat']),'file')
                    ROI_RSA_dir = [preprocessedpathstem subjects{crun} '/stats4_multi_3_nowritten2/TDTcrossnobis_ROI' mask_names{k}{i} '/RSA/spearman']; % Stupid coding error earlier in analysis led to misnamed directories
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
            if verLessThan('matlab', '9.2')
                legend('Controls','Patients','location','southeast')
            else
                legend('Controls','Patients','location','southeast','AutoUpdate','off')
            end
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



for j = 1:length(this_model_name)
    for k = 1:length(mask_names)
        for i = 1:length(mask_names{k})
            all_data = [];
            for crun = 1:nrun
                %ROI_RSA_dir = [preprocessedpathstem subjects{crun} '/stats4_multi_3_noabsent/TDTcrossnobis_ROI/RSA/spearman']; %Where are the results>
                ROI_RSA_dir = [preprocessedpathstem subjects{crun} '/stats4_multi_3_noabsent/TDTcrossnobis_ROI/' mask_names{k}{i} '/RSA/spearman'];
                if ~exist(fullfile(ROI_RSA_dir,['roi_effects_' this_model_name{j}{1} '.mat']),'file')
                    ROI_RSA_dir = [preprocessedpathstem subjects{crun} '/stats4_multi_3_noabsent/TDTcrossnobis_ROI' mask_names{k}{i} '/RSA/spearman']; % Stupid coding error earlier in analysis led to misnamed directories
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
            if verLessThan('matlab', '9.2')
                legend('Controls','Patients','location','southeast')
            else
                legend('Controls','Patients','location','southeast','AutoUpdate','off')
            end
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



for j = 1:length(this_model_name)
    for k = 1:length(mask_names)
        for i = 1:length(mask_names{k})
            all_data = [];
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
            if verLessThan('matlab', '9.2')
                legend('Controls','Patients','location','southeast')
            else
                legend('Controls','Patients','location','southeast','AutoUpdate','off')
            end
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