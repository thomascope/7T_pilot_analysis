function [all_sigma_pred,all_thresholds,controls_sigma_pred,controls_threshold,patients_sigma_pred,patients_threshold] = module_bayesian_behaviour(subjects,group,dates)
% A script for doing Bayesian behavioural modelling on PINFA data before
% and after scan. Expects you to have run the extraction scripts for the
% clarity rating and vocode report tasks beforehand, and saved them in a
% 'beh' structure (see analysebehaviourneutral_PINFA.m and
% Analyse_Vocode_Report_VESPAPre.m)
% Outputs two matrices for the standard deviation of the prior and the
% perceptual threshold, each containing a column per timepoint and row per
% subject

addpath('./plotting')

global meansarray report_data optimise_this_subject t % Inputs to search subfunction
global normalised_model_meansarray % Output from search subfunction
%global mu_input_mismatch mu_pred

% First read in data
clarity_data_dir = ['./NonScanTasks/Clarity rating/data/'];
report_data_dir = ['./NonScanTasks/Vocode report/data/'];
timepoints = {'Pre','Post'};
chance_report = 25; % Chance performance on the vocode report task - 4AFC

% For this dataset we have a problem, which is that vocode report was done
% with 3, 6 and 15 vocoder channels and clarity rating with 4, 8 and 16
% channels. This results in mean report ratios across all participants of
% 1:1.0736:1.3084 Pre, and 1:1.0313:1.2497 Post, while in the Nature Comms 2017 paper the mean report
% ratios were 1:1.2572:1.4031. I have therefore boosted the report 6 by a
% factor of 1.0905 to make the difficulty ratios from low-med-high more consistent.
report_ratio_factor = [1,1.0905,1];

for i = 1:length(subjects)
    for t = 1:length(timepoints)
        % First read in Vocode Report Performance
        temp_data_path = dir([report_data_dir 'beh_PA*' subjects{i} '*' timepoints{t} '_' dates{i}(3:end) '.mat']);
        if isempty(temp_data_path)
            disp(['Missing Vocode Report data for subject ' subjects{i} ' at timepoint ' timepoints{t} '. Please check.'])
            report_data(t,i,:) = nan(1,3);
        else
            assert(length(temp_data_path)==1,['More than one Vocode Report file found for subject ' subjects{i} ' at timepoint ' timepoints{t} '. Aborting.'])
            these_report_data = load([report_data_dir temp_data_path(1).name]);
            if group(i) == 1
                these_report_data = these_report_data.Condat; %Sorry about the capitalisation
            else
                these_report_data = these_report_data.patdat;
            end
            report_data(t,i,:) = mean(these_report_data.meansarray').*report_ratio_factor; %Average across distractor types, end up with 3-element vector per participant for low, med, high vocoder detail
        end
        
        % Then read in clarity ratings
        temp_data_path = dir([clarity_data_dir 'beh*' subjects{i} '*' timepoints{t} '.mat']);
        if isempty(temp_data_path)
            disp(['Missing Clarity Rating data for subject ' subjects{i} ' at timepoint ' timepoints{t} '. Please check.'])
            clarity_data(t,i,:) = nan(1,9);
        else
            assert(length(temp_data_path)==1,['More than one Vocode Report file found for subject ' subjects{i} ' at timepoint ' timepoints{t} '. Aborting.'])
            these_clarity_data = load([clarity_data_dir temp_data_path(1).name]);
            if group(i) == 1
                these_clarity_data = these_clarity_data.condat;
            else
                these_clarity_data = these_clarity_data.patdat;
            end
            % Structure of expected meansarray = [mismatch_rating_4_chan_patient(patnum), mismatch_rating_8_chan_patient(patnum), mismatch_rating_16_chan_patient(patnum), neutral_rating_4_chan_patient(patnum), neutral_rating_8_chan_patient(patnum), neutral_rating_16_chan_patient(patnum), match_rating_4_chan_patient(patnum), match_rating_8_chan_patient(patnum), match_rating_16_chan_patient(patnum)];
            clarity_data(t,i,:) = [mean(these_clarity_data.mismatch4), mean(these_clarity_data.mismatch8), mean(these_clarity_data.mismatch16), mean(these_clarity_data.neutral4), mean(these_clarity_data.neutral8), mean(these_clarity_data.neutral16),mean(these_clarity_data.match4), mean(these_clarity_data.match8), mean(these_clarity_data.match16)];
        end
    end
end

% % Work out the priors - ?Needed?
% mu_pred = 10; %Arbitrarily place prediction and input on x-scale. When 'match', mu_input will be the same, when mismatch, mu_input will be different
%         mu_input_match = 10; % Location of the prior and posterior in Bayesian space, total mismatch so far apart
%         mu_input_mismatch = 20;

connum = 0;
patnum = 0;
for optimise_this_subject = 1:length(subjects)
    if group(optimise_this_subject) == 1
        connum = connum+1;
    else
        patnum = patnum+1;
    end
    
    for t = 1:length(timepoints)
        meansarray = squeeze(clarity_data(t,optimise_this_subject,:))';
        denormed_meansarray = (meansarray-1)./3; %De-normalise
        
        if any(isnan(meansarray)) || any(isnan(squeeze(report_data(t,optimise_this_subject,:))))
            all_sigma_pred(t,optimise_this_subject) = NaN;
            all_thresholds(t,optimise_this_subject) = NaN;
            if group(optimise_this_subject) == 1
                all_meansarray_controls(t,connum,:,:) = nan(1,18);
                
                controls_sigma_pred(t,connum) = NaN;
                controls_threshold(t,connum) = NaN;
            else
                all_meansarray_patients(t,patnum,:,:) = nan(1,18);
                
                patients_sigma_pred(t,patnum) = NaN;
                patients_threshold(t,patnum) = NaN;
                
            end
            continue
            
        else
            [model_arguments, fval] = patternsearch(@search_subfunction_bayes_modelling_PINFA,[3,0],[],[],[],[],[eps, 0],[3, 1]);
            sigma_pred = model_arguments(1);
            threshold = model_arguments(2);
            
            eval(['modelled_error_for_subject' num2str(optimise_this_subject) ' = search_subfunction_bayes_modelling_PINFA(model_arguments)'])
            
            %             grouped_meansarray = [meansarray([1,7]);meansarray([2,8]);meansarray([3,9])];
            %             normalised_grouped_meansarray = [normalised_model_meansarray([1,7]);normalised_model_meansarray([2,8]);normalised_model_meansarray([3,9])];
            %
            %             all_meansarray = [normalised_grouped_meansarray,grouped_meansarray];
            all_meansarray = [normalised_model_meansarray,meansarray];
            
            stesarray = zeros(size(all_meansarray)); %Dummy stes for now.
            %             figure
            %             set(gcf,'position',[100,100,1200,800])
            %             barweb(all_meansarray,stesarray,[],{'4 channels';'8 channels';'16 channels'},['Compared Clarity Ratings by Prime Type and Vocoder Channels for Subject ' subjects{optimise_this_subject}],[],'Mean Clarity Rating',[],[],{'Model_Match','Model_Mismatch','Actual_Match','Actual_Mismatch'}) ;
            %             legend('Model Match','Model Mismatch','Actual Match','Actual Mismatch','location','NorthWest');
            %             set(gca,'ylim',[1,4]);
            
            all_sigma_pred(t,optimise_this_subject) = sigma_pred;
            all_thresholds(t,optimise_this_subject) = threshold;
            
            if group(optimise_this_subject) == 1
                all_meansarray_controls(t,connum,:,:) = all_meansarray;
                
                controls_sigma_pred(t,connum) = sigma_pred;
                controls_threshold(t,connum) = threshold;
            else
                all_meansarray_patients(t,patnum,:,:) = all_meansarray;
                
                patients_sigma_pred(t,patnum) = sigma_pred;
                patients_threshold(t,patnum) = threshold;
                
            end
        end
    end
end

for t = 1:length(timepoints)
    figure
    stesarray = reshape(stesarray,3,6);
    barweb(reshape(nanmean(squeeze(all_meansarray_patients(t,:,:)),1),3,6),stesarray,[],{'4 channels';'8 channels';'16 channels'},['Compared Clarity Ratings by Prime Type and Vocoder Channels for All Patients'],[],'Mean Clarity Rating',[],[],{}) ;
    set(gca,'ylim',[1,4]);
    figure
    barweb(reshape(nanmean(squeeze(all_meansarray_controls(t,:,:)),1),3,6),stesarray,[],{'4 channels';'8 channels';'16 channels'},['Compared Clarity Ratings by Prime Type and Vocoder Channels for All Controls'],[],'Mean Clarity Rating',[],[],{}) ;
    set(gca,'ylim',[1,4]);
    
    patient_means_forline = reshape(nanmean(squeeze(all_meansarray_patients(t,:,:)),1),3,6);
    control_means_forline = reshape(nanmean(squeeze(all_meansarray_controls(t,:,:)),1),3,6);
    
    figure
    set(gcf,'position',[100,100,1200,800])
    lineplot = tight_subplot(1,2,[0 0],[.1 .1],[.1 .1]);
    axes(lineplot(1));
    set(gca,'ylim',[1,4])
    errorbar(control_means_forline(:,1),stesarray(:,1),'r--','linewidth',3);
    hold on
    errorbar(control_means_forline(:,3),stesarray(:,3),'k--','linewidth',3);
    errorbar(control_means_forline(:,4),stesarray(:,4),'r','linewidth',3);
    %errorbar(control_means_forline(:,5),stesarray(:,5),'g','linewidth',3);
    errorbar(control_means_forline(:,6),stesarray(:,6),'k','linewidth',3);
    set(gca,'ylim',[1,4],'LineWidth', 2, 'Xtick', [1 2 3], 'XTickLabel',[4,8,16],'Fontsize',[14],'FontName','Tahoma','YAxisLocation','right')
    secondlegend = legend('Model MisMatch','Model Match','MisMatch','Match','location','SouthEast');
    set(secondlegend,'FontSize',18);
    title(['Controls ' timepoints{t}],'Color','k','fontsize',20)
    ylabel('Clarity Rating')
    xlabel('Vocode Channels')
    
    axes(lineplot(2));
    set(gca,'ylim',[1,4])
    errorbar(patient_means_forline(:,1),stesarray(:,1),'r--','linewidth',3);
    hold on
    errorbar(patient_means_forline(:,3),stesarray(:,3),'k--','linewidth',3);
    errorbar(patient_means_forline(:,4),stesarray(:,4),'r','linewidth',3);
    errorbar(patient_means_forline(:,6),stesarray(:,6),'k','linewidth',3);
    set(gca,'ylim',[1,4],'LineWidth', 2, 'Xtick', [1 2 3], 'XTickLabel',[4,8,16],'Fontsize',[14],'FontName','Tahoma','YAxisLocation','right')
    set(secondlegend,'FontSize',18);
    title(['Patients ' timepoints{t}],'Color','k','fontsize',20)
    ylabel('Clarity Rating')
    xlabel('Vocode Channels')
    img = getframe(gcf);
    try %This is an absolute path, so clearly won't work anywhere except my laptop
        addpath('C:\Users\Thoma\Documents\Academic work\MATLAB\ojwoodford-export_fig-216b30e')
        eval(['export_fig Model_Fits_NoNeutral_'  timepoints{t} '.png -transparent'])
        eval(['export_fig Model_Fits_NoNeutral_'  timepoints{t} '.pdf -transparent'])
    catch
    end
    
    
    figure
    boxplot(all_sigma_pred(t,:),group,'plotstyle','compact','symbol','k.','medianstyle','line','colors',hsv2rgb([0 0.6 0.6; 0.3 0.6 0.6; 0.6 0.6 0.6]),'jitter',0,'widths',0.5)
    title('Standard Deviation of Prior')
    xtix = {'\bfControls','\bfnfvPPA'};   % Your labels
    xtixloc = [1 2];      % Your label locations
    %     set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc,'FontSize',15,'ylim',[0 3]);
    set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc,'FontSize',15);
    set(findobj(gca,'Type','text'),'FontSize',80);
    set(findobj(gcf,'Tag','Median'),'Color',[0 0 0],'LineWidth',4);
    set(findobj(gcf,'Tag','Outliers'),'MarkerEdgeColor',[0 0 0],'MarkerSize',30);
    set(findobj(gcf,'Tag','Box'),'linewidth',20);
    set(findobj(gcf,'Tag','Whisker'),'linewidth',4);
    ylabel('\bfA.U.')
    drawnow
    try %This is an absolute path, so clearly won't work anywhere except my laptop
        tc_sigstar([1,2],ranksum(patients_sigma_pred,controls_sigma_pred))
        eval(['export_fig Standard_Deviation_of_Prior'  timepoints{t} '.png -transparent'])
        eval(['export_fig Standard_Deviation_of_Prior'  timepoints{t} '.pdf -transparent'])
    catch
    end
    
    figure
    boxplot(all_thresholds(t,:),group,'plotstyle','compact','symbol','k.','medianstyle','line','colors',hsv2rgb([0 0.6 0.6; 0.3 0.6 0.6; 0.6 0.6 0.6]),'jitter',0,'widths',0.5)
    title('Perceptual Threshold')
    xtix = {'\bfControls','\bfnfvPPA'};   % Your labels
    xtixloc = [1 2];      % Your label locations
    %     set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc,'FontSize',15,'ylim',[0 0.4]);
    set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc,'FontSize',15);
    set(findobj(gca,'Type','text'),'FontSize',80);
    set(findobj(gcf,'Tag','Median'),'Color',[0 0 0],'LineWidth',4);
    set(findobj(gcf,'Tag','Outliers'),'MarkerEdgeColor',[0 0 0],'MarkerSize',30);
    set(findobj(gcf,'Tag','Box'),'linewidth',20);
    set(findobj(gcf,'Tag','Whisker'),'linewidth',4);
    ylabel('\bfA.U.')
    drawnow
    try %This is an absolute path, so clearly won't work anywhere except my laptop
        eval(['export_fig Perceptual_Threshold'  timepoints{t} '.png -transparent'])
        eval(['export_fig Perceptual_Threshold'  timepoints{t} '.pdf -transparent'])
    catch
    end
    
    nonpara_pval_sigma_pred(t) = ranksum(patients_sigma_pred(t,:),controls_sigma_pred(t,:))
    nonpara_pval_threshold(t) = ranksum(patients_threshold(t,:),controls_threshold(t,:))
    nonpara_pval_ratios(t) = ranksum(patients_sigma_pred(t,:)./patients_threshold(t,:),controls_sigma_pred(t,:)./controls_threshold(t,:))
    
    [~, para_pval_sigma_pred(t)] = ttest2(patients_sigma_pred(t,:),controls_sigma_pred(t,:),'vartype','unequal')
    [~, para_pval_threshold(t)] = ttest2(patients_threshold(t,:),controls_threshold(t,:),'vartype','unequal')
    
end
