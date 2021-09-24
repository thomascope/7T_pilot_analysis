function plot_behavioural_data(subjects, dates, group, graph_this)

cd('./behavioural_data')
control_response_averages = [];
control_rt_averages = [];
control_rt_medians = [];
control_AFCs = [];
control_running_average = [];
patient_response_averages = [];
patient_rt_medians = [];
patient_rt_averages = [];
patient_rt_medians = [];
patient_AFCs = [];
patient_running_average = [];
all_response_averages = [];
all_rt_averages = []; 
all_rt_medians = [];
AFCs = []; 
running_average = [];
control_running_average_bycond = [];
patient_running_average_bycond = [];
control_normalised_running_average = [];
patient_normalised_running_average = [];

for crun = 1:length(subjects)
    
    [all_response_averages(crun,:), all_rt_averages(crun,:), all_rt_medians(crun,:), AFCs(crun,:), running_average(crun,:), running_average_bycond(crun,:,:), normalised_running_average(crun,:), response_order(crun)] = AFC_graph_this_subject_2021(subjects{crun}, dates{crun}, graph_this);
    
    if nansum(nansum(all_response_averages(crun,:))) == 0 
        continue % Ignore the patient who did not press any buttons
    end
    if group(crun) == 1 % Controls
        control_response_averages(end+1,:) = all_response_averages(crun,:);
        control_rt_averages(end+1,:) = all_rt_averages(crun,:);
        control_rt_medians(end+1,:) = all_rt_medians(crun,:);
        control_AFCs(end+1,:) = AFCs(crun,:);
        control_running_average(end+1,:) = running_average(crun,:);
        control_running_average_bycond(end+1,:,:) = running_average_bycond(crun,:,:);
        control_normalised_running_average(end+1,:) = normalised_running_average(crun,:);
    elseif group(crun) == 2 % Patients
        patient_response_averages(end+1,:) = all_response_averages(crun,:);
        patient_rt_averages(end+1,:) = all_rt_averages(crun,:);
        patient_rt_medians(end+1,:) = all_rt_medians(crun,:);
        patient_AFCs(end+1,:) = AFCs(crun,:);
        patient_running_average(end+1,:) = running_average(crun,:);
        patient_running_average_bycond(end+1,:,:) = running_average_bycond(crun,:,:);
        patient_normalised_running_average(end+1,:) = normalised_running_average(crun,:);
    end
end
cd('../')

%% First plot perceptual/task learning
figure
set(gcf,'Position',[100 100 1600 800]);
subplot(2,1,1)
title('5 Trial Moving Average')
hold on
plot(1:size(patient_running_average,2),nanmean(patient_running_average(patient_AFCs==2,:),1),'k--');
plot(1:size(patient_running_average,2),nanmean(patient_running_average(patient_AFCs==4,:),1),'r--');
plot(1:size(control_running_average,2),nanmean(control_running_average(control_AFCs==2,:),1),'k-');
plot(1:size(control_running_average,2),nanmean(control_running_average(control_AFCs==4,:),1),'r-');
legend({'Patient 2AFC','Patient 4AFC','Control 2AFC','Control 4AFC'},'location','SouthEast')
xlabel('Trial Number')
ylabel('Percent Correct')

subplot(2,1,2)
hold on
LM{1} = fitlm(1:size(patient_running_average,2),nanmean(patient_running_average(patient_AFCs==2,:),1));
all_r_squareds(1) = LM{1}.Rsquared.Ordinary;
all_model_ps(1) = LM{1}.anova.pValue(1);
this_reg_line = plot(LM{1});
this_reg_line(1).MarkerEdgeColor = 'k';
this_reg_line(1).Marker = 'none';
this_reg_line(2).LineStyle = '--';
this_reg_line(2).Color = 'k';
this_reg_line(3).Color = 'k';
this_reg_line(4).Color = 'k';
LM{2} = fitlm(1:size(patient_running_average,2),nanmean(patient_running_average(patient_AFCs==4,:),1));
all_r_squareds(2) = LM{2}.Rsquared.Ordinary;
all_model_ps(2) = LM{2}.anova.pValue(1);
this_reg_line = plot(LM{2});
this_reg_line(1).MarkerEdgeColor = 'r';
this_reg_line(1).Marker = 'none';
this_reg_line(2).LineStyle = '--';
this_reg_line(2).Color = 'r';
this_reg_line(3).Color = 'r';
this_reg_line(4).Color = 'r';
LM{3} = fitlm(1:size(control_running_average,2),nanmean(control_running_average(control_AFCs==2,:),1));
all_r_squareds(3) = LM{3}.Rsquared.Ordinary;
all_model_ps(3) = LM{3}.anova.pValue(1);
this_reg_line = plot(LM{3});
this_reg_line(1).MarkerEdgeColor = 'k';
this_reg_line(1).Marker = 'none';
this_reg_line(2).LineStyle = '-';
this_reg_line(2).Color = 'k';
this_reg_line(3).Color = 'k';
this_reg_line(4).Color = 'k';
LM{4} = fitlm(1:size(control_running_average,2),nanmean(control_running_average(control_AFCs==4,:),1));
all_r_squareds(4) = LM{4}.Rsquared.Ordinary;
all_model_ps(4) = LM{4}.anova.pValue(1);
this_reg_line = plot(LM{4});
this_reg_line(1).MarkerEdgeColor = 'r';
this_reg_line(1).Marker = 'none';
this_reg_line(2).LineStyle = '-';
this_reg_line(2).Color = 'r';
this_reg_line(3).Color = 'r';
this_reg_line(4).Color = 'r';
title('Linear Model Fits')
legend('off')
xlabel('Trial Number')
ylabel('Percent Correct')
suptitle('Running Average Performance')
for i = 1:4
    text(65,predict(LM{i},66),['p=' num2str(all_model_ps(i),2) ', r2=' num2str(all_r_squareds(i),2)])
end

figure
set(gcf,'Position',[100 100 1600 800]);
subplot(2,1,1)
title('5 Trial Moving Average')
hold on
plot(1:size(patient_normalised_running_average,2),nanmean(patient_normalised_running_average(patient_AFCs==2,:),1),'k--');
plot(1:size(patient_normalised_running_average,2),nanmean(patient_normalised_running_average(patient_AFCs==4,:),1),'r--');
plot(1:size(control_normalised_running_average,2),nanmean(control_normalised_running_average(control_AFCs==2,:),1),'k-');
plot(1:size(control_normalised_running_average,2),nanmean(control_normalised_running_average(control_AFCs==4,:),1),'r-');
legend({'Patient 2AFC','Patient 4AFC','Control 2AFC','Control 4AFC'},'location','SouthEast')
xlabel('Trial Number')
ylabel('Percent Correct')

subplot(2,1,2)
hold on
LM{1} = fitlm(1:size(patient_normalised_running_average,2),nanmean(patient_normalised_running_average(patient_AFCs==2,:),1));
all_r_squareds(1) = LM{1}.Rsquared.Ordinary;
all_model_ps(1) = LM{1}.anova.pValue(1);
this_reg_line = plot(LM{1});
this_reg_line(1).MarkerEdgeColor = 'k';
this_reg_line(1).Marker = 'none';
this_reg_line(2).LineStyle = '--';
this_reg_line(2).Color = 'k';
this_reg_line(3).Color = 'k';
this_reg_line(4).Color = 'k';
LM{2} = fitlm(1:size(patient_normalised_running_average,2),nanmean(patient_normalised_running_average(patient_AFCs==4,:),1));
all_r_squareds(2) = LM{2}.Rsquared.Ordinary;
all_model_ps(2) = LM{2}.anova.pValue(1);
this_reg_line = plot(LM{2});
this_reg_line(1).MarkerEdgeColor = 'r';
this_reg_line(1).Marker = 'none';
this_reg_line(2).LineStyle = '--';
this_reg_line(2).Color = 'r';
this_reg_line(3).Color = 'r';
this_reg_line(4).Color = 'r';
LM{3} = fitlm(1:size(control_normalised_running_average,2),nanmean(control_normalised_running_average(control_AFCs==2,:),1));
all_r_squareds(3) = LM{3}.Rsquared.Ordinary;
all_model_ps(3) = LM{3}.anova.pValue(1);
this_reg_line = plot(LM{3});
this_reg_line(1).MarkerEdgeColor = 'k';
this_reg_line(1).Marker = 'none';
this_reg_line(2).LineStyle = '-';
this_reg_line(2).Color = 'k';
this_reg_line(3).Color = 'k';
this_reg_line(4).Color = 'k';
LM{4} = fitlm(1:size(control_normalised_running_average,2),nanmean(control_normalised_running_average(control_AFCs==4,:),1));
all_r_squareds(4) = LM{4}.Rsquared.Ordinary;
all_model_ps(4) = LM{4}.anova.pValue(1);
this_reg_line = plot(LM{4});
this_reg_line(1).MarkerEdgeColor = 'r';
this_reg_line(1).Marker = 'none';
this_reg_line(2).LineStyle = '-';
this_reg_line(2).Color = 'r';
this_reg_line(3).Color = 'r';
this_reg_line(4).Color = 'r';
title('Linear Model Fits')
legend('off')
xlabel('Trial Number')
ylabel('Percent Correct')
suptitle('Normalised Running Average Performance')
for i = 1:4
    text(65,predict(LM{i},66),['p=' num2str(all_model_ps(i),2) ', r2=' num2str(all_r_squareds(i),2)])
end

condition_order = {'Match 3','Mismatch 3','Match 15','Mismatch 15'};
for j = 1:length(condition_order)
    figure
    set(gcf,'Position',[100 100 1600 800]);
    subplot(2,1,1)
    title('5 Trial Moving Average')
    hold on
    plot(1:size(patient_running_average_bycond,2),nanmean(patient_running_average_bycond(patient_AFCs==2,:,j),1),'k--');
    plot(1:size(patient_running_average_bycond,2),nanmean(patient_running_average_bycond(patient_AFCs==4,:,j),1),'r--');
    plot(1:size(control_running_average_bycond,2),nanmean(control_running_average_bycond(control_AFCs==2,:,j),1),'k-');
    plot(1:size(control_running_average_bycond,2),nanmean(control_running_average_bycond(control_AFCs==4,:,j),1),'r-');
    legend({'Patient 2AFC','Patient 4AFC','Control 2AFC','Control 4AFC'},'location','SouthEast')
    xlabel('Trial Number')
    ylabel('Percent Correct')
    
    subplot(2,1,2)
    hold on
    LM{1} = fitlm(1:size(patient_running_average_bycond,2),nanmean(patient_running_average_bycond(patient_AFCs==2,:,j),1));
    all_r_squareds(1) = LM{1}.Rsquared.Ordinary;
    all_model_ps(1) = LM{1}.anova.pValue(1);
    this_reg_line = plot(LM{1});
    this_reg_line(1).MarkerEdgeColor = 'k';
    this_reg_line(1).Marker = 'none';
    this_reg_line(2).LineStyle = '--';
    this_reg_line(2).Color = 'k';
    this_reg_line(3).Color = 'k';
    this_reg_line(4).Color = 'k';
    LM{2} = fitlm(1:size(patient_running_average_bycond,2),nanmean(patient_running_average_bycond(patient_AFCs==4,:,j),1));
    all_r_squareds(2) = LM{2}.Rsquared.Ordinary;
    all_model_ps(2) = LM{2}.anova.pValue(1);
    this_reg_line = plot(LM{2});
    this_reg_line(1).MarkerEdgeColor = 'r';
    this_reg_line(1).Marker = 'none';
    this_reg_line(2).LineStyle = '--';
    this_reg_line(2).Color = 'r';
    this_reg_line(3).Color = 'r';
    this_reg_line(4).Color = 'r';
    LM{3} = fitlm(1:size(control_running_average_bycond,2),nanmean(control_running_average_bycond(control_AFCs==2,:,j),1));
    all_r_squareds(3) = LM{3}.Rsquared.Ordinary;
    all_model_ps(3) = LM{3}.anova.pValue(1);
    this_reg_line = plot(LM{3});
    this_reg_line(1).MarkerEdgeColor = 'k';
    this_reg_line(1).Marker = 'none';
    this_reg_line(2).LineStyle = '-';
    this_reg_line(2).Color = 'k';
    this_reg_line(3).Color = 'k';
    this_reg_line(4).Color = 'k';
    LM{4} = fitlm(1:size(control_running_average_bycond,2),nanmean(control_running_average_bycond(control_AFCs==4,:,j),1));
    all_r_squareds(4) = LM{4}.Rsquared.Ordinary;
    all_model_ps(4) = LM{4}.anova.pValue(1);
    this_reg_line = plot(LM{4});
    this_reg_line(1).MarkerEdgeColor = 'r';
    this_reg_line(1).Marker = 'none';
    this_reg_line(2).LineStyle = '-';
    this_reg_line(2).Color = 'r';
    this_reg_line(3).Color = 'r';
    this_reg_line(4).Color = 'r';
    title('Linear Model Fits')
    legend('off')
    xlabel('Trial Number')
    ylabel('Percent Correct')
    suptitle(['Running Average Performance for condition ' condition_order{j}])
    for i = 1:4
        text(15,predict(LM{i},15),['p=' num2str(all_model_ps(i),2) ', r2=' num2str(all_r_squareds(i),2)])
    end
end

%Linear mixed effects model
crun = 1;
lme_table = table(crun*ones(size(running_average,2),1),group(crun)*ones(size(running_average,2),1),AFCs(crun)*ones(size(running_average,2),1),response_order(crun).this_cue_types',response_order(crun).this_vocoder_channels',[1:size(running_average,2)]',running_average(crun,:)','VariableNames',{'Subject','Diagnosis','AFCs','Congruency','Clarity','Trial','Running_Average'});

for crun = 2:length(subjects)
    this_subj_table = table(crun*ones(size(running_average,2),1),group(crun)*ones(size(running_average,2),1),AFCs(crun)*ones(size(running_average,2),1),response_order(crun).this_cue_types',response_order(crun).this_vocoder_channels',[1:size(running_average,2)]',running_average(crun,:)','VariableNames',{'Subject','Diagnosis','AFCs','Congruency','Clarity','Trial','Running_Average'});
    lme_table = [lme_table;this_subj_table];
end
lme = fitlme(lme_table,'Running_Average ~ Diagnosis+AFCs+Congruency+Clarity+Clarity:Congruency+Trial+(Trial|Subject)');

%% Now plot all individual results
figure
%cmap = colormap(parula(2));
colormap(jet)
set(gcf,'Position',[100 100 1600 800]);
subplot(3,1,1)
hold on
for this_x = 1:size(patient_response_averages,2)
    scatter(repmat(this_x+0.1,size(patient_response_averages,1),1)+rand(size(patient_response_averages,1),1)/30-(1/60),patient_response_averages(:,this_x),16,patient_AFCs/2)
    errorbar(this_x+0.2,mean(patient_response_averages(patient_AFCs==2,this_x)),std(patient_response_averages(patient_AFCs==2,this_x))/sqrt(sum(patient_AFCs==2)),'kx')
    errorbar(this_x+0.2,mean(patient_response_averages(patient_AFCs==4,this_x)),std(patient_response_averages(patient_AFCs==4,this_x))/sqrt(sum(patient_AFCs==4)),'rx')
    scatter(repmat(this_x-0.1,size(control_response_averages,1),1)+rand(size(control_response_averages,1),1)/30-(1/60),control_response_averages(:,this_x),16,control_AFCs/2)
    errorbar(this_x-0.2,mean(control_response_averages(control_AFCs==2,this_x)),std(control_response_averages(control_AFCs==2,this_x))/sqrt(sum(control_AFCs==2)),'kx')
    errorbar(this_x-0.2,mean(control_response_averages(control_AFCs==4,this_x)),std(control_response_averages(control_AFCs==4,this_x))/sqrt(sum(control_AFCs==4)),'rx')
end
xlim([0 6])
ylim([0 100])
title('Percent Correct')
set(gca,'XTick',[0:1:6])
set(gca,'XTickLabel',{'','Match 3','Mismatch 3','','Match 15','Mismatch 15',''},'XTickLabelRotation',15)

RM_table = [table(group',AFCs,'VariableNames',{'Diagnosis','AFCs'}),array2table(all_response_averages(:,[1:2,4:5]))];
factorNames = {'Congruency','Clarity'};
all_congruencies = {'Match','Mismatch','Match','Mismatch'};
all_clarities = {'3','3','15','15'};
withindesign = table(all_congruencies',all_clarities','VariableNames',factorNames);
rm = fitrm(RM_table,'Var1-Var4~Diagnosis+AFCs','WithinDesign',withindesign);
accuracy_ranovatbl = ranova(rm, 'WithinModel','Congruency*Clarity');

subplot(3,1,2)
hold on
for this_x = 1:size(patient_response_averages,2)
    scatter(repmat(this_x+0.1,size(patient_rt_averages,1),1)+rand(size(patient_rt_averages,1),1)/30-(1/60),patient_rt_averages(:,this_x),16,patient_AFCs/2+1)
        errorbar(this_x+0.2,mean(patient_rt_averages(patient_AFCs==2,this_x)),std(patient_rt_averages(patient_AFCs==2,this_x))/sqrt(sum(patient_AFCs==2)),'kx')
    errorbar(this_x+0.2,mean(patient_rt_averages(patient_AFCs==4,this_x)),std(patient_rt_averages(patient_AFCs==4,this_x))/sqrt(sum(patient_AFCs==4)),'rx')
    scatter(repmat(this_x-0.1,size(control_rt_averages,1),1)+rand(size(control_rt_averages,1),1)/30-(1/60),control_rt_averages(:,this_x),16,control_AFCs/2+1)
        errorbar(this_x-0.2,mean(control_rt_averages(control_AFCs==2,this_x)),std(control_rt_averages(control_AFCs==2,this_x))/sqrt(sum(control_AFCs==2)),'kx')
    errorbar(this_x-0.2,mean(control_rt_averages(control_AFCs==4,this_x)),std(control_rt_averages(control_AFCs==4,this_x))/sqrt(sum(control_AFCs==4)),'rx')
end
xlim([0 6])
title('Mean RT')
set(gca,'XTick',[0:1:6])
set(gca,'XTickLabel',{'','Match 3','Mismatch 3','','Match 15','Mismatch 15',''},'XTickLabelRotation',15)

RM_table = [table(group',AFCs,'VariableNames',{'Diagnosis','AFCs'}),array2table(all_rt_averages(:,[1:2,4:5]))];
factorNames = {'Congruency','Clarity'};
all_congruencies = {'Match','Mismatch','Match','Mismatch'};
all_clarities = {'3','3','15','15'};
withindesign = table(all_congruencies',all_clarities','VariableNames',factorNames);
rm = fitrm(RM_table,'Var1-Var4~Diagnosis+AFCs','WithinDesign',withindesign);
rt_ranovatbl = ranova(rm, 'WithinModel','Congruency*Clarity');

subplot(3,1,3)
hold on
for this_x = 1:size(patient_rt_medians,2)
    scatter(repmat(this_x+0.1,size(patient_rt_medians,1),1)+rand(size(patient_rt_medians,1),1)/30-(1/60),patient_rt_medians(:,this_x),16,patient_AFCs/2)
        errorbar(this_x+0.2,mean(patient_rt_medians(patient_AFCs==2,this_x)),std(patient_rt_medians(patient_AFCs==2,this_x))/sqrt(sum(patient_AFCs==2)),'kx')
    errorbar(this_x+0.2,mean(patient_rt_medians(patient_AFCs==4,this_x)),std(patient_rt_medians(patient_AFCs==4,this_x))/sqrt(sum(patient_AFCs==4)),'rx')
    scatter(repmat(this_x-0.1,size(control_rt_medians,1),1)+rand(size(control_rt_medians,1),1)/30-(1/60),control_rt_medians(:,this_x),16,control_AFCs/2)
        errorbar(this_x-0.2,mean(control_rt_medians(control_AFCs==2,this_x)),std(control_rt_medians(control_AFCs==2,this_x))/sqrt(sum(control_AFCs==2)),'kx')
    errorbar(this_x-0.2,mean(control_rt_medians(control_AFCs==4,this_x)),std(control_rt_medians(control_AFCs==4,this_x))/sqrt(sum(control_AFCs==4)),'rx')
end
xlim([0 6])
title('Median RT')
set(gca,'XTick',[0:1:6])
set(gca,'XTickLabel',{'','Match 3','Mismatch 3','','Match 15','Mismatch 15',''},'XTickLabelRotation',15)

%% Now do paired contrasts between conditions
figure
%cmap = colormap(parula(2));
colormap(jet)
set(gcf,'Position',[100 100 1600 800]);
subplot(3,1,1)
hold on
%First perceptual detail
for this_x = 1:3
    scatter(repmat(this_x+0.1,size(patient_response_averages,1),1)+rand(size(patient_response_averages,1),1)/30-(1/60),patient_response_averages(:,this_x+3)-patient_response_averages(:,this_x),16,patient_AFCs/2)
    errorbar(this_x+0.2,mean(patient_response_averages(patient_AFCs==2,this_x+3)-patient_response_averages(patient_AFCs==2,this_x)),std(patient_response_averages(patient_AFCs==2,this_x+3)-patient_response_averages(patient_AFCs==2,this_x))/sqrt(sum(patient_AFCs==2)),'kx')
    errorbar(this_x+0.2,mean(patient_response_averages(patient_AFCs==4,this_x+3)-patient_response_averages(patient_AFCs==4,this_x)),std(patient_response_averages(patient_AFCs==4,this_x+3)-patient_response_averages(patient_AFCs==4,this_x))/sqrt(sum(patient_AFCs==4)),'rx')
    scatter(repmat(this_x-0.1,size(control_response_averages,1),1)+rand(size(control_response_averages,1),1)/30-(1/60),control_response_averages(:,this_x+3)-control_response_averages(:,this_x),16,control_AFCs/2)
    errorbar(this_x-0.2,mean(control_response_averages(control_AFCs==2,this_x+3)-control_response_averages(control_AFCs==2,this_x)),std(control_response_averages(control_AFCs==2,this_x+3)-control_response_averages(control_AFCs==2,this_x))/sqrt(sum(control_AFCs==2)),'kx')
    errorbar(this_x-0.2,mean(control_response_averages(control_AFCs==4,this_x+3)-control_response_averages(control_AFCs==4,this_x)),std(control_response_averages(control_AFCs==4,this_x+3)-control_response_averages(control_AFCs==4,this_x))/sqrt(sum(control_AFCs==4)),'rx')
end
%Then cue congruency
scatter(repmat(4+0.1,size(patient_response_averages,1),1)+rand(size(patient_response_averages,1),1)/30-(1/60),patient_response_averages(:,1)-patient_response_averages(:,2),16,patient_AFCs/2)
errorbar(4+0.2,mean(patient_response_averages(patient_AFCs==2,1)-patient_response_averages(patient_AFCs==2,2)),std(patient_response_averages(patient_AFCs==2,1)-patient_response_averages(patient_AFCs==2,2))/sqrt(sum(patient_AFCs==2)),'kx')
errorbar(4+0.2,mean(patient_response_averages(patient_AFCs==4,1)-patient_response_averages(patient_AFCs==4,2)),std(patient_response_averages(patient_AFCs==4,1)-patient_response_averages(patient_AFCs==4,2))/sqrt(sum(patient_AFCs==4)),'rx')
scatter(repmat(4-0.1,size(control_response_averages,1),1)+rand(size(control_response_averages,1),1)/30-(1/60),control_response_averages(:,1)-control_response_averages(:,2),16,control_AFCs/2)
errorbar(4-0.2,mean(control_response_averages(control_AFCs==2,1)-control_response_averages(control_AFCs==2,2)),std(control_response_averages(control_AFCs==2,1)-control_response_averages(control_AFCs==2,2))/sqrt(sum(control_AFCs==2)),'kx')
errorbar(4-0.2,mean(control_response_averages(control_AFCs==4,1)-control_response_averages(control_AFCs==4,2)),std(control_response_averages(control_AFCs==4,1)-control_response_averages(control_AFCs==4,2))/sqrt(sum(control_AFCs==4)),'rx')

scatter(repmat(5+0.1,size(patient_response_averages,1),1)+rand(size(patient_response_averages,1),1)/30-(1/60),patient_response_averages(:,4)-patient_response_averages(:,5),16,patient_AFCs/2)
errorbar(5+0.2,mean(patient_response_averages(patient_AFCs==2,4)-patient_response_averages(patient_AFCs==2,5)),std(patient_response_averages(patient_AFCs==2,4)-patient_response_averages(patient_AFCs==2,5))/sqrt(sum(patient_AFCs==2)),'kx')
errorbar(5+0.2,mean(patient_response_averages(patient_AFCs==4,4)-patient_response_averages(patient_AFCs==4,5)),std(patient_response_averages(patient_AFCs==4,4)-patient_response_averages(patient_AFCs==4,5))/sqrt(sum(patient_AFCs==4)),'rx')
scatter(repmat(5-0.1,size(control_response_averages,1),1)+rand(size(control_response_averages,1),1)/30-(1/60),control_response_averages(:,4)-control_response_averages(:,5),16,control_AFCs/2)
errorbar(5-0.2,mean(control_response_averages(control_AFCs==2,4)-control_response_averages(control_AFCs==2,5)),std(control_response_averages(control_AFCs==2,4)-control_response_averages(control_AFCs==2,5))/sqrt(sum(control_AFCs==2)),'kx')
errorbar(5-0.2,mean(control_response_averages(control_AFCs==4,4)-control_response_averages(control_AFCs==4,5)),std(control_response_averages(control_AFCs==4,4)-control_response_averages(control_AFCs==4,5))/sqrt(sum(control_AFCs==4)),'rx')

plot([0 6],[0 0],'k--')
xlim([0 6])
title('Percent Correct Difference')
set(gca,'XTick',[0:1:6])
set(gca,'XTickLabel',{'','Match 15-3','Mismatch 15-3','','Match-Mismatch 3','Match-Mismatch 15',''},'XTickLabelRotation',15)

subplot(3,1,2)
hold on
%First perceptual detail
for this_x = 1:3
    scatter(repmat(this_x+0.1,size(patient_rt_averages,1),1)+rand(size(patient_rt_averages,1),1)/30-(1/60),patient_rt_averages(:,this_x+3)-patient_rt_averages(:,this_x),16,patient_AFCs/2)
    errorbar(this_x+0.2,mean(patient_rt_averages(patient_AFCs==2,this_x+3)-patient_rt_averages(patient_AFCs==2,this_x)),std(patient_rt_averages(patient_AFCs==2,this_x+3)-patient_rt_averages(patient_AFCs==2,this_x))/sqrt(sum(patient_AFCs==2)),'kx')
    errorbar(this_x+0.2,mean(patient_rt_averages(patient_AFCs==4,this_x+3)-patient_rt_averages(patient_AFCs==4,this_x)),std(patient_rt_averages(patient_AFCs==4,this_x+3)-patient_rt_averages(patient_AFCs==4,this_x))/sqrt(sum(patient_AFCs==4)),'rx')
    scatter(repmat(this_x-0.1,size(control_rt_averages,1),1)+rand(size(control_rt_averages,1),1)/30-(1/60),control_rt_averages(:,this_x+3)-control_rt_averages(:,this_x),16,control_AFCs/2)
    errorbar(this_x-0.2,mean(control_rt_averages(control_AFCs==2,this_x+3)-control_rt_averages(control_AFCs==2,this_x)),std(control_rt_averages(control_AFCs==2,this_x+3)-control_rt_averages(control_AFCs==2,this_x))/sqrt(sum(control_AFCs==2)),'kx')
    errorbar(this_x-0.2,mean(control_rt_averages(control_AFCs==4,this_x+3)-control_rt_averages(control_AFCs==4,this_x)),std(control_rt_averages(control_AFCs==4,this_x+3)-control_rt_averages(control_AFCs==4,this_x))/sqrt(sum(control_AFCs==4)),'rx')
end
%Then cue congruency
scatter(repmat(4+0.1,size(patient_rt_averages,1),1)+rand(size(patient_rt_averages,1),1)/30-(1/60),patient_rt_averages(:,1)-patient_rt_averages(:,2),16,patient_AFCs/2)
errorbar(4+0.2,mean(patient_rt_averages(patient_AFCs==2,1)-patient_rt_averages(patient_AFCs==2,2)),std(patient_rt_averages(patient_AFCs==2,1)-patient_rt_averages(patient_AFCs==2,2))/sqrt(sum(patient_AFCs==2)),'kx')
errorbar(4+0.2,mean(patient_rt_averages(patient_AFCs==4,1)-patient_rt_averages(patient_AFCs==4,2)),std(patient_rt_averages(patient_AFCs==4,1)-patient_rt_averages(patient_AFCs==4,2))/sqrt(sum(patient_AFCs==4)),'rx')
scatter(repmat(4-0.1,size(control_rt_averages,1),1)+rand(size(control_rt_averages,1),1)/30-(1/60),control_rt_averages(:,1)-control_rt_averages(:,2),16,control_AFCs/2)
errorbar(4-0.2,mean(control_rt_averages(control_AFCs==2,1)-control_rt_averages(control_AFCs==2,2)),std(control_rt_averages(control_AFCs==2,1)-control_rt_averages(control_AFCs==2,2))/sqrt(sum(control_AFCs==2)),'kx')
errorbar(4-0.2,mean(control_rt_averages(control_AFCs==4,1)-control_rt_averages(control_AFCs==4,2)),std(control_rt_averages(control_AFCs==4,1)-control_rt_averages(control_AFCs==4,2))/sqrt(sum(control_AFCs==4)),'rx')

scatter(repmat(5+0.1,size(patient_rt_averages,1),1)+rand(size(patient_rt_averages,1),1)/30-(1/60),patient_rt_averages(:,4)-patient_rt_averages(:,5),16,patient_AFCs/2)
errorbar(5+0.2,mean(patient_rt_averages(patient_AFCs==2,4)-patient_rt_averages(patient_AFCs==2,5)),std(patient_rt_averages(patient_AFCs==2,4)-patient_rt_averages(patient_AFCs==2,5))/sqrt(sum(patient_AFCs==2)),'kx')
errorbar(5+0.2,mean(patient_rt_averages(patient_AFCs==4,4)-patient_rt_averages(patient_AFCs==4,5)),std(patient_rt_averages(patient_AFCs==4,4)-patient_rt_averages(patient_AFCs==4,5))/sqrt(sum(patient_AFCs==4)),'rx')
scatter(repmat(5-0.1,size(control_rt_averages,1),1)+rand(size(control_rt_averages,1),1)/30-(1/60),control_rt_averages(:,4)-control_rt_averages(:,5),16,control_AFCs/2)
errorbar(5-0.2,mean(control_rt_averages(control_AFCs==2,4)-control_rt_averages(control_AFCs==2,5)),std(control_rt_averages(control_AFCs==2,4)-control_rt_averages(control_AFCs==2,5))/sqrt(sum(control_AFCs==2)),'kx')
errorbar(5-0.2,mean(control_rt_averages(control_AFCs==4,4)-control_rt_averages(control_AFCs==4,5)),std(control_rt_averages(control_AFCs==4,4)-control_rt_averages(control_AFCs==4,5))/sqrt(sum(control_AFCs==4)),'rx')

plot([0 6],[0 0],'k--')
xlim([0 6])
title('Mean RT Difference')
set(gca,'XTick',[0:1:6])
set(gca,'XTickLabel',{'','Match 15-3','Mismatch 15-3','','Match-Mismatch 3','Match-Mismatch 15',''},'XTickLabelRotation',15)

subplot(3,1,3)
hold on
%First perceptual detail
for this_x = 1:3
    scatter(repmat(this_x+0.1,size(patient_rt_medians,1),1)+rand(size(patient_rt_medians,1),1)/30-(1/60),patient_rt_medians(:,this_x+3)-patient_rt_medians(:,this_x),16,patient_AFCs/2)
    errorbar(this_x+0.2,mean(patient_rt_medians(patient_AFCs==2,this_x+3)-patient_rt_medians(patient_AFCs==2,this_x)),std(patient_rt_medians(patient_AFCs==2,this_x+3)-patient_rt_medians(patient_AFCs==2,this_x))/sqrt(sum(patient_AFCs==2)),'kx')
    errorbar(this_x+0.2,mean(patient_rt_medians(patient_AFCs==4,this_x+3)-patient_rt_medians(patient_AFCs==4,this_x)),std(patient_rt_medians(patient_AFCs==4,this_x+3)-patient_rt_medians(patient_AFCs==4,this_x))/sqrt(sum(patient_AFCs==4)),'rx')
    scatter(repmat(this_x-0.1,size(control_rt_medians,1),1)+rand(size(control_rt_medians,1),1)/30-(1/60),control_rt_medians(:,this_x+3)-control_rt_medians(:,this_x),16,control_AFCs/2)
    errorbar(this_x-0.2,mean(control_rt_medians(control_AFCs==2,this_x+3)-control_rt_medians(control_AFCs==2,this_x)),std(control_rt_medians(control_AFCs==2,this_x+3)-control_rt_medians(control_AFCs==2,this_x))/sqrt(sum(control_AFCs==2)),'kx')
    errorbar(this_x-0.2,mean(control_rt_medians(control_AFCs==4,this_x+3)-control_rt_medians(control_AFCs==4,this_x)),std(control_rt_medians(control_AFCs==4,this_x+3)-control_rt_medians(control_AFCs==4,this_x))/sqrt(sum(control_AFCs==4)),'rx')
end
%Then cue congruency
scatter(repmat(4+0.1,size(patient_rt_medians,1),1)+rand(size(patient_rt_medians,1),1)/30-(1/60),patient_rt_medians(:,1)-patient_rt_medians(:,2),16,patient_AFCs/2)
errorbar(4+0.2,mean(patient_rt_medians(patient_AFCs==2,1)-patient_rt_medians(patient_AFCs==2,2)),std(patient_rt_medians(patient_AFCs==2,1)-patient_rt_medians(patient_AFCs==2,2))/sqrt(sum(patient_AFCs==2)),'kx')
errorbar(4+0.2,mean(patient_rt_medians(patient_AFCs==4,1)-patient_rt_medians(patient_AFCs==4,2)),std(patient_rt_medians(patient_AFCs==4,1)-patient_rt_medians(patient_AFCs==4,2))/sqrt(sum(patient_AFCs==4)),'rx')
scatter(repmat(4-0.1,size(control_rt_medians,1),1)+rand(size(control_rt_medians,1),1)/30-(1/60),control_rt_medians(:,1)-control_rt_medians(:,2),16,control_AFCs/2)
errorbar(4-0.2,mean(control_rt_medians(control_AFCs==2,1)-control_rt_medians(control_AFCs==2,2)),std(control_rt_medians(control_AFCs==2,1)-control_rt_medians(control_AFCs==2,2))/sqrt(sum(control_AFCs==2)),'kx')
errorbar(4-0.2,mean(control_rt_medians(control_AFCs==4,1)-control_rt_medians(control_AFCs==4,2)),std(control_rt_medians(control_AFCs==4,1)-control_rt_medians(control_AFCs==4,2))/sqrt(sum(control_AFCs==4)),'rx')

scatter(repmat(5+0.1,size(patient_rt_medians,1),1)+rand(size(patient_rt_medians,1),1)/30-(1/60),patient_rt_medians(:,4)-patient_rt_medians(:,5),16,patient_AFCs/2)
errorbar(5+0.2,mean(patient_rt_medians(patient_AFCs==2,4)-patient_rt_medians(patient_AFCs==2,5)),std(patient_rt_medians(patient_AFCs==2,4)-patient_rt_medians(patient_AFCs==2,5))/sqrt(sum(patient_AFCs==2)),'kx')
errorbar(5+0.2,mean(patient_rt_medians(patient_AFCs==4,4)-patient_rt_medians(patient_AFCs==4,5)),std(patient_rt_medians(patient_AFCs==4,4)-patient_rt_medians(patient_AFCs==4,5))/sqrt(sum(patient_AFCs==4)),'rx')
scatter(repmat(5-0.1,size(control_rt_medians,1),1)+rand(size(control_rt_medians,1),1)/30-(1/60),control_rt_medians(:,4)-control_rt_medians(:,5),16,control_AFCs/2)
errorbar(5-0.2,mean(control_rt_medians(control_AFCs==2,4)-control_rt_medians(control_AFCs==2,5)),std(control_rt_medians(control_AFCs==2,4)-control_rt_medians(control_AFCs==2,5))/sqrt(sum(control_AFCs==2)),'kx')
errorbar(5-0.2,mean(control_rt_medians(control_AFCs==4,4)-control_rt_medians(control_AFCs==4,5)),std(control_rt_medians(control_AFCs==4,4)-control_rt_medians(control_AFCs==4,5))/sqrt(sum(control_AFCs==4)),'rx')

plot([0 6],[0 0],'k--')
xlim([0 6])
title('Median RT Difference')
set(gca,'XTick',[0:1:6])
set(gca,'XTickLabel',{'','Match 15-3','Mismatch 15-3','','Match-Mismatch 3','Match-Mismatch 15',''},'XTickLabelRotation',15)

saveas(gcf,'Inscanner_Behavioural_Effects.png')
saveas(gcf,'Inscanner_Behavioural_Effects.pdf')



