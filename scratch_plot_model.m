model_run_date = '13-May-2021';
try
    load(['./modelparameters/modelparameters_' model_run_date '.mat'])
catch
    [all_sigma_pred,all_thresholds,controls_sigma_pred,controls_threshold,patients_sigma_pred,patients_threshold] = module_bayesian_behaviour(subjects,group,dates);
    save(['./modelparameters/modelparameters_' date '.mat'],'all_sigma_pred','all_thresholds','controls_sigma_pred','controls_threshold','patients_sigma_pred','patients_threshold');
end

figure
scatter(ones(1,length(controls_sigma_pred)),controls_sigma_pred(1,:))
hold on
scatter(2*ones(1,length(patients_sigma_pred)),patients_sigma_pred(1,:))
scatter(4*ones(1,length(controls_sigma_pred)),controls_sigma_pred(2,:))
scatter(5*ones(1,length(patients_sigma_pred)),patients_sigma_pred(2,:))
scatter(7*ones(1,length(controls_sigma_pred)),nanmean(controls_sigma_pred))
scatter(8*ones(1,length(patients_sigma_pred)),nanmean(patients_sigma_pred))
xlim([0 9])
set(gca,'XTick',[1:8],'XtickLabel',({'C Pre','P Pre','','C Post','P Post','','C Mean','P Mean'}),'XTickLabelRotation',45)
title('Standard Deviation of Prior')

figure
scatter(ones(1,length(controls_threshold)),controls_threshold(1,:))
hold on
scatter(2*ones(1,length(patients_threshold)),patients_threshold(1,:))
scatter(4*ones(1,length(controls_threshold)),controls_threshold(2,:))
scatter(5*ones(1,length(patients_threshold)),patients_threshold(2,:))
scatter(7*ones(1,length(controls_threshold)),nanmean(controls_threshold))
scatter(8*ones(1,length(patients_threshold)),nanmean(patients_threshold))
xlim([0 9])
set(gca,'XTick',[1:8],'XtickLabel',({'C Pre','P Pre','','C Post','P Post','','C Mean','P Mean'}),'XTickLabelRotation',45)
title('Perceptual Threshold')