all_subj_names = strcat(subjects{:},' ');
all_subj_names = strrep(all_subj_names,'P7',' P7');
cmd = ['aparcstats2table --subjects' all_subj_names ' --hemi lh --meas thickness --tablefile ./freesurfer_stats/all_subj_thicknesses'];
mkdir('./freesurfer_stats/')
system(cmd);
all_thicknessess = readtable('./freesurfer_stats/all_subj_thicknesses','Delimiter','tab','ReadRowNames',true);

all_patient_thickness_table = table();
all_control_thickness_table = table();
for crun = 1:length(subjects)
    if group(crun) == 1
        all_control_thickness_table(end+1,:) = all_thicknessess(subjects{crun},:);
    else
        all_patient_thickness_table(end+1,:) = all_thicknessess(subjects{crun},:);
    end
end
figure
hold on
title('Thicknesses')
h = [];
p = [];
for i = 1:size(all_patient_thickness_table,2)
    [h(i),p(i)] = ttest2(all_patient_thickness_table{:,i},all_control_thickness_table{:,i});
    scatter(i*ones(1,size(all_patient_thickness_table,1))-0.1,all_patient_thickness_table{:,i},'b')
    scatter(i*ones(1,size(all_control_thickness_table,1))+0.1,all_control_thickness_table{:,i},'r')
    if h(i)
        plot(i,3.4,'k*')
        errorbar(i-0.1,mean(all_patient_thickness_table{:,i}),std(all_patient_thickness_table{:,i})/sqrt(size(all_patient_thickness_table,1)),'-s','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue','LineWidth',1,'color','black')
        errorbar(i+0.1,mean(all_control_thickness_table{:,i}),std(all_control_thickness_table{:,i})/sqrt(size(all_control_thickness_table,1)),'-s','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',1,'color','black')
    else
        errorbar(i-0.1,mean(all_patient_thickness_table{:,i}),std(all_patient_thickness_table{:,i})/sqrt(size(all_patient_thickness_table,1)),'-s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',1,'color','black')
        errorbar(i+0.1,mean(all_control_thickness_table{:,i}),std(all_control_thickness_table{:,i})/sqrt(size(all_control_thickness_table,1)),'-s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',1,'color','black')
        
    end
    
end
set(gca,'XTick',1:i)
set(gca,'XTickLabel',all_patient_thickness_table.Properties.VariableNames,'XTickLabelRotation',90,'TickLabelInterpreter','None')

significant_regions_thickness = all_patient_thickness_table.Properties.VariableNames(logical(h));

% cmd = ['asegstats2table --subjects' all_subj_names ' --meas mean --tablefile ./freesurfer_stats/all_subj_mean_intensity'];
% mkdir('./freesurfer_stats/')
% system(cmd);
% all_mean_intensitys = readtable('./freesurfer_stats/all_subj_mean_intensity','Delimiter','tab','ReadRowNames',true);
%
% all_patient_intensity_table = table();
% all_control_intensity_table = table();
% for crun = 1:length(subjects)
%     if group(crun) == 1
%         all_control_intensity_table(end+1,:) = all_mean_intensitys(subjects{crun},:);
%     else
%         all_patient_intensity_table(end+1,:) = all_mean_intensitys(subjects{crun},:);
%     end
% end
% figure
% hold on
% title('mean_intensity','Interpreter','none')
% h = [];
% p = [];
% for i = 1:size(all_patient_intensity_table,2)
%     [h(i),p(i)] = ttest2(all_patient_intensity_table{:,i},all_control_intensity_table{:,i});
%     scatter(i*ones(1,size(all_patient_intensity_table,1))-0.1,all_patient_intensity_table{:,i},'b')
%     scatter(i*ones(1,size(all_control_intensity_table,1))+0.1,all_control_intensity_table{:,i},'r')
%        if ~isnan(h(i)) && h(i)
%            plot(i,3.4,'k*')
%        end
% end
% set(gca,'XTick',1:i)
% set(gca,'XTickLabel',all_patient_intensity_table.Properties.VariableNames,'XTickLabelRotation',90,'TickLabelInterpreter','None')
%
% significant_regions_intensity = all_patient_intensity_table.Properties.VariableNames(logical(h));
%
% cmd = ['aparcstats2table --subjects' all_subj_names ' --hemi lh --tablefile ./freesurfer_stats/all_subj_areaes'];
% mkdir('./freesurfer_stats/')
% system(cmd);
% all_areaess = readtable('./freesurfer_stats/all_subj_areaes','Delimiter','tab','ReadRowNames',true);
%
% all_patient_area_table = table();
% all_control_area_table = table();
% for crun = 1:length(subjects)
%     if group(crun) == 1
%         all_control_area_table(end+1,:) = all_areaess(subjects{crun},:);
%     else
%         all_patient_area_table(end+1,:) = all_areaess(subjects{crun},:);
%     end
% end
% figure
% hold on
% title('areas')
% h = [];
% p = [];
% for i = 1:size(all_patient_area_table,2)-1 % Exclude WM last
%     [h(i),p(i)] = ttest2(all_patient_area_table{:,i},all_control_area_table{:,i});
%     scatter(i*ones(1,size(all_patient_area_table,1))-0.1,all_patient_area_table{:,i},'b')
%     scatter(i*ones(1,size(all_control_area_table,1))+0.1,all_control_area_table{:,i},'r')
%        if ~isnan(h(i)) && h(i)
%            plot(i,11000,'k*')
%        end
% end
% set(gca,'XTick',1:i)
% set(gca,'XTickLabel',all_patient_area_table.Properties.VariableNames,'XTickLabelRotation',90,'TickLabelInterpreter','None')
%
% significant_regions_area = all_patient_area_table.Properties.VariableNames(logical(h));
