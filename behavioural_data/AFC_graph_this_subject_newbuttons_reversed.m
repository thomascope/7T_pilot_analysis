function [ all_response_averages, all_rt_averages] = AFC_graph_this_subject_newbuttons_reversed( subject, date )

all_files = ls(['*' subject '*' num2str(date) '.mat']);
%all_files = ['AFC_7T_' subject '_Run_' num2str(runnum) '_' num2str(date) '.mat'];

all_runs_resps = []; 
all_runs_rts = [];
all_runs_resp_corr = [];
for i = size(all_files,1):-1:1
    load(all_files(i,:))
    all_runs_resps = [all_runs_resps, resp(resp~=0)];
    resp_corr = abs(resp(resp~=0)-5)==all_trial_targets';
    all_runs_resp_corr = [all_runs_resp_corr, resp_corr];
    all_runs_rts = [all_runs_rts, all_rts(all_rts~=0)];
end

nanpadded_runs_corr = [all_runs_resp_corr NaN(1,(size(response_order.this_vocoder_channels,2)-size(all_runs_resp_corr,2)))];
nanpadded_runs_resps =  [all_runs_resps NaN(1,(size(response_order.this_vocoder_channels,2)-size(all_runs_resps,2)))];
nanpadded_runs_rts = [all_runs_rts NaN(1,(size(response_order.this_vocoder_channels,2)-size(all_runs_rts,2)))];
all_response_averages = 100*[nanmean(nanpadded_runs_corr(response_order.this_vocoder_channels==1&response_order.this_cue_types==1)) nanmean(nanpadded_runs_corr(response_order.this_vocoder_channels==1&response_order.this_cue_types==2)) nanmean(nanpadded_runs_corr(response_order.this_vocoder_channels==1&response_order.this_cue_types==3)) nanmean(nanpadded_runs_corr(response_order.this_vocoder_channels==2&response_order.this_cue_types==1)) nanmean(nanpadded_runs_corr(response_order.this_vocoder_channels==2&response_order.this_cue_types==2)) nanmean(nanpadded_runs_corr(response_order.this_vocoder_channels==2&response_order.this_cue_types==3))]
figure
bar(all_response_averages)
ylim([0 100])
title('Percent Correct')
set(gca,'XTickLabel',{'Match 3','Mismatch 3','Neutral 3','Match 16','Mismatch 16','Neutral 16'},'XTickLabelRotation',15)
all_rt_averages = [nanmean(nanpadded_runs_rts(response_order.this_vocoder_channels==1&response_order.this_cue_types==1)) nanmean(nanpadded_runs_rts(response_order.this_vocoder_channels==1&response_order.this_cue_types==2)) nanmean(nanpadded_runs_rts(response_order.this_vocoder_channels==1&response_order.this_cue_types==3)) nanmean(nanpadded_runs_rts(response_order.this_vocoder_channels==2&response_order.this_cue_types==1)) nanmean(nanpadded_runs_rts(response_order.this_vocoder_channels==2&response_order.this_cue_types==2)) nanmean(nanpadded_runs_rts(response_order.this_vocoder_channels==2&response_order.this_cue_types==3))]
figure
bar(all_rt_averages)
title('Mean RT')
set(gca,'XTickLabel',{'Match 3','Mismatch 3','Neutral 3','Match 16','Mismatch 16','Neutral 16'},'XTickLabelRotation',15)
all_rt_averages = [nanmedian(nanpadded_runs_rts(response_order.this_vocoder_channels==1&response_order.this_cue_types==1)) nanmedian(nanpadded_runs_rts(response_order.this_vocoder_channels==1&response_order.this_cue_types==2)) nanmedian(nanpadded_runs_rts(response_order.this_vocoder_channels==1&response_order.this_cue_types==3)) nanmedian(nanpadded_runs_rts(response_order.this_vocoder_channels==2&response_order.this_cue_types==1)) nanmedian(nanpadded_runs_rts(response_order.this_vocoder_channels==2&response_order.this_cue_types==2)) nanmedian(nanpadded_runs_rts(response_order.this_vocoder_channels==2&response_order.this_cue_types==3))]
figure
bar(all_rt_averages)
title('Median RT')
set(gca,'XTickLabel',{'Match 3','Mismatch 3','Neutral 3','Match 16','Mismatch 16','Neutral 16'},'XTickLabelRotation',15)

end

