function [all_response_averages, all_rt_averages, all_rt_medians, AFCs] = AFC_graph_this_subject_2021(subject, date, graph_this)

if ~exist('graph_this','var')
    graph_this = 0;
end

all_files = dir(['*' subject '*' num2str(date) '.mat']);
if isempty(all_files)
    all_files = dir(['*' subject '*' date '.mat']);
end
if isempty(all_files)
    error(['No files found for ' subject])
end
%all_files = ['AFC_7T_' subject '_Run_' num2str(runnum) '_' num2str(date) '.mat'];

all_runs_resps = []; 
all_runs_rts = [];
all_runs_resp_corr = [];
for i = size(all_files,1):-1:1
    load(all_files(i).name)
    all_runs_resps = [all_runs_resps, resp(resp~=0)];
    if AFCs == 2
        resp(resp~=0) = 2-mod(resp(resp~=0),2);
        resp_corr = resp(resp~=0)==all_trial_targets';
        if mean(resp_corr) < 0.45 %Likely buttons reversed by subject
            resp(resp~=0) = mod(resp(resp~=0),2)+1;
            resp_corr = resp(resp~=0)==all_trial_targets';
        end
    elseif AFCs == 4
        resp_corr = resp(resp~=0)==all_trial_targets';
        if mean(resp_corr) < 0.3 %Likely buttons reversed by subject
            resp_corr = 5-resp(resp~=0)==all_trial_targets';
        end
    end
    all_runs_resp_corr = [all_runs_resp_corr, resp_corr];
    all_runs_rts = [all_runs_rts, all_rts(all_rts~=0)];
end

if AFCs == 2
    poor_thresh = 0.6;
elseif AFCs == 4
    poor_thresh = 0.4;
end
if sum(all_runs_resp_corr)/length(all_runs_resp_corr) < poor_thresh
    warning(['Poor performance for subject ' subject ' not fixed by reversing buttons - you will want to check their data'])
end

nanpadded_runs_corr = [all_runs_resp_corr NaN(1,(size(response_order.this_vocoder_channels,2)-size(all_runs_resp_corr,2)))];
nanpadded_runs_resps =  [all_runs_resps NaN(1,(size(response_order.this_vocoder_channels,2)-size(all_runs_resps,2)))];
nanpadded_runs_rts = [all_runs_rts NaN(1,(size(response_order.this_vocoder_channels,2)-size(all_runs_rts,2)))];
all_response_averages = 100*[nanmean(nanpadded_runs_corr(response_order.this_vocoder_channels==1&response_order.this_cue_types==1)) nanmean(nanpadded_runs_corr(response_order.this_vocoder_channels==1&response_order.this_cue_types==2)) nanmean(nanpadded_runs_corr(response_order.this_vocoder_channels==1&response_order.this_cue_types==3)) nanmean(nanpadded_runs_corr(response_order.this_vocoder_channels==2&response_order.this_cue_types==1)) nanmean(nanpadded_runs_corr(response_order.this_vocoder_channels==2&response_order.this_cue_types==2)) nanmean(nanpadded_runs_corr(response_order.this_vocoder_channels==2&response_order.this_cue_types==3))]
if graph_this
    figure
    bar(all_response_averages)
    ylim([0 100])
    title('Percent Correct')
    set(gca,'XTickLabel',{'Match 3','Mismatch 3','','Match 15','Mismatch 15',''},'XTickLabelRotation',15)
    all_rt_averages = [nanmean(nanpadded_runs_rts(response_order.this_vocoder_channels==1&response_order.this_cue_types==1)) nanmean(nanpadded_runs_rts(response_order.this_vocoder_channels==1&response_order.this_cue_types==2)) nanmean(nanpadded_runs_rts(response_order.this_vocoder_channels==1&response_order.this_cue_types==3)) nanmean(nanpadded_runs_rts(response_order.this_vocoder_channels==2&response_order.this_cue_types==1)) nanmean(nanpadded_runs_rts(response_order.this_vocoder_channels==2&response_order.this_cue_types==2)) nanmean(nanpadded_runs_rts(response_order.this_vocoder_channels==2&response_order.this_cue_types==3))];
    figure
    bar(all_rt_averages)
    title('Mean RT')
    set(gca,'XTickLabel',{'Match 3','Mismatch 3','','Match 15','Mismatch 15',''},'XTickLabelRotation',15)
    all_rt_medians = [nanmedian(nanpadded_runs_rts(response_order.this_vocoder_channels==1&response_order.this_cue_types==1)) nanmedian(nanpadded_runs_rts(response_order.this_vocoder_channels==1&response_order.this_cue_types==2)) nanmedian(nanpadded_runs_rts(response_order.this_vocoder_channels==1&response_order.this_cue_types==3)) nanmedian(nanpadded_runs_rts(response_order.this_vocoder_channels==2&response_order.this_cue_types==1)) nanmedian(nanpadded_runs_rts(response_order.this_vocoder_channels==2&response_order.this_cue_types==2)) nanmedian(nanpadded_runs_rts(response_order.this_vocoder_channels==2&response_order.this_cue_types==3))];
    figure
    bar(all_rt_medians)
    title('Median RT')
    set(gca,'XTickLabel',{'Match 3','Mismatch 3','','Match 15','Mismatch 15',''},'XTickLabelRotation',15)
else
    all_rt_averages = [nanmean(nanpadded_runs_rts(response_order.this_vocoder_channels==1&response_order.this_cue_types==1)) nanmean(nanpadded_runs_rts(response_order.this_vocoder_channels==1&response_order.this_cue_types==2)) nanmean(nanpadded_runs_rts(response_order.this_vocoder_channels==1&response_order.this_cue_types==3)) nanmean(nanpadded_runs_rts(response_order.this_vocoder_channels==2&response_order.this_cue_types==1)) nanmean(nanpadded_runs_rts(response_order.this_vocoder_channels==2&response_order.this_cue_types==2)) nanmean(nanpadded_runs_rts(response_order.this_vocoder_channels==2&response_order.this_cue_types==3))];
    all_rt_medians = [nanmedian(nanpadded_runs_rts(response_order.this_vocoder_channels==1&response_order.this_cue_types==1)) nanmedian(nanpadded_runs_rts(response_order.this_vocoder_channels==1&response_order.this_cue_types==2)) nanmedian(nanpadded_runs_rts(response_order.this_vocoder_channels==1&response_order.this_cue_types==3)) nanmedian(nanpadded_runs_rts(response_order.this_vocoder_channels==2&response_order.this_cue_types==1)) nanmedian(nanpadded_runs_rts(response_order.this_vocoder_channels==2&response_order.this_cue_types==2)) nanmedian(nanpadded_runs_rts(response_order.this_vocoder_channels==2&response_order.this_cue_types==3))];
end

