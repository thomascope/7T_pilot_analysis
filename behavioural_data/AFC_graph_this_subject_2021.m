function [all_response_averages, all_rt_averages, all_rt_medians, AFCs, running_average, running_average_bycond, normalised_running_average, response_order] = AFC_graph_this_subject_2021(subject, date, graph_this)

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
all_trialtypes = {}; % For debugging/checking
response_trialtypes = {};
%for i = size(all_files,1):-1:1 % Was incrementing backwards here - I think this was a double negative mistake. Added check below to compare recorded display output against the recorded cue and channel number
for i = 1:size(all_files,1)
    resp = [];
    resp_corr = [];
    all_rts = [];
    load(all_files(i).name)
    num_expected_resps = size(response_order.all_trial_codes,2);
    
    if AFCs == 2
        resp(resp~=0) = 2-mod(resp(resp~=0),2);
        resp_corr = resp(resp~=0)==all_trial_targets';
        resp_corr = double(resp_corr);
        resp_corr(isnan(resp(resp~=0))) = NaN;
        if nanmean(resp_corr) < 0.45 %Likely buttons reversed by subject
            resp(resp~=0) = mod(resp(resp~=0),2)+1;
            resp_corr = resp(resp~=0)==all_trial_targets';
            resp_corr = double(resp_corr);
            resp_corr(isnan(resp(resp~=0))) = NaN;
        end
    elseif AFCs == 4
        resp_corr = resp(resp~=0)==all_trial_targets';
        resp_corr = double(resp_corr);
        resp_corr(isnan(resp(resp~=0))) = NaN;
        if nanmean(resp_corr) < 0.3 %Likely buttons reversed by subject
            resp_corr = 5-resp(resp~=0)==all_trial_targets';
            resp_corr = double(resp_corr);
            resp_corr(isnan(resp(resp~=0))) = NaN;
        end
    end
    all_trialtypes = [all_trialtypes, trialtype];
    these_response_trialtypes = trialtype(strncmp('Response',trialtype,8));
    resp_corr = double(resp_corr);
    if length(these_response_trialtypes) < num_expected_resps %In case of trunkated runs
        these_response_trialtypes(length(these_response_trialtypes)+1:num_expected_resps) = {NaN};
        resp_corr(length(resp_corr)+1:num_expected_resps) = NaN;
        all_rts(length(all_rts)+1:num_expected_resps) = NaN;
        resp(length(resp)+1:num_expected_resps) = NaN;
    end
    all_runs_resps = [all_runs_resps, resp(resp~=0)];
    all_runs_resp_corr = [all_runs_resp_corr, resp_corr];
    all_runs_rts = [all_runs_rts, all_rts(all_rts~=0)];
    response_trialtypes = [response_trialtypes, these_response_trialtypes];
end

%Check the trial encoding is as expected
for i = 1:length(response_trialtypes)
    if ~ isnan(response_trialtypes{i})
    split_colon = strsplit(response_trialtypes{i},': ');
    split_underscores = strsplit(split_colon{2},'_');
    derived_response_cue{i} = split_underscores{1};
    derived_vocoder_channels(i) = str2num(split_underscores{2});
    if strncmp(split_underscores{1},'Match',5)
        derived_response_cue_code(i) = 1;
    elseif strncmp(split_underscores{1},'MisMatch',8)
        derived_response_cue_code(i) = 2;
    else
        error(['Unknown cue type ' split_underscores{1}])
    end
    if derived_vocoder_channels(i) == 3
        derived_vocoder_channel_code(i) = 1;
    elseif derived_vocoder_channels(i) == 15
        derived_vocoder_channel_code(i) = 2;
    else
        error(['Unknown vocoder number ' split_underscores{2}])
    end
    else
        derived_vocoder_channel_code(i) = NaN;
        derived_response_cue_code(i) = NaN;
    end
end

assert(all(derived_vocoder_channel_code == response_order.this_vocoder_channels(1:length(derived_vocoder_channel_code)) | isnan(derived_vocoder_channel_code)),'Vocoder channels do not match in recorded and derived formats')
assert(all(derived_response_cue_code == response_order.this_cue_types(1:length(derived_response_cue_code)) | isnan(derived_response_cue_code)),'Cue types do not match in recorded and derived formats')

if AFCs == 2
    poor_thresh = 0.6;
elseif AFCs == 4
    poor_thresh = 0.4;
end
if nansum(all_runs_resp_corr)/length(all_runs_resp_corr) < poor_thresh
    warning(['Poor performance for subject ' subject ' not fixed by reversing buttons - you will want to check their data'])
end

running_average = movmean(all_runs_resp_corr,5,'omitnan');
running_average_bycond(:,1) = movmean(all_runs_resp_corr(response_order.this_vocoder_channels==1&response_order.this_cue_types==1),5,'omitnan');
running_average_bycond(:,2) = movmean(all_runs_resp_corr(response_order.this_vocoder_channels==1&response_order.this_cue_types==2),5,'omitnan');
running_average_bycond(:,3) = movmean(all_runs_resp_corr(response_order.this_vocoder_channels==2&response_order.this_cue_types==1),5,'omitnan');
running_average_bycond(:,4) = movmean(all_runs_resp_corr(response_order.this_vocoder_channels==2&response_order.this_cue_types==2),5,'omitnan');

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

%Now normalise running averages for single subject difficulty of trial types
normalised_all_runs_resp_corr = nan(size(all_runs_resp_corr));
normalised_all_runs_resp_corr(response_order.this_vocoder_channels==1&response_order.this_cue_types==1) = all_runs_resp_corr(response_order.this_vocoder_channels==1&response_order.this_cue_types==1)/all_response_averages(1)*100;
normalised_all_runs_resp_corr(response_order.this_vocoder_channels==1&response_order.this_cue_types==2) = all_runs_resp_corr(response_order.this_vocoder_channels==1&response_order.this_cue_types==2)/all_response_averages(2)*100;
normalised_all_runs_resp_corr(response_order.this_vocoder_channels==2&response_order.this_cue_types==1) = all_runs_resp_corr(response_order.this_vocoder_channels==2&response_order.this_cue_types==1)/all_response_averages(4)*100;
normalised_all_runs_resp_corr(response_order.this_vocoder_channels==2&response_order.this_cue_types==2) = all_runs_resp_corr(response_order.this_vocoder_channels==2&response_order.this_cue_types==2)/all_response_averages(5)*100;

normalised_running_average = movmean(normalised_all_runs_resp_corr,5,'omitnan');
