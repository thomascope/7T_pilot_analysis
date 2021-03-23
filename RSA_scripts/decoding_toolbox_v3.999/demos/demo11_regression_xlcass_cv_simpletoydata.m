% This script is a demo showing regression on toy data.
% The toy data are simple matlab matrices and no "real" fMRI or EEG data.
% At the end, the example also demonstrates how to train on partially
% different data than to test on (still using leave-one-run-out CV)
%
% Kai, 29.3.2018

clear all
dbstop if error % if something goes wrong

fpath = fileparts(fileparts(mfilename('fullpath')));
addpath(fpath)

% initialize TDT & cfg
cfg = decoding_defaults;

%% Set the output directory where data will be saved
% cfg.results.dir = % e.g. 'toyresults'
cfg.results.write = 0; % no results are written to disk

%% generate some toy data
% define number of "runs" and center means
nruns = 4; % lets simulate we have n runs

set1.mean = [5 2]; % should have the same dim as set1, otherwise it wont work (and would not make sense, either)
set2.mean = [0 0];
set3.mean = [0 0];
set4.mean = [-2 -5];

% add guassian noise
data1 = randn(nruns, 2) + repmat(set1.mean, nruns, 1);
data2 = randn(nruns, 2) + repmat(set2.mean, nruns, 1);
data3 = randn(nruns, 2) + repmat(set3.mean, nruns, 1);
data4 = randn(nruns, 2) + repmat(set4.mean, nruns, 1);

% alternative: uniform
% data1 = rand(nruns, length(set1.mean)) + repmat(set1.mean, nruns, 1);
% data2 = rand(nruns, length(set2.mean)) + repmat(set2.mean, nruns, 1);

% put all together in a data matrix
data = [data1; data2; data3; data4]; 

%% add data description
% save labels
cfg.files.label = [-1*ones(size(data1,1), 1); 
    0*ones(size(data2,1), 1);
    0*ones(size(data3,1), 1);
    1*ones(size(data4,1), 1)];

class_names = [repmat('A', nruns, 1); repmat('B', nruns, 1); repmat('C', nruns, 1); repmat('D', nruns, 1)];

% save run number
cfg.files.chunk = [1:nruns, 1:nruns, 1:nruns, 1:nruns]';

all_chunks = unique(cfg.files.chunk);
all_labels = unique(cfg.files.label);

% save a description
for ifile = 1:length(cfg.files.label)
    curr_label = cfg.files.label(ifile);
    curr_chunk = cfg.files.chunk(ifile);
    cfg.files.name(ifile) = {sprintf('class%s_label%i_run%i', class_names(ifile), curr_label, curr_chunk)};
end

% add an empty mask
cfg.files.mask = '';

%% plot the data (if 2d)
if size(data, 2) == 2
    resfig = figure('name', 'Data');
    scatter(data(:, 1), data(:, 2), 30, cfg.files.label);
end

%% Prepare data for passing
passed_data.data = data;
passed_data.mask_index = 1:size(data, 2); % use all voxels
passed_data.files = cfg.files;
passed_data.hdr = ''; % we don't need a header, because we don't write img-files as output (but mat-files)
passed_data.dim = [length(set1.mean), 1, 1]; % add dimension information of the original data
% passed_data.voxelsize = [1 1 1];


%% Add defaults for the remaining parameters that we did not specify
cfg = decoding_defaults(cfg);

% Set the analysis that should be performed (here we only want to do 1
% decoding)
cfg.analysis = 'wholebrain';
cfg.decoding.method = 'regression';
cfg.acknowledge_transres_SVM_pattern_for_regression_is_experimental = 1; % necessary to get pattern for REGRESSION (ok for classification, but for regression the implementation is still experimental, see occuring warnings)
cfg.results.output = {'corr', 'SVM_weights', 'SVM_pattern'}; % 'model_parameters' returns the full models (warning, might be huge, one model for each step)
                               % To get the filterweights or the pattern, use or 'SVM_weights', 'SVM_weights_plus_bias'. 
                               % 'SVM_pattern' or 'SVM_pattern_alldata'
                               % 
                               % Important:
                               % See Haufe et al, 2014, Neuroimage for the the important distinction 
                               % between Filters (aka "Weightmaps") and Patterns
                                                       
% cfg.results.output = {'accuracy', 'model_parameters', 'SVM_weights', 'SVM_weights_plusbias'}; % example for more possible outputs

%% Nothing needs to be changed below for a standard leave-one-run out cross validation analysis.
% Create a leave-one-run-out cross validation design:

cfg.design = make_design_cv(cfg); 
cfg.fighandles.plot_design = plot_design(cfg);
cfg.plot_selected_voxels = 1;

%% Acknowledge that label 0 occurs more often than -1 and 1
cfg.design.unbalanced_data = 'ok';

%% Change the design: Remove one class 0 in test and the other in train
% get index of class B
class_B_filt = strncmp(cfg.files.name, 'classB', length('classB'));
class_C_filt = strncmp(cfg.files.name, 'classC', length('classC'));

% create design train on B test on C
cfg_BC = cfg;
cfg_BC.results.filestart = 'resBC'; % CHANGE NAME of result files
cfg_BC.design.train(class_C_filt, :) = 0; % remove C from train set
cfg_BC.design.test(class_B_filt, :) = 0; % remove B from test set
plot_design(cfg_BC);

% create design train on C test on B
cfg_CB = cfg;
cfg_CB.results.filestart = 'resCB'; % CHANGE NAME of result files
cfg_CB.design.train(class_B_filt, :) = 0; 
cfg_CB.design.test(class_C_filt, :) = 0;
plot_design(cfg_CB);

%% Run decoding (all)

[results, cfg] = decoding(cfg, passed_data);

display('Result full CV (correlation):')
results.corr

%% Run decoding (train B test C)

[results_BC, cfg_BC] = decoding(cfg_BC, passed_data);

display('Result train B test C CV (correlation):')
results_BC.corr

%% Run decoding (train C test B)

[results_CB, cfg_CB] = decoding(cfg_CB, passed_data);

display('Result train C test B CV (correlation):')
results_CB.corr

%% Combine both if you like
cfg_BCCB = combine_designs(cfg_BC, cfg_CB); 
cfg_BCCB.results.filestart = 'resBCCB';
plot_design(cfg_BCCB);
[results_BCCB, cfg_BCCB] = decoding(cfg_BCCB, passed_data);
display('Result train B test C, then C B, both CV (averaged correlation [better: zcorr]):')
results_BCCB.corr
