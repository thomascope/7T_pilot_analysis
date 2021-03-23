% This script is a template who want to get an interpretable weight map of their wholebrain
% classifier (i.e. backtransformed pattern aka the "Haufe-method" (see
% Haufe et al., 2014) and is particularly well suited for people who have
% betas available from an SPM.mat. 
% If you don't have an SPM.mat available, then use
% decoding_template_nobetas.m and adjust the script to reflect wholebrain
% analysis and make_design_alldata

% The analysis involves two steps:
% (1) we need to make sure using cross-validation that our data generalizes
% well. If it doesn't it is difficult to interpret the weights because they
% are unstable.
% (2) then using all data we will generate a weight map. We use all data
% because we showed that using this type of classifier information is
% present (if we were using only a subset of the data, then we would have
% to calculate multiple maps and average them together - this is of course
% possible, as well, and theoretically even better).

% Set defaults
cfg = decoding_defaults;

% Set the analysis that should be performed (default is 'searchlight', we want 'wholebrain' or 'roi' for our weight maps)
cfg.analysis = 'wholebrain';

% Set the output directory where data will be saved, e.g. 'c:\exp\results\buttonpress'
cfg.results.dir = 

% Set the filepath where your SPM.mat and all related betas are, e.g. 'c:\exp\glm\model_button'
beta_loc = 

% Set the filename of your brain mask (or your ROI masks as cell matrix) 
% for searchlight or wholebrain e.g. 'c:\exp\glm\model_button\mask.img' OR 
% for ROI e.g. {'c:\exp\roi\roimaskleft.img', 'c:\exp\roi\roimaskright.img'}
% You can also use a mask file with multiple masks inside that are
% separated by different integer values (a "multi-mask")
cfg.files.mask = 

% Set the label names to the regressor names which you want to use for 
% decoding, e.g. 'button left' and 'button right'
% don't remember the names? -> run display_regressor_names(beta_loc)
labelname1 = 
labelname2 = 

%% Set additional parameters
% Set additional parameters manually if you want (see decoding.m or
% decoding_defaults.m). We need to set the output to the reconstructed
% pattern, i.e. weights transformed according to Haufe et al. (2014)

% cfg.results.output = {'accuracy_minus_chance','AUC_minus_chance'};

% cfg.verbose = 2; % you want all information to be printed on screen
% cfg.decoding.train.classification.model_parameters = '-s 0 -t 0 -c 1 -b 0 -q'; 

% Some other cool stuff
% Check out 
%   combine_designs(cfg, cfg2)
% if you like to combine multiple designs in one cfg.


%% Add additional output measures if you like
% See help decoding_transform_results for possible measures

% cfg.results.output = {'accuracy_minus_chance', 'AUC'}; % 'accuracy_minus_chance' by default



%% Nothing needs to be changed below for a standard leave-one-run out cross
%% validation analysis.

% The following function extracts all beta names and corresponding run
% numbers from the SPM.mat
regressor_names = design_from_spm(beta_loc);

% Extract all information for the cfg.files structure (labels will be [1 -1] )
cfg = decoding_describe_data(cfg,{labelname1 labelname2},[1 -1],regressor_names,beta_loc);

% This creates the leave-one-run-out cross validation design:
cfg.design = make_design_cv(cfg); 

% Run decoding
results = decoding(cfg);

%% Now inspect the results and if they look reasonable, run the analysis again for the reconstructed pattern

cfg.results.output = 'SVM_pattern'; % if you really want the weights, use 'SVM_weights'

cfg.design = make_design_alldata(cfg); % overwrite the existing design
cfg.design.nonindependence = 'ok'; % set non-independence to be ok

results = decoding(cfg); % re-run decoding analysis

% you can potentially devise an analysis that uses permutation testing on
% permuted labels for revealing parts of the pattern that are reliable


