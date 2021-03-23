% This script is a template for people who are interested in carrying out
% between-group decoding of brain imaging data. This is the case where each
% of your subjects carries only one condition. If you have one group and
% have multiple conditions per subject you want to compare between subject,
% see decoding_template_between_subject.
% If you have multiple groups where each subject has multiple conditions,
% our suggestion is to carry out the analysis within group and then compare
% the results (e.g. compare accuracies) between groups. 
% Since you cannot make use of the automatic extraction of image names,
% labels and decoding chunks, you need to enter the image names, labels and
% chunks (i.e. what data belong together) separately.
% Important: make sure that all subjects are in the same space (e.g.
% spatially-normalized data), else the voxels are not comparable across
% participants.

% Make sure the decoding toolbox and your favorite software (SPM or AFNI)
% are on the Matlab path (e.g. addpath('/home/decoding_toolbox') )
addpath('$ADD FULL PATH TO TOOLBOX AS STRING OR MAKE THIS LINE A COMMENT IF IT IS ALREADY$')
addpath('$ADD FULL PATH TO TOOLBOX AS STRING OR MAKE THIS LINE A COMMENT IF IT IS ALREADY$')

% Set defaults
cfg = decoding_defaults;

% Set the analysis that should be performed (default is 'searchlight')
cfg.analysis = 'searchlight';

% Set the output directory where data will be saved, e.g. 'c:\exp\results\buttonpress'
cfg.results.dir = 

% Set the filename of your brain mask (or your ROI masks as cell array) 
% for searchlight or wholebrain e.g. 'c:\exp\glm\model_button\mask.img' OR 
% for ROI e.g. {'c:\exp\roi\roimaskleft.img', 'c:\exp\roi\roimaskright.img'}
% You can also use a mask file with multiple masks inside that are
% separated by different integer values (a "multi-mask")
cfg.files.mask = 

% Set the following field:
% Full path to file names (1xn cell array) (e.g.
% {'c:\exp\glm\model_button\im1.nii', 'c:\exp\glm\model_button\im2.nii', ... }
cfg.files.name =  
% and the other two fields if you use a make_design function (e.g. make_design_cv)
%
% (1) a nx1 vector to indicate what data you want to keep together for 
% cross-validation (typically only matched controls in between-group,
% because subjects are else independent). If you don't have an obvious way
% to create chunks, set all participants to 1 (e.g. cfg.files.chunk = ones(n,1) )
cfg.files.chunk =
%
% (2) any numbers as class labels, normally we use 1 and -1. Each file gets a
% label number (i.e. a nx1 vector)
cfg.files.label = 

% Set additional parameters manually if you want (see decoding.m or
% decoding_defaults.m). Below some example parameters that you might want 
% to use:

% cfg.searchlight.unit = 'mm';
% cfg.searchlight.radius = 12; % this will yield a searchlight radius of 12mm.
% cfg.searchlight.spherical = 1;
% cfg.verbose = 2; % you want all information to be printed on screen
% cfg.decoding.train.classification.model_parameters = '-s 0 -t 0 -c 1 -b 0 -q'; 

% Some other cool stuff
% Check out 
%   combine_designs(cfg, cfg2)
% if you like to combine multiple designs in one cfg.

% Decide whether you want to see the searchlight/ROI/... during decoding
cfg.plot_selected_voxels = 500; % 0: no plotting, 1: every step, 2: every second step, 100: every hundredth step...

% Add additional output measures if you like
% cfg.results.output = {'accuracy_minus_chance', 'AUC_minus_chance'}

% Assuming there are no matched controls between groups, the way in which
% data are split up is arbitrary. For that reason, you can repeatedly
% subsample from both groups, in this case 100 times. This creates the
% leave-one-pair-out cross validation design with 100 decoding steps:
cfg.design = make_design_boot(cfg,100,1); % the 1 keeps test data balanced, too
% If you have a balanced design with multiple chunks (e.g. matched samples), use this function:
% cfg.design = make_design_cv(cfg);

% If you used a bootstrap design, then you might speed up processing using
% this function:
cfg.design = sort_design(cfg.design);

% Run decoding
results = decoding(cfg);