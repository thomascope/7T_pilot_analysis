% function output = decoding_transform_results(method,decoding_out,chancelevel,cfg,data)
%
% This function calculates a lot of different result measures defined by
% METHOD. It is called and applied by decoding_generate_output.m
%
% FURTHER METHODS: Type
%   >> transres_  
% and hit tab key
%
% e.g., it also calls external transres_XX functions that implement other
% methods, e.g. "trans_model_parameters" if method = 'model_parameter'.
%
% Please note that not all methods are useful for all approaches in the
% toolbox and that we don't test for all meaningless combinations. For
% example, for regression approaches (e.g. SVR), accuracy doesn't make much
% sense, but rather corr or zcorr should be used.
%
% METHODS IMPLEMENTED HERE:
% accuracy: decoding accuracy
% accuracy_minus_chance: decoding accuracy minus chance level (useful for 
%   SPM 2nd level)
% decision_values: 'raw' decision values for all input patterns as returned 
%   by the method
% predicted_labels: predicted labels for  all input patterns as returned 
%   by the method
% true_labels: true labels for  all input patterns as returned 
%   by the method
% sensitivity: accuracy of first label
% sensitivity_minus_chance: sensitivity minus chance level
% specificity: accuracy of second label
% specificity_minus_chance: specificity minus chance level
% balanced_accuracy: mean accuracy, calculated for all labels separately
%   (useful when bias present)
% balanced_accuracy_minus_chance: balanced_accuracy minus chance level
% dprime: z(hit rate) - z(false alarm rate)
% loglikelihood: measure of bias to one label: -1/2*(zHIT_rate^2 - zFA_rate^2)
% AUC: Area under the ROC (Receiver Operator Characteristics) Curve times 
%   100 (i.e. values from 0 to 100), built from classifier decision values, 
%   not from sensitivity/specificity (for more than 2 classes, average of
%   all pairwise comparisons is calculated)
% AUC_pairwise: Same result as AUC but for more than 2 classes
% AUC_minus_chance: like AUC, but minus chance level(useful for SPM 2nd level)
% AUC_matrix: for more than 2 classes all pairwise comparisons in a matrix
% corr: Correlation (useful e.g. for regression approaches, e.g. SVR)
% zcorr: Fisher-z-transformed correlation (necessary for averaging
%   correlations, e.g. across subjects)
%
% The function also allows adding new result transformation functions by 
% calling
%
%   output = transres_METHOD(decoding_out,chancelevel,cfg,model,data);
%
% where METHOD will be replaced by the provided method name.
% E.g., if you want to write your own result transformation function
% "yourmethod", the method should be named "transres_yourmethod", take
% the above inputs as input (or better varargin), and provide your desired
% output measure as output. If the method also should return a chancelevel,
% output should be a struct with output.chancelevel and output.output (see
% e.g. transres_accuracy_matrix_minus_chance).
%
% IN
%   method: desired method name as string (see above)
%   decoding_out: struct with result from last decoding step
%           Typically, decoding_out contains the fields
%           decoding_out.predicted_labels
%           decoding_out.true_labels
%       which are both 1 x n_step double vectors containing the predicted
%       and the true labels, so that these can be compared. However, in
%       principle output contains whatever the decoding method puts out
%       (e.g. if you write your own method).
%   cfg: the standard decoding cfg struct that was used for the last
%        decoding
%
% OUT
%   output: can be either a single number or a struct ({}) that can contain
%       any type of data. For nomal application, output contains the fields
%           output.predicted_labels
%           output.true_labels
%       which are both 1 x n_step double vectors containing the predicted
%       and the true labels, so that these can be compared. However, in
%       principle output contains whatever the decoding method puts out
%       (e.g. if you write your own method).

% TODO: to improve speed when multiple methods are calculated, pass
% structure "loaded" with fields predicted_labels, true_labels, labels,
% decision_values and introduce check in each method if fields exist
% (rather than using isfield which is slow inintialize
% loaded.isloaded.decision_values = 0; and change to 1 as soon as it is
% initialized. Problem: Structure of function doesn't work with setwise. 

function output = decoding_transform_results(method,decoding_out,chancelevel,cfg,data)

%% If method is a string
if strcmpi(method, 'accuracy') || strcmpi(method, 'accuracy_minus_chance')
    predicted_labels =  vertcat(decoding_out.predicted_labels);
    true_labels = vertcat(decoding_out.true_labels);
    
    output = 100 * (1/size(predicted_labels,1)) * sum(predicted_labels == true_labels); % calculate mean (faster than Matlab function)
    
    if strcmpi(method, 'accuracy_minus_chance')
        output = output - chancelevel; % subtract chancelevel from all output entries
    end

elseif strcmpi(method, 'decision_values')
    output = {vertcat(decoding_out.decision_values)};
    
elseif strcmpi(method, 'predicted_labels')
    output = {vertcat(decoding_out.predicted_labels)};
    
elseif strcmpi(method, 'true_labels')
    output.true_labels = cell(size(decoding_out));
    for step_ind = 1:length(decoding_out)
        output.true_labels{step_ind} = decoding_out(step_ind).true_labels;
    end  
    
elseif strcmpi(method, 'sensitivity') || strcmpi(method, 'sensitivity_minus_chance') % where the first label is correct
    predicted_labels =  vertcat(decoding_out.predicted_labels);
    true_labels = vertcat(decoding_out.true_labels);
    
    labels = uniqueq(true_labels);
    n_labels = size(labels,1);
    if n_labels > 2
        error('Too many labels for sensitivity measure! Check input labels.')
    end
    labelfilt = true_labels == labels(1); % use first label (only works with two labels)
    output = 100 * mean(predicted_labels(labelfilt) == true_labels(labelfilt));
    
    if strcmpi(method, 'sensitivity_minus_chance')
        output = output - chancelevel; % subtract chancelevel from all output entries
    end
    
elseif strcmpi(method, 'specificity') || strcmpi(method, 'specificity_minus_chance') % where the other label is correct
    predicted_labels =  vertcat(decoding_out.predicted_labels);
    true_labels = vertcat(decoding_out.true_labels);
    
    labels = uniqueq(true_labels);
    n_labels = size(labels,1);
    if n_labels > 2
        error('Too many labels for sensitivity measure! Check input labels.')
    end
    labelfilt = true_labels == labels(end); % use last label (only works with two labels)
    output = 100 * mean(predicted_labels(labelfilt) == true_labels(labelfilt));
    if strcmpi(method, 'specificity_minus_chance')
        output = output - chancelevel; % subtract chancelevel from all output entries
    end
    
elseif strcmpi(method, 'balanced_accuracy') || strcmpi(method, 'balanced_accuracy_minus_chance')
    predicted_labels =  vertcat(decoding_out.predicted_labels);
    true_labels = vertcat(decoding_out.true_labels);
    
    labels = uniqueq(true_labels);
    n_labels = size(labels,1);
    for i_label = 1:n_labels
        labelfilt = true_labels == labels(i_label);
        if i_label == 1
            output = 100 * mean(predicted_labels(labelfilt) == true_labels(labelfilt));
        else
            output = output + 100 * mean(predicted_labels(labelfilt) == true_labels(labelfilt));
        end
    end
    output = (1/n_labels) * output;

    if strcmpi(method, 'balanced_accuracy_minus_chance')
        output = output - chancelevel; % subtract chancelevel from all output entries
    end
    
elseif strcmpi(method, 'dprime')
    predicted_labels =  vertcat(decoding_out.predicted_labels);
    true_labels = vertcat(decoding_out.true_labels);
    
    output = dprimestats(true_labels,predicted_labels);
    
elseif strcmpi(method, 'loglikelihood')
    predicted_labels =  vertcat(decoding_out.predicted_labels);
    true_labels = vertcat(decoding_out.true_labels);
    
    [dprime,output] = dprimestats(true_labels,predicted_labels); %#ok<ASGLU>
    
elseif strcmpi(method, 'AUC') || strcmpi(method, 'AUC_minus_chance') ...
    || strcmpi(method, 'AUC_pairwise') || strcmpi(method, 'AUC_pairwise_minus_chance')
    
    decision_values = vertcat(decoding_out.decision_values);
    true_labels = vertcat(decoding_out.true_labels);
    labels = uniqueq(true_labels);
    
    if length(labels) == 2
        output = 100*AUCstats(decision_values,true_labels,labels,0); % express in percent
    elseif length(labels) > 2
        if strcmpi(method, 'AUC') || strcmpi(method, 'AUC_minus_chance')
            warningv('DECODING_TRANSFORM_RESULTS:ReportAUCpairwise','More than 2 labels for AUC. Running all pairwise comparisons and averaging (using AUCstats_pairwise.m).')
        end
        output = 100*mean(AUCstats_matrix(decision_values,true_labels,labels));
    end
        
    if strcmpi(method, 'AUC_minus_chance') || strcmpi(method, 'AUC_pairwise_minus_chance')
        output = output - chancelevel; % center around 0
    end
    
elseif strcmpi(method, 'AUC_matrix') || strcmpi(method, 'AUC_matrix_minus_chance')

    decision_values = vertcat(decoding_out.decision_values);
    true_labels = vertcat(decoding_out.true_labels);
    labels = uniqueq(true_labels);
    n_label = length(labels);
    
    temp = 100*AUCstats_matrix(decision_values,true_labels,labels);
    
    % transform vector to matrix
    output = zeros(n_label,n_label);
    output(logical(tril(ones(n_label,n_label),-1))) = temp;
    output = output + output' + diag(nan(1,n_label));
    
    if strcmpi(method,'AUC_matrix_minus_chance')
        output = output - chancelevel; % center around 0
    end
    
elseif strcmpi(method, 'corr')
       
    n_steps = length(decoding_out);
    output_sep = zeros(1,n_steps);
    for i_step = 1:n_steps
       output_sep(i_step) = correl(decoding_out(i_step).predicted_labels,decoding_out(i_step).true_labels);
       % check for potential error sources if any output is nan
       if any(isnan(output_sep(i_step)))
          if all(decoding_out(i_step).true_labels == decoding_out(i_step).true_labels(1)) % corr cannot work with 1 data point
             error('decoding_transform_results:corr_with_only_one_datapoint', 'decoding_transform_results.m tried to calculate the output measure cfg.results.output=''corr'' in a cv step in which has one test datapoint. Calculating a correlation requires at least two datapoints. Suggested solution: Make sure to have at least at least two samples in the test set of each cv step.'); 
          elseif all(decoding_out(i_step).true_labels == decoding_out(i_step).true_labels(1))
             error('decoding_transform_results:corr_with_all_true_labels_equal', 'decoding_transform_results.m tried to calculate the output measure cfg.results.output=''corr'' in a cv step in which all true values (labels) are equal. Calculating a correlation requires at least two different x-values (values of the independent variables). Suggested solution: Make sure to have at least at least two test samples with different true values/labels in each cv step.'); 
          elseif any(isnan(decoding_out(i_step).true_labels))
             error('decoding_transform_results:corr_with_true_label_nan', 'decoding_transform_results.m tried to calculate the output measure cfg.results.output=''corr'' but has at least one label that is nan. This does not work. Make sure the labels of your data are not none.');
          else
             warningv('decoding_transform_results:corr_returns_nan', 'decoding_transform_results.m tried to calculate the output measure cfg.results.output=''corr'' and returned nan. This might be because one part of your data (e.g. a voxel value) is nan. Other potential reasons were checked and are not the reason: more than 1 data point, different true labels (values of the independent variable), and no true label is nan. If nans occur in your data, you can ignore this warning');
          end
       end
    end
    output = tanh(mean(atanh(output_sep))); % z-transform and back to average correlation
    
elseif strcmpi(method, 'zcorr')

    n_steps = length(decoding_out);
    output_sep = zeros(1,n_steps);
    for i_step = 1:n_steps
       output_sep(i_step) = correl(decoding_out(i_step).predicted_labels,decoding_out(i_step).true_labels);
       % check for potential error sources if any output is nan
       if any(isnan(output_sep(i_step)))
          if all(decoding_out(i_step).true_labels == decoding_out(i_step).true_labels(1)) % corr cannot work with 1 data point
             error('decoding_transform_results:zcorr_with_only_one_datapoint', 'decoding_transform_results.m tried to calculate the output measure zcfg.results.output=''corr'' in a cv step in which has one test datapoint. Calculating a correlation requires at least two datapoints. Suggested solution: Make sure to have at least at least two samples in the test set of each cv step.'); 
          elseif all(decoding_out(i_step).true_labels == decoding_out(i_step).true_labels(1))
             error('decoding_transform_results:zcorr_with_all_true_labels_equal', 'decoding_transform_results.m tried to calculate the output measure zcfg.results.output=''corr'' in a cv step in which all true values (labels) are equal. Calculating a correlation requires at least two different x-values (values of the independent variables). Suggested solution: Make sure to have at least at least two test samples with different true values/labels in each cv step.');           
          else
             warningv('decoding_transform_results:zcorr_returns_nan', 'decoding_transform_results.m tried to calculate the output measure cfg.results.output=''zcorr'' and returned nan. This might be because one part of your data (e.g. a voxel value) is nan. Other potential reasons were checked and are not the reason: more than 1 data point, different true labels (values of the independent variable), and no true label is nan. If nans occur in your data, you can ignore this warning');
          end
       end
    end
    output = mean(atanh(output_sep)); % z-transform
    
elseif ischar(method) % all other methods
    
    fhandle = str2func(['transres_' method]);
    output = feval(fhandle,decoding_out,chancelevel,cfg,data);
    % e.g. if method = 'yourmethod', this calls:
    %  output = transres_yourmethod(decoding_out,chancelevel,cfg,data);
    
elseif isobject(method)
    % use passed handle directly and return
    output = method.apply(decoding_out,chancelevel,cfg,data);
    
else
    error('Dont know how to handle method %s', method)
    
end
    
end
