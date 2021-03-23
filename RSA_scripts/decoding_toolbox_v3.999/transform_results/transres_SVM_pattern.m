% output = transres_SVM_pattern(decoding_out, chancelevel, cfg, data)
% 
% Calculates the pattern according to Haufe et al (2014), Neuroimage. This
% is done by first getting the weights in source space (primal problem), if
% a linear SVM was used (for non-linear methods no weights can be
% calculated for the primal problem). The bias term is not needed for this.
% To use it, use
%
%   cfg.results.output = {'SVM_pattern'}
%
% 
% Caution: This function uses cfg.design, so it needs a design and assumes
% you are in the main analysis (and not in e.g. feature_selection). It
% further assumes that all input models are related to their decoding step
% (i.e. model{1} is from iteration 1, etc.)
%
% Note: The functions currently does NOT work for
%   cfg.decoding.method = 'regression'
%
% OUTPUT
%   1x1 cell array of cell arrays for each output(step), with the pattern
%   as a 1xn_features numeric output.
%
%   
% Martin, 2014-01-15

% Update Kai: 2020-01-17: Forbid regression, requires code change & testing
% Update MH: 2016-12-10: Fixed bug that would only allow running one model
% Update MH 2016-08-22: Introduced compatibility with more than two classes
% and with kernel methods

function output = transres_SVM_pattern(decoding_out, chancelevel, cfg, data)

if strcmp(cfg.decoding.method, 'regression')
    % if you change the code so regression works, please adapt the header, too
    if isfield(cfg, 'acknowledge_transres_SVM_pattern_for_regression_is_experimental') && ...
            cfg.acknowledge_transres_SVM_pattern_for_regression_is_experimental
        warning('transres_SVM_pattern:not_implemented_for_regression', 'transres_SVM_pattern:Pattern reconstruction has currently not been extensively tested for cfg.decoding.method=''regression''. The implementation is probably very similar to the current one, but we have not used it so far and libsvm does not return the labels of the data for regression as is done for classification, so the code needs modification and testing. Feel free to extensively test it and the code, test it, and then be so kind to send us the result :)');
    else
        error('transres_SVM_pattern:not_implemented_for_regression', 'transres_SVM_pattern:Pattern reconstruction has currently not been extensively tested for cfg.decoding.method=''regression''. The implementation is probably very similar to the current one, but we have not used it so far and libsvm does not return the labels of the data for regression as is done for classification, so the code needs modification and testing. Feel free to extensively test it and the code, test it, and then be so kind to send us the result :) You can set cfg.acknowledge_transres_SVM_pattern_for_regression_is_experimental=1 to agree to use it.');
    end
end

%% check that input data has not been changed without the user knowing it
check_datatrans(mfilename, cfg); 

%% Get weights (implementation from libsvm website)

w = transres_SVM_weights(decoding_out, chancelevel, cfg, data);

% Unpack model
model = [decoding_out.model];

n_models = length(model);
output{1} = cell(n_models,1);
for i_model = 1:n_models
    
    weights = w{1}{i_model};
    currlabel = cfg.design.label(:,i_model);

%% Get pattern    

    pattern = zeros(size(weights));
    
    ulabel = sort(model(i_model).Label);
    n_label = length(ulabel);
    
    ct = 0;
    for i_label = 1:n_label
        for j_label = i_label+1:n_label
            ct = ct+1;
            label_ind = currlabel == ulabel(i_label) | currlabel == ulabel(j_label);
            data_train = data(cfg.design.train(:, i_model) > 0 & label_ind, :);
            [n_samples, n_dim] = size(data_train);
            curr_weights = weights(:,ct);
            
            if n_dim^2<10^7 % if pattern doesn't have a very large number of voxels
                pattern(:,ct) = cov(data_train)*curr_weights / cov(curr_weights'*data_train'); % like cov(X)*W * inv(W'*X')
            else % else do row by row (not much slower, even if we chunk it no dramatic speed-up)
                warningv('TRANSRES_SVM_PATTERN:pattern_calculation_slow','Pattern is very large, so its estimation will be very slow (up to minutes)!')
                scale_param = cov(curr_weights'*data_train');
                pattern_unscaled = zeros(n_dim,1);
                fprintf(repmat(' ',1,20))
                backstr = repmat('\b',1,20);
                for i = 1:n_dim % now calculate columnwise
                    if i == 1 || ~mod(i,round(n_dim/50)) || i == n_dim
                        fprintf([backstr '%03.0f percent finished'],100*i/n_dim)
                    end
                    data_train(:,i) = data_train(:,i) - mean(data_train(:,i)); % remove mean columnwise
                    data_cov = (data_train(:,i)'*data_train)/(n_samples-1);
                    pattern_unscaled(i,1) = data_cov * curr_weights;
                end
                fprintf('\ndone.\n')
                pattern(:,ct) = pattern_unscaled / scale_param; % like cov(X)*W * inv(W'*X')
            end
        end
    end
    output{1}{i_model} = pattern;
end



