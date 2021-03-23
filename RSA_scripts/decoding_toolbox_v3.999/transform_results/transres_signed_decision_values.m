% function output = transres_signed_decision_values(decoding_out, chancelevel, varargin)
%
% Calculate signed decision values. This function outputs decision values
% of the classifier with negative sign for incorrect predictions and
% positive sign for correct predictions. The output can be understood as
% accuracy weighted by decision values, with an expected value of 0. In
% other words, it uses the size of the decision value as evidence weight.
% The more obviously correct samples will receive a higher weight than less
% obviously correct samples. This can be useful also when there are only
% very few test samples available to get a more continuous results measure.
%
% For more than 2 classes, this automatically implements the pairwise
% approach (multiclass seems awkward).
%
% To use this transformation, use
%
%   cfg.results.output = {'signed_decision_values'}
%
% Antonius Wiehler & Martin Hebart, 2015-09-03

function output = transres_signed_decision_values(decoding_out, chancelevel, varargin)

persistent M occurrence keepind

all_labels = uniqueq(vertcat(decoding_out.true_labels));

if all_labels == 2
    decision_values = vertcat(decoding_out.decision_values);  % decision values of classifier
    predicted_labels = vertcat(decoding_out.predicted_labels);  % predicted labels by classifier
    true_labels = vertcat(decoding_out.true_labels);  % observed labels
    % prepare predictions for multiplication (based on accuracy of each sample):
    correct_mult = 2 * (predicted_labels == true_labels) - 1;
    % wrong predictions will receive a negative weight (-1),
    % correct predictions a positive weight (1)
    % decision values are signed by correctness of prediction
    output = sum(abs(decision_values) .* correct_mult);
    return
end

n_step = length(decoding_out);

n_label = size(all_labels,1);

if size(keepind,1) ~= n_label % if this is the first iteration or there was any change
    
    % we will run a check if all training labels in all runs are the same and have the same identity and count (not necessarily order) as the test labels
    test_labels = decoding_out(1).true_labels;
    if ~isequal(all_labels,uniqueq(test_labels))
        error('Some runs have different labels than others. transres_accuracy_matrix cannot deal with this case yet (sorry, it''s really difficult to code, see code for a detailed description). Please set up a design that does all pairwise classifications and use output sensitivity and specificity.')
    end
    prev = sort(decoding_out(1).model.Label);
    for i_step = 2:length(decoding_out)
        s = sort(decoding_out(i_step).model.Label);
        test_labels = decoding_out(i_step).true_labels;
        if (length(prev) ~= length(s)) || any(prev ~= s) || (n_label == length(s) && any(all_labels ~= s))
            error('Number and/or identity of training labels and/or test labels is not the same across all steps. Unfortunately transres_accuracy_matrix cannot deal with this case yet.  Please set up a design that does all pairwise classifications and use output sensitivity and specificity.')
        end
        if ~isequal(all_labels,uniqueq(test_labels))
            error('Some runs have different labels than others. transres_accuracy_matrix cannot deal with this case yet (sorry, it''s really difficult to code, see code for a detailed description). Please set up a design that does all pairwise classifications and use output sensitivity and specificity.')
        end
        prev = s;
    end
    
    % For speeding everything up, we will use a matrix that we define here that
    % extracts accuracies from the multiclass classification decision values
    M = zeros(n_label,n_label*(n_label-1)/2); % init
    
    % we just need the position of 1, 2, etc. in the lower diagonal matrix,
    % with one row where each of them appears
    [x,y] = meshgrid(1:n_label);
    % pick lower diagonal
    keepind = tril(true(n_label),-1);
    ind = [x(keepind) y(keepind)];
    
    % Now we set each column to 1 where there is the respective number in the index
    for i = 1:n_label
        M(i,ind(:,1)==i) = 1;
        M(i,ind(:,2)==i) = -1;
    end
    
    % to get an idea what this is you can also figure, imagesc(M) and inspect each row
    occurrence = logical(M); % we need this for indexing
    
end

% loop over all cross-validation iterations
sdv = zeros(n_step,1);
for i_step = 1:n_step
    
    test_labels = decoding_out(i_step).true_labels;
    
    % get the unique labels and their order for the model
    ulabel = decoding_out(i_step).model.Label;
    % get the decision values
    dv  = decoding_out(i_step).decision_values;
    
    if isequal(all_labels(:),test_labels(:)) % each label occurs only once and in the correct order
        sdv(i_step) = sum((sum(dv.*M)./sum(occurrence)));
    else
        % get the index of occurrence of all instances
        [~,idx] = ismember(test_labels,ulabel);
        
        % This step calculates all accuracies in the matrix by multiplying
        % the DVs with a sign matrix which will make all correct DVs
        % positive, all incorrect ones negative and will set all irrelevant
        % ones to 0; it will then look what percentage was positive.
        % It respects the order of occurrence of each label and will
        % rearrange the matrix M accordingly (or increase its size if
        % there are multiple occurrences of each)
        sdv(i_step) = sum(sum(dv.*M(idx,:))./sum(occurrence(idx,:)));
        
    end
    
end

output = sum(sdv)/n_step; % get mean signed-decision value by dividing sum by number of occurrences and arrange in right format
    
    
