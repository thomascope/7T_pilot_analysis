% function output = transres_accuracy_matrix(decoding_out, chancelevel, varargin)
%
% Get a matrix of the accuracies of all pairwise comparisons (i.e. chance =
% 50 %) from a multiclass classification. Libsvm doesn't provide this
% directly. The diagonal is not defined and set to NaN.
% IMPORTANT:
% This function is NOT a confusion matrix and does not provide a matrix of multiclass
% accuracies (where chance = 1/n_class). For the confusion matrix use
% transres_confusion_matrix 
% Reporting the multiclass accuracy of each class is not implemented yet.
%
%
% OUT
% The output will be an NxN matrix where n is the number of unique labels.
% The accuracucy for classifiying label1 vs label2 is in 
% acc(label1, label2). Labels are in increasing order.
%
% NOTE: If labels are missing, e.g. if only the labels [1, 2, 4] exist, the
% value for classifying 2 vs 4 is in acc(2, 3) (or in acc(3,2) as
% classification is symmetric), not in acc(2, 4) (or acc(4,2)).
%
% This code runs faster if all labels are in the same order in all decoding
% steps (e.g. runs).
%
% To use this transformation, use
%
%   cfg.results.output = {'accuracy_matrix'}
%
% Martin Hebart 2016-03-09
% Kai, 2020-06-17: Major update to canonical dv format from libsvm_test.m
%
% See also decoding_transform_results 
%   transres_accuracy_pairwise transres_confusion_matrix
%   transres_accuracy_pairwise_minus_chance transres_accuracy_matrix_minus_chance

% TODO: allow using subset of accuracy matrix

function output = transres_accuracy_matrix(decoding_out, chancelevel, varargin)

% A more general form of this code would be to allow using only a subset
% of the matrix in each iteration (e.g. have only some test labels). This
% is not as easy as it may sound. The most elegant solution would probably
% be to create an accuracy matrix of size n_label x n_label x n_step and in
% each iteration use NaNs for the missing values, and average across the
% matrices using our own version of nanmean.
%
% In addition, we would have to use a different M matrix for each run
% (which should not add much computational overhead since M is a persistent
% variable anyway (same applies to occurrence and keepind)
%
% We would also have to be able to relate the original labels to the
% subindex of the labels that we have in a current iteration.

persistent M occurrence keepind persistentlabels 

n_step = length(decoding_out);
all_labels = vertcat(decoding_out.true_labels);
all_ulabels = uniqueq(all_labels);
n_label = size(all_ulabels,1);

%% Create index matrix and keep it during iterations (updated when labels change)
if ~isequal(persistentlabels, all_labels) % if this is the first iteration or there was any change
    persistentlabels = all_labels; % remember labels to avoid recomputation
    
    % we will run a check if all training labels in all runs are the same and have the same identity and count (not necessarily order) as the test labels
    test_labels = decoding_out(1).true_labels;
    if ~isequal(all_ulabels,uniqueq(test_labels))
        error('Some runs have different labels than others. transres_accuracy_matrix cannot deal with this case yet (sorry, it''s really difficult to code, see code for a detailed description). Please set up a design that does all pairwise classifications and use output sensitivity and specificity.')
    end
    prev = sort(decoding_out(1).model.Label);
    for i_step = 2:length(decoding_out)
        s = sort(decoding_out(i_step).model.Label);
        test_labels = decoding_out(i_step).true_labels;
        if (length(prev) ~= length(s)) || any(prev ~= s) || (n_label == length(s) && any(all_ulabels ~= s))
            error('Number and/or identity of training labels and/or test labels is not the same across all steps. Unfortunately transres_accuracy_matrix cannot deal with this case yet.  Please set up a design that does all pairwise classifications and use output sensitivity and specificity.')
        end
        if ~isequal(all_ulabels,uniqueq(test_labels))
            error('Some runs have different labels than others. transres_accuracy_matrix cannot deal with this case yet (sorry, it''s really difficult to code, see code for a detailed description). Please set up a design that does all pairwise classifications and use output sensitivity and specificity.')
        end
    end
    
    % For speeding everything up, we will use a matrix that we define here that
    % extracts accuracies from the multiclass classification decision values
    M = zeros(n_label,(n_label*(n_label-1))/2); % init with nchoosek(n_label,2)
    
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
    occurrence = logical(M); % we need this for indexing (M: cols in ascendig label order)
    
end

%% Get pairwise accuracies from decision value matrix dv
% Note: expecting dv in canonical form, see e.g. tdt function libsvm_test.m
acc = zeros(n_label,n_label);
a = zeros(n_label,n_label);

% loop over all cross-validation iterations
for i_step = 1:n_step
    
    test_labels = decoding_out(i_step).true_labels;
    
    % get the unique labels for the model
    % note: order of labels is not necessary anymore, because dv are 
    % ordered and brought to canonical form by libsvm_test.m
    modellabel = sort(decoding_out(i_step).model.Label); % Note: why sort(): libsvm values are not ordered, but decision values ARE now ordered by tdt function libsvm_test.m, so we can order the labels here

    % get the decision value matrix dv. 
    % dimension of dv: ntestsamples x npairwiseclassifiers.
    % dv contains in each column the decision value of a binary classifier, 
    % i.e. has as many columns as pairwise comparisions are possible.
    % E.g. for 3 classes with labels [1, 2, 4] three columns exist, for 
    % classifiers that classify 1vs2, 1vs4, 2vs4 (note: each classifiers
    % _test.m function should take care that the dv has this order, see
    % e.g. libsvm_test.m).
    % Each row then contains the decision value for each of these 
    % classifiers for each test sample, of course independent of the true 
    % label of the test sample.
    % Different options exist to get a class decision from these pairs,
    % which is what the dv matrix is normally used for. Here we use it
    % however to get a pairwise accuracy matrix from this.
   
    dv  = decoding_out(i_step).decision_values;
    
    if isequal(all_ulabels(:),test_labels(:)) % each label occurs only once and in the correct order
        % dv now sorted by libsvm_test.m, not necessary to check label order anymore
        acc(keepind) = acc(keepind) + (sum(dv.*M>0)./sum(occurrence))';
    else
        % get the index of occurrence of all instances
        [~,idx] = ismember(test_labels,sort(modellabel)); % tdt function libsvm_test.m brings dv in canonical form now
        
        % create an index for later sorting of rows and columns
        [~,sortind] = sort(modellabel); % labels should be sorted now

        % This step calculates all accuracies in the matrix by multiplying
        % the DVs with a sign matrix which will make all correct DVs
        % positive, all incorrect ones negative and will set all irrelevant
        % ones to 0; it will then look what percentage was positive.
        % For that, it creates a new matrix which takes for each test label 
        % the row corresponding row of matrix M, creating a 
        % ntestsample x nclasses matrix
  
        a(keepind) = sum(dv.*M(idx,:)>0)./sum(occurrence(idx,:));

        % sum up with previous steps and sort labels to correct order (it
        % will look weird, but once we have added the transpose to the
        % matrix all is fine again)
        % Note: sorting labels is not necessary anymore, now sorted by tdt
        % function libsvm_test.m
        acc = acc+a(sortind,sortind);
        
        % code to display the individual values for bugtracking purposes
        %         idx
        %         disp('M(idx,:)')
        %         M(idx,:)
        %         disp('dv.*M(idx,:)')
        %         dv.*M(idx,:)
        %         disp('sum(dv.*M(idx,:)>0)')
        %         sum(dv.*M(idx,:)>0)
        %         disp('occurrence(idx,:)')
        %         occurrence(idx,:)
        %         sum(occurrence(idx,:))
        %         a
    end
    
end

accuracy = acc/n_step; % get mean accuracy by dividing sum by number of occurrences and arrange in right format
output = {100*(accuracy + accuracy') + diag(nan(1,n_label))}; % add transpose (accuracies are symmetrical) and convert to percent