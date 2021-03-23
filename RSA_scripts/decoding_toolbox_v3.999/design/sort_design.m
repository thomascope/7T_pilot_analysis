% function [design, sortind, sortind_inv] = sort_design(design, sortind)
%
% This function is used to sort designs in a way to speed up decoding
% analyses. If training data happen to be identical in two of the many
% cross-validation iterations, then re-training of data is not necessary.
% decoding.m can recognize this and skip repeated training.
%
% Input:
%   design: with fields train, test, label and set
%   sortind (optional): If sorting should be done manually. This may be
%       useful e.g.to revert the sorting process later if requested.
%
% Output:
%   design: sorted
%   sortind: index used for sorting (if interesting)
%   sortind_inv: index necessary to invert sorting to original (if interesting)

% by Martin Hebart

% Martin 2018/01/23 Fixed bug that prevented this function from working
%
% Kai, 2017/03/28 Using sortrows() instead of sort(), no idea what it did
%   before. Also added new sorting using the first occurence of a column.
% 
%   

function [design, sortind, sortind_inv] = sort_design(design, sortind)

tr = design.train;

%% Get a sorting order if none is provided
if ~exist('sortind','var')

    % Option 1: Sort by first occurence
    [sortind, sortind_inv] = sort_columns_by_first_occurence(tr);

    % Option 2: use standard sortrows
%     [dummy, sortind] = sortrows(tr'); % sort the columns of the training matrix (which equals the rows of the transposed matrix)
%         % we dont keep the sorted training matrix (here dummy) but only the sortind to make sure that we resort all design fields in the same way, including the training set
%     sortind = sortind(end:-1:1); % reverse sorting order, will typically create training matrices that look more familiar to our standard designs because they start in the upper left corner, not in the lower left

    % Original option:
    % that did not work (remove at some point if its not useful)
    % Martins first version that for some reason does not work on Kais computer
    % If sortrows below does not work in newer Matlabversions, than that might
    % be the difference. 
    %     sortind = 1:size(tr,2);
    %     for i = size(tr,1):-1:1
    %         [ignore,subind] = sort(tr(end+1-i,:));
    %         sortind = sortind(subind);
    %     end

    %% check if there are any repetitions
    trcheck = tr(:,sortind);
    d = diff(trcheck,1,2);
    if all(sum(abs(d)))
        sortind = 1:size(tr,2); % if no repetitions exist, just use original index
        dispv(1,'sort_design: Did not re-sort design, because no training set occurs multiple times, thus no additional speed-up through sorting possible...')
    else
        dispv(1,'sort_design: Re-sorting design for additional speed-up...')
    end

end

%% Check that sortind will resort all columns once
% verify that all columns where assigned exactly once
if ~isequal(sort(sortind), 1:size(tr, 2))
    error('Sortind not valid because it does not contain each column number of the design once.')
end


%% Apply sortind to all other relevant fields
design.train = design.train(:,sortind);
design.test  = design.test(:,sortind);
design.label = design.label(:,sortind);
design.set   = design.set(sortind);

%% Proactive check: Check if any other field exists with the same length and through a warning
fnames = fieldnames(design);
for fn_ind = 1:length(fnames)
    curr_fname = fnames{fn_ind};
    if ~ismember(curr_fname, {'train', 'test', 'label', 'set'}) && (length(design.(curr_fname)) == size(tr, 2) ||  size(design.(curr_fname), 2) == size(tr, 2))
        warningv(1, ['sort_design:design.' curr_fname '_not_resorted_but_has_same_length'], ['sort_design: The field design.' curr_fname ' was not resorted but might need resorting because it has the same number of elements. Possible fixes: 1. Resort this fields with subind, 2. Notify the developers to add it'])
    end
end

%% Generate 
reverse = 1:size(tr,2);
sortind_inv(sortind) = reverse;

end

%% Helper function [sortind, sortind_inv] = sort_columns_by_first_occurence(m)
%
% Will provide sorting order if columns should be sorted by first
% occurence.
% Look more familiar to the user than using sortrows, and is equally 
% efficient.
%
% Kai, 2017/03/28

function [sortind, sortind_inv] = sort_columns_by_first_occurence(m)

%% resort a matrix by unique columns, using the first occurence of each column

% subplot(3,1,1);
% imagesc(m)

% initialize sortind_inv with 0
sortind_inv = zeros(1, size(m, 2)); % sortind_inv specifies at then end where to move each column to in the sorted matrix
sort_counter = 0;

for col_ind = 1:size(m, 2)
    if sortind_inv(col_ind) ~= 0
        % nothing to do, this column has already been assigned a sortind_inv
    else
        % a new column, assign it the next number
        sort_counter = sort_counter + 1;
        sortind_inv(col_ind) = sort_counter;
%         display(sprintf('Col %i is new, assigning sortind_inv %i', col_ind, sort_counter))
        
        % get remaining unassigned columns
        unassigned_cols = find(sortind_inv==0);
        
        % check if any of these matches the current column
        for curr_unassigned_col = unassigned_cols
            if all(m(:, curr_unassigned_col) == m(:, col_ind))
                % column matches, add index
                sort_counter = sort_counter + 1;
                sortind_inv(curr_unassigned_col) = sort_counter;
%                 display(sprintf('Col %i matches col %i, assigning sortind_inv %i', curr_unassigned_col, col_ind, sort_counter))
            end
        end
    end
end

% get index to be applied for sorting (sort the inverse sorting)
[dummy, sortind] = sort(sortind_inv);
% verify that all columns where assigned exactly once
if ~isequal(sort(sortind), 1:size(m, 2))
    error('Sorting did not succeed, please check the code')
end

% %% 
% subplot(3,1,2);
% m_sorted = m(:, sortind);
% imagesc(m_sorted);
% 
% subplot(3,1,3);
% m_sorted_inv = m_sorted(:, sortind_inv);
% imagesc(m_sorted_inv);
end