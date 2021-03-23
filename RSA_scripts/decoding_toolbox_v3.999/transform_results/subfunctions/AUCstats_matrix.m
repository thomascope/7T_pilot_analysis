% function AUC = AUCstats_matrix(decision_values, true_labels, labels)
% computes the area unter the curve (AUC) for the n-by-1 vector of decision
% values, corresponding n-by-1 vector of true labels. This is done for
% multiple pairwise comparisons. The order of assigned labels is given by
% the variable "labels".
%
% AUC - Area under the ROC curve (between 0 and 1), chance = 0.5

% 2016/06/01 Martin Hebart: rewrote from AUCstats to include multiple pairwise classification comparisons

function AUC = AUCstats_matrix(decision_values, true_labels, labels)

% sort labels
labels = sort(labels);

% sort values
[decision_values,ind] = sort(decision_values);
true_labels = true_labels(ind);

n_label = length(labels);

label_position = false(length(decision_values),n_label);
for i_label = 1:length(labels)
    label_position(:,i_label) = true_labels == labels(i_label);
end
n_ulabel = sum(label_position);
sensitivity = bsxfun(@times,1./n_ulabel,cumsum(label_position));

% now run all pairwise AUCs
AUC = zeros(nchoosek(n_label,2),1);
ct = 0; % init counter

for i_label = 1:n_label
    for j_label = i_label+1:n_label
        ct = ct+1;
        curr_labelposition = label_position(:,i_label) | label_position(:,j_label);
        curr_dv = decision_values(curr_labelposition);
        % handle ties and add [0 0] to start and [1 1] to end
        curr_sensitivity = [0 0; sensitivity(curr_labelposition,[i_label j_label]); 1 1];
        rmind = curr_dv(1:end-1) == curr_dv(2:end);
        curr_sensitivity(rmind,:) = [];
        
        % compute the area under the "curve"
        AUC(ct) = 1 - (1/2)*sum((curr_sensitivity(2:end,2) - curr_sensitivity(1:end-1,2)).*(curr_sensitivity(2:end,1) + curr_sensitivity(1:end-1,1)));
    end
end


