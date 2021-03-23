% function output = transres_accuracy_pairwise_minus_chance(decoding_out, chancelevel, varargin)
%
% For more than two classes, rather than getting the mean multiclass
% accuracy across all classes (i.e. chance = 1/n_label), report the mean
% accuracy of all pairwise comparisons (i.e. chance = 1/2).
%
% The function uses transres_accuracy_pairwise, so see there for more
% information on output etc.
%
% NOTE: This code does not return the accuracy for each pair, but the mean
% accuraccy across all pairs. If you want the accuracy for each pair, use
% transres_accuracy_matrix. It is also no confusion matrix (see 
% transres_confusion_matrix for that).
%
% OUT
% The output will be one number: the average accuracy across all possible 
% pairwise comparisons minus the pairwise chancelevel (50 percent).
%
% This code runs faster if all labels are in the same order in all decoding
% steps (e.g. runs).
%
% To use this transformation, use
%
%   cfg.results.output = {'accuracy_pairwise'}
%
% Martin Hebart 2016-03-09
%
% See also transres_accuracy_matrix transres_accuracy_pairwise transres_confusion_matrix
%   transres_accuracy_matrix_minus_chance

% Hist: Kai, 2020-06-17: introduced chance level, now correct

% TODO: allow using subset of accuracy matrix

function output = transres_accuracy_pairwise_minus_chance(decoding_out, chancelevel, varargin)

% use new way to also return chancelevel
chancelevel = 50; % reset chancelevel: for pairwise accuracies always 50 percent
output.output = transres_accuracy_pairwise(decoding_out,chancelevel,varargin{:});
output.output = output.output - chancelevel; % because our chancelevel is always 50 %
output.chancelevel = chancelevel;