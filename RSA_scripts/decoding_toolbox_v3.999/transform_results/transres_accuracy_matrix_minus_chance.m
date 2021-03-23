% function output = transres_accuracy_matrix_minus_chance(decoding_out, chancelevel, varargin)
%
% Get a matrix of the accuracies of all pairwise comparisons (i.e. chance =
% 50 %). Will also return the used chancelevel.
%
% The function uses transres_accuracy_matrix, so see there for more
% information on output etc.
%
% IMPORTANT:
% This function is NOT a confusion matrix and does not provide a matrix of 
% multiclass accuracies (where chance = 1/n_class). For the confusion 
% matrix use transres_confusion_matrix 
%
% To use this transformation, use
%
%   cfg.results.output = {'accuracy_matrix_minus_chance'}
%
% Martin Hebart 2016-03-09
%
% See also decoding_transform_results 
%   transres_accuracy_matrix transres_accuracy_pairwise transres_confusion_matrix
%   transres_accuracy_pairwise_minus_chance 

% Hist: Kai, 2020-06-17: introduced chance level, now correct

function output = transres_accuracy_matrix_minus_chance(decoding_out, chancelevel, varargin)

% use new way to also return chancelevel
chancelevel = 50; % reset chancelevel: for pairwise accuracies always 50 percent
output.output = transres_accuracy_matrix(decoding_out,chancelevel,varargin{:});
output.output = {output.output{1} - chancelevel};
output.chancelevel = chancelevel;
