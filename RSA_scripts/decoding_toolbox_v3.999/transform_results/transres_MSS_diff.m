% function output = transres_MSS_diff(decoding_out, chancelevel, varargin)
%
% Get the sum of squares difference between predicted and true labels. This
% function is particularly useful for regression approaches, but can also
% be used meaningfully for classification when neighboring classes are also
% more similar to each other. Smaller difference is better, but the only
% way to assess the quality of the final result is a permutation test.
%
% To use this transformation, use
%
%   cfg.results.output = {'MSS_diff'}
%
% Martin Hebart 2018-07-26
%
% See also decoding_transform_results

function output = transres_MSS_diff(decoding_out, chancelevel, varargin)

% run comparison for each step separately and average
n_steps = length(decoding_out);
output_sep = zeros(1,n_steps);
for i_step = 1:n_steps
    x = (decoding_out(i_step).predicted_labels-decoding_out(i_step).true_labels); % difference
    output_sep(i_step) = 1/size(x,1) * (x'*x); % mean sum of squares difference
end
output = output_sep/n_steps;