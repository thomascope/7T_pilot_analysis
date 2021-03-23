% output = transres_other_meandist(decoding_out, varargin)
% 
% This function averages other output that was generated across decoding
% steps (in the field decoding_out.opt), but works only if this other
% output is numerical. This is useful if you have a similarity matrix per
% cross-validation step and want to average across those
%
% To use it, use
%
%   cfg.results.output = {'other_meandist'}
%
% Martin, 2017-06-30

function output = transres_other_meandist(decoding_out, varargin)

currdim = ndims(decoding_out(1).opt)+1;

% average across the optional output and across entries in the dissimilarity matrix
output = {mean(squareformq(mean(cat(currdim, decoding_out.opt),currdim)))};