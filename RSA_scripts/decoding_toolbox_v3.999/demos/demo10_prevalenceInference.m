% This function tells you the location for demos for prevalence_inference.
% You should find them in decoding_toolbox/statistics/prevalence_inference

if isempty(which('decoding_defaults')), error('Please add TDT'), end

disp(['See demos in' newline ...
    '   ' fileparts(which('prevalenceTDT')) newline ...
    'on how to use the prevalence inference analysis'])
