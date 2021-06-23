function output = es_removeBetween_rotated(input,rotations)
% Remove between subjects variability from data matrix (nsubjects X nconditions or nsubjects X nconditions X nmeasures)
% for plotting suitable error bars for within-subject designs.
% Third "nmeasures" dimension useful for high-dimensional data such as EEG where conditions have multiple timepoints of data.
% See Loftus and Masson (1994)
%
% Written by Ed Sohoglu 2011. Modified 2014.

input = permute(input,rotations); %Permute to match expected dimensions

nconditions = size(input,2);
nsubjects = size(input,1);

means_all_conditions = nanmean(input,2);
grand_mean = nanmean(means_all_conditions,1);
grand_mean = repmat(grand_mean,nsubjects,1);
adjustment_factors = grand_mean - means_all_conditions;
adjustment_factors = repmat(adjustment_factors,1,nconditions);
output = input + adjustment_factors;
output = ipermute(output,rotations); %Permute back to input dimensions