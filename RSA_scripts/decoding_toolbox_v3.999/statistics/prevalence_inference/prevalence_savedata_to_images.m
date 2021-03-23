% prevalence_savedata_to_images(outputfilename, mask, vol, results, result_indcs, source_and_ref)
%
% TDT helper function to save data from the prevalence analysis to images. 
% Will create one image for each field of results (see prevalenceCore.m).
% Called by prevalence.m from TDT.
%
% IN
%   outputfilename: Filename with full path to be put before the three
%       images.
%   mask: Volume with true where data should be written to. Will be written
%       to _mask.nii.
%   vol: should contain a 4x4 transformation/rotation matrix. 
%   results: Struct, output of prevalenceCore.m. One file will be written 
%       for each field of results, e.g. _gamma0 for the gamma0 parameters.
%       For more description of the information, see prevalenceCore.m
% OPT
%   result_indcs: Indices that are taken from all fields in results, e.g.
%       to save ROI images for different ROIs.
%   source_and_ref: what the written images should have in their
%       description where the image comes from, e.g.
%           ' from prevalence inference, Allefeld et al 2016, Neuroimage'
% OUT
%   Written files, see above.
%
% Kai, 2016/08/01 (adapted to core function)

function prevalence_savedata_to_images(outputfilename, mask, vol, results, result_indcs, source_and_ref)

disp(['Saving images to ' outputfilename])
fnames = fieldnames(results);

for fname_ind = 1:length(fnames)
    curr_fname = fnames{fname_ind};
    curr_outputfile = [outputfilename '_' curr_fname '.nii'];
    
    % take all indices if no indices have been provided
    if ~exist('result_indcs', 'var') || isempty(result_indcs)
        result_indcs = 1:length(results.(curr_fname));
    end
    
        % create output image
    data = nan(size(mask));
    data(mask) = results.(curr_fname)(result_indcs);
    disp(['  Saving ' curr_outputfile])
    if ~exist('source_and_ref', 'var')
        source_and_ref = ' from prevalence inference, Allefeld et al 2016, Neuroimage'; % default for backward compatibility reasons
    end
    saveMRImage(data, curr_outputfile, vol.mat, [curr_fname, source_and_ref])
end

% also save the mask
saveMRImage(uint8(mask), [outputfilename '_mask.nii'], vol.mat, 'prevalence map mask')