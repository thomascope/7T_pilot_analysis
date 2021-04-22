function show_VBM_files(prefix)

load('image_paths.mat')

try
split_stem = regexp(all_imagepaths, '/mprage_125/', 'split');

these_imagepaths = cell(size(all_imagepaths));
for crun = 1:length(all_imagepaths)
    these_imagepaths(crun) = cellstr([split_stem{crun}{1} '/mprage_125/' prefix split_stem{crun}{2}]);
    if ~exist(char(these_imagepaths(crun)),'file')
        these_imagepaths(crun) = cellstr([split_stem{crun}{1} '/mprage_125/' prefix split_stem{crun}{2}(1:end-4) '_Template.nii']); %Deal with the difference in naming convention depending if part of the dartel templating or not
    end
end

spm_check_registration(char(these_imagepaths))
catch
   disp('Warning, controls only')
   
   split_stem = regexp(Larger_control_full_imagepath, '/mprage_125/', 'split');
   
   these_imagepaths = cell(size(Larger_control_full_imagepath));
   for crun = 1:length(Larger_control_full_imagepath)
       these_imagepaths(crun) = cellstr([split_stem{crun}{1} '/mprage_125/' prefix split_stem{crun}{2}]);
       if ~exist(char(these_imagepaths(crun)),'file')
           these_imagepaths(crun) = cellstr([split_stem{crun}{1} '/mprage_125/' prefix split_stem{crun}{2}(1:end-4) '_Template.nii']); %Deal with the difference in naming convention depending if part of the dartel templating or not
       end
   end
   
   spm_check_registration(char(these_imagepaths))
end