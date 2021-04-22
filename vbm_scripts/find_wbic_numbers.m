cd('/imaging/tc02/vespa/scans/PNFA_VBM/tom')
load('image_paths.mat')
for i = 1:45
test = strsplit(all_imagepaths{i},'mprage_125.nii')
thisnum(i) = str2num(test{1}(end-5:end-1))
end
thisnum = thisnum'