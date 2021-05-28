function module_make_FSGD(subjects, group, group_names, age_lookup, datafolder, studyname)
% A script for making a Freesurfer FSGD file, defining the data for group
% analysis, then makes a contrast for group1 > group2, accounting for age
% group_names = {'Control','nfvPPA'}
% age_lookup = readtable('Pinfa_ages.csv');
% datafolder = [preprocessedpathstem '/freesurfer_skullstripped/'];
% studyname = 'PINFA'
FSGDdir = [datafolder '/FSGD/'];
if ~exist(FSGDdir,'dir')
    mkdir(FSGDdir)
end

jobfile = [FSGDdir studyname '.fsgd'];

fileID = fopen(jobfile,'w');

%Mandatory preamble
fprintf(fileID,['GroupDescriptorFile 1\n']);
fprintf(fileID,['Title ' studyname '\n']);

%Name groups
for this_group = 1:length(group_names)
fprintf(fileID,['Class ' group_names{this_group} '\n']);
end

fprintf(fileID,['Variables Age\n']);

%Now define subjects
for this_subj = 1:length(subjects)
    this_age = age_lookup.Age(strcmp(age_lookup.Study_ID,subjects{this_subj}));
    fprintf(fileID,['Input ' subjects{this_subj} ' ' group_names{group(this_subj)} ' ' num2str(this_age) '\n']);
end

fclose(fileID);

% Now make contrasts
Contrastdir = [datafolder '/Contrast/'];
if ~exist(Contrastdir,'dir')
    mkdir(Contrastdir)
end

jobfile = [Contrastdir group_names{1} '-' group_names{2} '.mtx'];
fileID = fopen(jobfile,'w');
% fprintf(fileID,['1 -1 0 0']); %for DODS
fprintf(fileID,['1 -1 0']); %for DOSS
fclose(fileID);