function module_make_FSGD_covariate(subjects, group, group_names, age_lookup, datafolder, studyname,covariates)
% A script for making a Freesurfer FSGD file, defining the data for group
% analysis, then makes a contrast for group1 > group2, accounting for age
% group_names = {'Control','nfvPPA'}
% age_lookup = readtable('Pinfa_ages.csv');
% datafolder = [preprocessedpathstem '/freesurfer_skullstripped/'];
% studyname = 'PINFA'

assert(length(subjects)==size(covariates,1),'Must be one covariate row per subject')

FSGDdir = [datafolder '/FSGD/'];
if ~exist(FSGDdir,'dir')
    mkdir(FSGDdir)
end

jobfile = [FSGDdir studyname '_grouped.fsgd'];

fileID = fopen(jobfile,'w');

%Mandatory preamble
fprintf(fileID,['GroupDescriptorFile 1\n']);
fprintf(fileID,['Title ' studyname '\n']);

%Name groups
for this_group = 1:length(group_names)
    fprintf(fileID,['Class ' group_names{this_group} '\n']);
end

%Now define covariates
hasnans = zeros(1,length(subjects));
this_covariate_string = {};
covariate_names = [];
for i = 1:size(covariates,2)
    covariate_names = [covariate_names ' Covariate_' num2str(i)];
    for j = 1:size(covariates,1)
        if i == 1
            this_covariate_string{j} = num2str(covariates(j,i));
        else
            this_covariate_string{j} = [this_covariate_string{j} ' ' num2str(covariates(j,i))];
        end
        if isnan(covariates(j,i))
            hasnans(j) = 1; % FSGD cannot have nans
        end
    end
end

fprintf(fileID,['Variables Age' covariate_names '\n']);

%Now define subjects
for this_subj = 1:length(subjects)
    this_age = age_lookup.Age(strcmp(age_lookup.Study_ID,subjects{this_subj}));
    if ~hasnans(this_subj)
        fprintf(fileID,['Input ' subjects{this_subj} ' ' group_names{group(this_subj)} ' ' num2str(this_age) ' ' this_covariate_string{this_subj} '\n']);
    end
end

fclose(fileID);

% Now make contrasts
Contrastdir = [datafolder '/Contrast/'];
if ~exist(Contrastdir,'dir')
    mkdir(Contrastdir)
end

for i = 1:size(covariates,2)
    jobfile = [Contrastdir 'Grouped_Covariate_' num2str(i) '.mtx'];
    fileID = fopen(jobfile,'w');
    % fprintf(fileID,['0 0 0 0 1 1']); %for DODS
    fprintf(fileID,['0 0 0 -1']); %for DOSS
    fclose(fileID);
end

%% Now ungrouped
jobfile = [FSGDdir studyname '_ungrouped.fsgd'];

fileID = fopen(jobfile,'w');

%Mandatory preamble
fprintf(fileID,['GroupDescriptorFile 1\n']);
fprintf(fileID,['Title ' studyname '\n']);

%Name groups - here gender
fprintf(fileID,['Class M\n']);
fprintf(fileID,['Class F\n']);

%Now define covariates
this_covariate_string = {};
covariate_names = [];
for i = 1:size(covariates,2)
    covariate_names = [covariate_names ' Covariate_' num2str(i)];
    for j = 1:size(covariates,1)
        if i == 1
            this_covariate_string{j} = num2str(covariates(j,i));
        else
            this_covariate_string{j} = [this_covariate_string{j} ' ' num2str(covariates(j,i))];
        end
        if isnan(covariates(j,i))
            hasnans(j) = 1;
        end
    end
end

fprintf(fileID,['Variables Age' covariate_names '\n']);

%Now define subjects
for this_subj = 1:length(subjects)
    this_age = age_lookup.Age(strcmp(age_lookup.Study_ID,subjects{this_subj}));
    this_sex = age_lookup.Sex{strcmp(age_lookup.Study_ID,subjects{this_subj})};
    if ~hasnans(this_subj)
        fprintf(fileID,['Input ' subjects{this_subj} ' ' this_sex ' ' num2str(this_age) ' ' this_covariate_string{this_subj} '\n']);
    end
end

fclose(fileID);

% Now make contrasts
Contrastdir = [datafolder '/Contrast/'];
if ~exist(Contrastdir,'dir')
    mkdir(Contrastdir)
end

for i = 1:size(covariates,2)
    jobfile = [Contrastdir 'Ungrouped_Covariate_' num2str(i) '.mtx'];
    fileID = fopen(jobfile,'w');
    fprintf(fileID,['0 0 0 -1']); %for DOSS
    fclose(fileID);
end

%% Now patients only
jobfile = [FSGDdir studyname '_patients.fsgd'];

fileID = fopen(jobfile,'w');

%Mandatory preamble
fprintf(fileID,['GroupDescriptorFile 1\n']);
fprintf(fileID,['Title ' studyname '\n']);

%Name groups - here gender
fprintf(fileID,['Class M\n']);
fprintf(fileID,['Class F\n']);

%Now define covariates
this_covariate_string = {};
covariate_names = [];
for i = 1:size(covariates,2)
    covariate_names = [covariate_names ' Covariate_' num2str(i)];
    for j = 1:size(covariates,1)
        if i == 1
            this_covariate_string{j} = num2str(covariates(j,i));
        else
            this_covariate_string{j} = [this_covariate_string{j} ' ' num2str(covariates(j,i))];
        end
        if isnan(covariates(j,i))
            hasnans(j) = 1;
        end
    end
end

fprintf(fileID,['Variables Age' covariate_names '\n']);

%Now define subjects
for this_subj = 1:length(subjects)
    this_age = age_lookup.Age(strcmp(age_lookup.Study_ID,subjects{this_subj}));
    this_sex = age_lookup.Sex{strcmp(age_lookup.Study_ID,subjects{this_subj})};
    if ~hasnans(this_subj) && group(this_subj)==2
        fprintf(fileID,['Input ' subjects{this_subj} ' ' this_sex ' ' num2str(this_age) ' ' this_covariate_string{this_subj} '\n']);
    end
end

fclose(fileID);

% Now make contrasts
Contrastdir = [datafolder '/Contrast/'];
if ~exist(Contrastdir,'dir')
    mkdir(Contrastdir)
end

for i = 1:size(covariates,2)
    jobfile = [Contrastdir 'Patients_Covariate_' num2str(i) '.mtx'];
    fileID = fopen(jobfile,'w');
    fprintf(fileID,['0 0 0 -1']); %for DOSS
    fclose(fileID);
end