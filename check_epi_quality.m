%Now manually review data quality
all_masks = {};
for i = 1:length(subjects)
    all_masks{i} = ['/imaging/mlr/users/tc02/PINFA_preprocessed_2021/' subjects{i} '/stats3_3/mask.nii,1'];
end
spm_check_registration(char(all_masks))

%Manually review, list subjects with deficient masks:
Bad_subj_masks = {'P7P01'
    'P7P16'
    'P7C17'
    'P7C18'};

%Now manually review each run quality
all_badsubj_epis = {};
for i = 1:length(subjects)
    if ~strcmp(subjects{i},Bad_subj_masks)
        continue
    else
        theseepis = find(strncmp(blocksout{i},'Run',3));
        for sess = 1:length(theseepis)
            all_badsubj_epis{end+1} = ['/imaging/mlr/users/tc02/PINFA_preprocessed_2021/' subjects{i} '/s3rtopup_' blocksin{i}{theseepis(sess)} ',1'];
        end
    end
end
spm_check_registration(char(all_badsubj_epis))

%Manually review, list bad runs
Bad_run_subj_pairs = {'P7P16' , 4;
'P7C17', 4
'P7C18', 3}; %P7P01 has a frontal dropout on all scans, but probably outside of area of interest, so let him stand.