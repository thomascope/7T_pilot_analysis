function feature = wordToFeature_3d(word, feature_mat_words)
PRINT_DIAGNOSTICS = 0;

if PRINT_DIAGNOSTICS
    word
end

% Matrix multiplication
for i = 1:length(feature_mat_words)
feature{i} = word * feature_mat_words{i};
end


if PRINT_DIAGNOSTICS
    feature 
end
end