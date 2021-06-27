function word = featureToWord_reverse(feature, feature_mat_words)
PRINT_DIAGNOSTICS = 0;

if PRINT_DIAGNOSTICS
    feature
end

% Transpose multiplication / correlation without normalization
word = feature * inv(feature_mat_words);

if PRINT_DIAGNOSTICS
    word 
end
end
