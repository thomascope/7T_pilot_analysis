function [countLoopPE, PE_feature_Accumulated, posterior_word_Iterative, residual_noise] = ...
    predictiveCodingLoop_Sohoglumodel(sensory_inputProb, prior_word_Iterative, feature_mat_words, ...
    prior_update_weight, STOPcriterion)

% initialize values for iterative loop
nWords      = size(feature_mat_words, 1);
nFeatures      = size(feature_mat_words, 2);
countLoopPE = 0;
sumPE       = 10;
PE_feature_Accumulated = zeros(1,nFeatures);

prior_features_Iterative = wordToFeature(prior_word_Iterative, feature_mat_words);

% estimate precisions of prior and sensory input
precision_word_prior = std(prior_word_Iterative/sum(prior_word_Iterative));
precision_sensory    = std(sensory_inputProb/sum(sensory_inputProb));

while sumPE > STOPcriterion
    
    % compute prediction error for features (here 37 dimensions for full set)
    PE_feature = sensory_inputProb - prior_features_Iterative;
    
    % compute precition of the Prediction Error
    % Precision of PE depends on the precisions of its constituents:
    % precision_sensory = low + precision_word_prior = low
    %     => precision PE = low
    % precision_sensory = high + precision_word_prior = low (or reverse)
    %     => precision PE = medium
    % precision_sensory = high + precision_word_prior = high
    %     => precision PE = high
    precision_PE = precision_sensory + precision_word_prior;
    
    prior_features_Iterative = prior_features_Iterative+(PE_feature*prior_update_weight*precision_PE);
    
    % get prediction error on word level from first iteration
    if countLoopPE < 1
        PE_feature_Accumulated = PE_feature_Accumulated + abs(PE_feature);
    end
    
    % check stop criterion for iteration
    sumPE = sum(abs(PE_feature));
    countLoopPE = countLoopPE + 1;
    store_sumPE(countLoopPE) = sumPE;
%     stop loop, after x iterations
    if countLoopPE >= 500
        break
    end
end
[residual_noise, posterior_word_Iterative] = min(sum(abs(feature_mat_words - prior_features_Iterative),2));
end
