function [sum_iterative_error, this_correlation,heard_correct_word, amount_above_average] = PINFA_blank_simplified(sigma_pred, threshold, sensory_detail)
% An adaptation of Helen Blank's 2016 computational modelling, integrated
% with Thomas Cope's 2017 Nature Comms Bayesian modelling, which was in
% turn based on Ed Sohoglu's Bayesian Prediction Error modelling.

condition{1} = 'A heard, Mismatch prior (B>A), Low sensory detail';
condition_code(1,:) = [1,1];
condition{2} = 'A heard, Mismatch prior (B>A), High sensory detail';
condition_code(2,:) = [1,2];
condition{3} = 'A heard, Matching prior (A>B), Low sensory detail';
condition_code(3,:) = [2,1];
condition{4} = 'A heard, Matching prior (A>B), High sensory detail';
condition_code(4,:) = [2,2];

%% Define experimental feature set
wordlist = {'bard','barge','lard','large','pit','pick','kit','kick','debt','deck','net','neck','robe','road','lobe','load'};

phonemelist = {'b#d','b#_','l#d','l#_','pIt','pIk','kIt','kIk','dEt','dEk','nEt','nEk','r5b','r5d','l5b','l5d'};
nWords = numel(phonemelist);

% transform words into features
feature_mat_words = PINFA_word2phoneme(phonemelist);
nFeatures = size(feature_mat_words, 2);

%% set priors in word space
nIterations = 20; %How many iterations is reasonable in the timescale? 
prior_update_weight = (sigma_pred^2)/100; %Prior variance divided by 100
% neutral = all words equally probable - not used here, but left for
% clarity of what we're trying to do and possible future model adaptation
prior_neutral_word = ones(nWords) * 1/nWords;
% Exactly 50/50 chance of a match - set up a matrix with 50% on diagonal
% and 50%/nWords on the off diagonal.
diagonal_value     = 0.5;
off_diagonal_value = (1-diagonal_value) / (nWords-1);
prior_match_word   = MatProb(nWords, diagonal_value, off_diagonal_value);
% MisMatch - Here there was a consistent relationship where the 8th next
% word was predicted so it's a circular shift. In other experimental
% paradigms this may not be the case.
prior_mismatch_word = circshift(prior_match_word,[8,0]);

%% set priors in feature space
prior_match_feature = wordToFeature(prior_match_word,feature_mat_words);
prior_mismatch_feature = wordToFeature(prior_mismatch_word,feature_mat_words);
prior_neutral_feature = wordToFeature(prior_neutral_word,feature_mat_words);
prior_match_feature = sqrt(prior_match_feature); % Otherwise always a very high prediction error in the match clear case - I think this is reasonable to assume there may be a squared relationship
prior_mismatch_feature = sqrt(prior_mismatch_feature);
prior_neutral_feature = sqrt(prior_neutral_feature);

all_priors{1} = prior_mismatch_feature;
all_priors{2} = prior_match_feature;
all_priors{3} = prior_neutral_feature;

rng('default') %Seed the random number generator for reproducibility
input_noise = abs(randn(size(prior_match_feature)));
input_noise = input_noise/max(max(input_noise)); %Scale 0-1
% 
% x = [-10:1:10];
% xdisplay = [-20 20];
PE_features = [];
sum_iterative_error = [];
posterior_features = [];
all_PE_features = [];
for i=1:length(condition)
    
    this_prior_features = all_priors{condition_code(i,1)};
    this_sensory_detail = sensory_detail(condition_code(i,2));
    
    for this_word = 1:nWords
        %this_input = (feature_mat_words(this_word,:)*this_sensory_detail)+(input_noise(this_word,:)*(1-this_sensory_detail)); %Assume that the noise represents a multiplication of the closed experimental set - reasonable I think, as one would not expect to confuse with never-heard features
        this_input = (feature_mat_words(this_word,:)*this_sensory_detail)+(input_noise(this_word,:)*(1-this_sensory_detail)); %Assume that the noise represents a multiplication of the closed experimental set - reasonable I think, as one would not expect to confuse with never-heard features
        PE_features(i,this_word,:) = this_input-this_prior_features(this_word,:);
        %Iterative loop - see how close we get in n iterations
        prior_features_updated = this_prior_features(this_word,:);
        for this_iteration = 1:nIterations
            this_PE_features = (this_input-prior_features_updated);
            if this_iteration ==1
                all_PE_features(i,this_word,:) = this_PE_features;
                sum_iterative_error(i,this_word) = sum(abs(this_PE_features));
            else
                all_PE_features(i,this_word,:) = squeeze(all_PE_features(i,this_word,:))'+this_PE_features;
                sum_iterative_error(i,this_word) = sum_iterative_error(i,this_word)+sum(abs(this_PE_features));
            end
            prior_features_updated = prior_features_updated+(this_PE_features*prior_update_weight);
            prior_features_updated = max(prior_features_updated,0); % Constrain priors between 0 and 1
            prior_features_updated = min(prior_features_updated,1);
        end
        posterior_features(i,this_word,:) = prior_features_updated;
        all_correlations = corr(squeeze(posterior_features(i,this_word,:)),feature_mat_words');
        [this_correlation(i,this_word),heard_word(i,this_word)] = max(all_correlations);
        heard_correct_word(i,this_word) = double(heard_word(i,this_word)==this_word);
        amount_above_average(i,this_word) = max(all_correlations)-((sum(all_correlations)-max(all_correlations))/length(all_correlations));
    end
end

       
%     condition{i}(1)
%     condition{i}(2)
%     % set prior activation levels for phoneme layer
%     prior_phonemes{1} = pa_all{1}(i);
%     prior_phonemes{2} = pa_all{2}(i);
%     
%     % set model input (to acoustic-phonetic layer)
%     input{1} = la_all{1}(i);
%     input{2} = la_all{2}(i);
%     
%     % compute phoneme-to-feature weights
%     phoneme2feature_weights{1} = pdf('norm',x,pm_all{1}(i),ps_all{1}(i)); % category A
%     phoneme2feature_weights{2} = pdf('norm',x,pm_all{2}(i),ps_all{2}(i)); % category B
%     
%     % compute input-to-feature weights
%     input2feature_weights{1} = pdf('norm',x,lm_all{1}(i),ls_all{1}(i)); % category A
%     input2feature_weights{2} = pdf('norm',x,lm_all{2}(i),ls_all{2}(i)); % category B
%     
%     % compute activation levels for acoustic-phonetic feature layer
%     like_features{1} = input2feature_weights{1}*input{1};
%     like_features{2} = input2feature_weights{2}*input{2};
%     prior_features{1} = phoneme2feature_weights{1}*prior_phonemes{1};
%     prior_features{2} = phoneme2feature_weights{2}*prior_phonemes{2};
%     PE_features{1} = like_features{1}-prior_features{1};
%     PE_features{2} = like_features{2}-prior_features{2};
%     
%     % update phoneme predictions with phoneme prediction errors
%     PE_phonemes{1} = (PE_features{1}+PE_features{2})*pinv(phoneme2feature_weights{1}); % sum feature prediction errors over categories and map to phoneme layer
%     PE_phonemes{2} = (PE_features{1}+PE_features{2})*pinv(phoneme2feature_weights{2});
%     prior_phonemes_updated{1} = prior_phonemes{1}+PE_phonemes{1}*update(i); %potentially here account for the lack of distraction in the case of mismatches - assume that participants are able to tell if the prime is incongruent
%     prior_phonemes_updated{2} = prior_phonemes{2}+PE_phonemes{2}*update(i);
%     prior_phonemes_updated_scaled{1} = prior_phonemes_updated{1} * 1/(prior_phonemes_updated{1}+prior_phonemes_updated{2}); % normalize to 1
%     prior_phonemes_updated_scaled{2} = prior_phonemes_updated{2} * 1/(prior_phonemes_updated{1}+prior_phonemes_updated{2}); % normalize to 1
%     PE_phonemes_scaled{1} = prior_phonemes_updated_scaled{1}-prior_phonemes{1};  % normalize to 1
%     PE_phonemes_scaled{2} = prior_phonemes_updated_scaled{2}-prior_phonemes{2}; % normalize to 1
%     
%     % store prediction error and identification of A vs B
%     error_all(i) = sum(abs(PE_features{1}+PE_features{2}));
%     identification_all(i) = log(prior_phonemes_updated{1})-log(prior_phonemes_updated{2});
%     raw_above_threshold(i) = (input{1}+(pa_all{1}(i)*(input{1}/(sigma_pred^2))))-threshold; % XXX New - AUC of sensory detail + weighting*AUC of sensory detail*Precision of prior (1/prior variance) - perceptual threshold
%     amount_above_threshold(i) = max(raw_above_threshold(i),0); % XXX If below perceptual threshold then zero
%                     
%     %plot activation levels - not necessary for optimisation and significantly slows it down
% %     figure(100+i);
% %     subplot(3,3,2); bar(x,like_features{1},'b'); hold on; bar(x,like_features{2},'r'); ylim([-0.2 0.2]); xlim([xdisplay(1) xdisplay(2)]); title('Sensory input','FontSize',12); set(gca,'FontSize',12);
% %     subplot(3,3,5); bar(x,prior_features{1},'b'); hold on; bar(x,prior_features{2},'r'); ylim([-0.2 0.2]); xlim([xdisplay(1) xdisplay(2)]);  title('Prediction','FontSize',12); set(gca,'FontSize',12);
% %     subplot(3,3,8); bar(x,PE_features{1},'b'); hold on; bar(x,PE_features{2},'r'); ylim([-0.2 0.2]); xlim([xdisplay(1) xdisplay(2)]);  title('Prediction error','FontSize',12); set(gca,'FontSize',12);
% %     subplot(3,3,1); bar(1,prior_phonemes{1},'b'); hold on; bar(2,prior_phonemes{2},'r'); ylim([0 1]); title('Prediction','FontSize',12); set(gca,'xtick',[1:2]); set(gca,'xticklabel',{'A' 'B'}); set(gca,'FontSize',12);
% %     subplot(3,3,4); bar(1,PE_phonemes_scaled{1},'b'); hold on; bar(2,PE_phonemes_scaled{2},'r'); ylim([-1 1]); title('Prediction error','FontSize',12); set(gca,'xtick',[1:2]); set(gca,'xticklabel',{'A' 'B'}); set(gca,'FontSize',12);
% %     subplot(3,3,7); bar(1,prior_phonemes_updated_scaled{1},'b'); hold on; bar(2,prior_phonemes_updated_scaled{2},'r'); ylim([0 1]); title('Prediction (updated)','FontSize',12); set(gca,'xtick',[1:2]); set(gca,'xticklabel',{'A' 'B'}); set(gca,'FontSize',12);
% %     subplot(3,3,3); bar(1,input{1},'b'); hold on; bar(2,input{2},'r'); ylim([0 1]); title('Model input','FontSize',12); set(gca,'xtick',[1:2]); set(gca,'xticklabel',{'A' 'B'}); set(gca,'FontSize',12);
% %     subplot(3,3,6); bar(1,(input{1}+(pa_all{1}(i)*(input{1}/(sigma_pred^2)))),'b'); hold on; ylim([0 1]); plot([0 2],[threshold threshold],'r'); title('Posterior','FontSize',12); set(gca,'xtick',[]); set(gca,'xticklabel',{}); set(gca,'FontSize',12);
% %     subplot(3,3,9); bar(1,normalised_model_meansarray(i),'b'); hold on; ylim([0 4]); plot([0 2],threshold,'r'); title('Clarity Rating','FontSize',12); set(gca,'xtick',[]); set(gca,'xticklabel',{''}); set(gca,'FontSize',12);
% 
%     % save figure of activation levels for current condition
%     %print(['Activations_condition' num2str(i) '_Sohoglu2012'],'-depsc');
%         
% end
% 
% 
% 
% %% generate sensory input and noise
% % loop through different clarity levels:
% countClarity = 1;
% IterationCounter = [];
% for clarity_level = [lowClarity highClarity]
%     
%     %% generation of normalized sensory input as probabilites (sum = 1, [0;1]
%     % loop through words to transform each word into prabability for its feature groups
%     for w = 1:nPhonemes
%         sensory_inputProb(w,:) = normalizeFeatureGroup2Prob(feature_mat_words(w,:),clarity_level);
%     end
%     
%     %% START MODEL
%     % dimensions = n(words) x n(clarity levels) x n(features)/n(words)
%     nWords = nPhonemes;
%     for v = 1:nWords
%         
%         % measurement noise for words, normally distributed
%         noiseWord_match_m1   = randn(size(prior_match_word))*2;
%         noiseWord_mismatch_m1 = randn(size(prior_match_word))*2;
%         noisefeature_match_m1   = randn(size(feature_mat_words))*2;
%         noisefeature_mismatch_m1   = randn(size(feature_mat_words))*2;
%         
%         
%         %% mismatch
%         % Run Predictive Coding loop or Sharpening loop
%         [...
%             IterationCounter(1, countClarity, v), ...
%             mismatch_feature_Accumulated(v, countClarity, :), ...
%             posterior_mismatch_word_Iterative(v, countClarity), residual_noise_mismatch(v, countClarity)] = ...
%             loopFunction(sensory_inputProb(v,:), prior_mismatch_word(v,:), ...
%             feature_mat_words, prior_update_weight, STOPcriterion); %
%         
% %         % Prepare behavioural decoding
% %         % 1. transform to probability
% %         prior_mismatch_word_Iterative = softmax(prior_mismatch_word_Iterative, temperature);
% %         % 2. store posterior for all conditions/repetitions
% %         posterior_mismatch(v, countClarity, :) = prior_mismatch_word_Iterative'; %#ok<AGROW> % store final perceptual representation
% %         
%         %% match
%         % Run Predictive Coding loop or Sharpening loop
%         [...
%             IterationCounter(2, countClarity, v), ...
%             match_feature_Accumulated(v, countClarity, :), ...
%             posterior_match_word_Iterative(v, countClarity), residual_noise_match(v, countClarity)] = ...
%             loopFunction(sensory_inputProb(v,:), prior_match_word(v,:), ...
%             feature_mat_words, prior_update_weight, STOPcriterion); %#ok<AGROW>
%         
%         % add measurement noise to "PC error signal"/"sharp signal"
%         % for RSA analysis:
%         match_feature_Accumulated(v, countClarity, :) = ...
%             squeeze(match_feature_Accumulated(v, countClarity, :)) + noisefeature_match_m1(v,:)';
%         
%         mismatch_feature_Accumulated(v, countClarity, :) = ...
%             squeeze(mismatch_feature_Accumulated(v, countClarity, :)) + noisefeature_mismatch_m1(v,:)';
%         
% %         % Prepare behavioural decoding
% %         % 1. transform to probability
% %         posterior_match_word_Iterative = softmax(posterior_match_word_Iterative, temperature);
% %         % 2. store posterior for all conditions/repetitions
% %         posterior_match(v, countClarity, :) = posterior_match_word_Iterative'; %#ok<AGROW> % store final perceptual representation
%     end
%     % increase clarity counter
%     countClarity  = countClarity +1;
% end
% 
% 
% %% Compute univariate results
% % 5 dimensions:
% % 1. 2 conditions: mismatch/Match,
% % 2. 2 x Clarity levels,
% % 3. 24 words (v)
% % 4. 21 subjects (rep1)
% % 5. 6 runs (rep2)
% % avegrage over runs, subjects and words
% % to get number of iterations per noise level and condition
% IterationCounterMean = mean(IterationCounter, 3);
% 
% %% Plot univariate results
% %     % 1. plot univariate results based on iterations
% %     figure;bar([IterationCounterMean(1,1), IterationCounterMean(2,1), ...
% %         IterationCounterMean(1,2), IterationCounterMean(2,2)]);
% %     hold on; 
% %     title(['Univariate: ' figureTitle ' - number of iterations']);
% %     set(gca, 'Xtick',1:4)
% %     set(gca, 'XTickLabel',{'mismatch 4 channel', 'match 4 channel', ...
% %         'mismatch 12 channel', 'match 12 channel'})
% %     ylabel('number of iterations')
% 
% 
% % if showBarFig
% %     
% %     IterationCounterMeanPerCondition = squeeze(mean(mean(IterationCounter, 3), 5));
% %     IterationCounterse = se(squeeze(...
% %         [IterationCounterMeanPerCondition(1,1,:)...
% %         IterationCounterMeanPerCondition(1,2,:)...
% %         IterationCounterMeanPerCondition(2,1,:)...
% %         IterationCounterMeanPerCondition(2,2,:)])');
% %     
% %     % 1. plot univariate results based on iterations
% %     figure;bar([IterationCounterMean(1,1), IterationCounterMean(2,1), ...
% %         IterationCounterMean(1,2), IterationCounterMean(2,2)]);
% %     hold on; errorbar([IterationCounterMean(1,1), IterationCounterMean(2,1), ...
% %         IterationCounterMean(1,2), IterationCounterMean(2,2)], ...
% %         [IterationCounterse(1,1), IterationCounterse(1,3), ...
% %         IterationCounterse(1,2), IterationCounterse(1,4)], '.');
% %     title(['Univariate: ' figureTitle ' - number of iterations']);
% %     set(gca, 'Xtick',1:4)
% %     set(gca, 'XTickLabel',{'mismatch 4 channel', 'match 4 channel', ...
% %         'mismatch 12 channel', 'match 12 channel'})
% %     ylabel('number of iterations')
% % 
% % end
% % 
% % %% Decode behavior
% % 
% % for r1 = 1:rep1
% %     %for r2 = 1:rep2
% %     behaviour_mismatch_LowClarity(r1)  = decodeBehavior(squeeze(mean(posterior_mismatch(:,1,:,r1),5)), behaviour_noise); %#ok<AGROW>
% %     behaviour_mismatch_HighClarity(r1) = decodeBehavior(squeeze(mean(posterior_mismatch(:,2,:,r1),5)), behaviour_noise); %#ok<AGROW>
% %     behaviour_match_LowClarity(r1)    = decodeBehavior(squeeze(mean(posterior_match(:,1,:,r1),5)), behaviour_noise); %#ok<AGROW>
% %     behaviour_match_HighClarity(r1)   = decodeBehavior(squeeze(mean(posterior_match(:,2,:,r1),5)), behaviour_noise); %#ok<AGROW>
% % end
% 
% %% Plot behavioural results
% % if showBarFig
% %     % SE after Loftus and Masson (1994) = the same for all conditions
% %     se4condition = se([behaviour_mismatch_LowClarity', behaviour_match_LowClarity', ...
% %         behaviour_mismatch_HighClarity', behaviour_match_HighClarity']);
% %     figure; bar([mean(behaviour_mismatch_LowClarity), mean(behaviour_match_LowClarity),...
% %         mean(behaviour_mismatch_HighClarity),mean(behaviour_match_HighClarity)]);
% %     hold on; errorbar([mean(behaviour_mismatch_LowClarity), mean(behaviour_match_LowClarity), ...
% %         mean(behaviour_mismatch_HighClarity), mean(behaviour_match_HighClarity)], ...
% %         [se4condition(1,1), se4condition(1,2), ...
% %         se4condition(1,3), se4condition(1,4)], '.');
% %     title(['Behavioural: ' figureTitle]);
% %     set(gca, 'Xtick',1:4)
% %     set(gca, 'XTickLabel',{'mismatch 4 channel', 'match 4 channel', ...
% %         'mismatch 12 channel', 'match 12 channel'})
% %     ylim([0 1.1])
% %     ylabel('% correct responses')
% % end
% % 
% % %% simulate RSA
% % experimental_setup.rep1 = rep1;
% % experimental_setup.rep2 = rep2;
% % experimental_setup.nPhonemes = nPhonemes;
% % 
% % results.match_word   = match_feature_Accumulated;
% % results.mismatch_word = mismatch_feature_Accumulated;
% % % within condition (i.e.: match-match and mismatch-mismatch)
% % % distance matrix with condition order: 4ch-run1, 12ch-run1, 4ch-run2, 12ch-run2
% % [simulated_RSA_results, distanceMatMatch_4ch, distanceMatMatch_12ch, ...
% %     distanceMatmismatch_4ch, distanceMatmismatch_12ch] = ...
% %     PINFA_simulateRSAwithinConditionTest_spearman(experimental_setup, results, figureTitle, showBarFig, showImageFig);
% % distanceMatMatch   = [distanceMatMatch_4ch, distanceMatMatch_12ch];
% % distanceMatmismatch = [distanceMatmismatch_4ch, distanceMatmismatch_12ch];
% % 
% % %% compute goodness of fit
% % % initialize overall error to 0
% % MSS_error = 0;
% % fittedString = '';
% % 
% % % compute discrepancy between actual and simulated univariate
% % % 1. collect simulation data
% % simulated_univariate_results = [IterationCounterMean(1,1), IterationCounterMean(2,1), ...
% %     IterationCounterMean(1,2), IterationCounterMean(2,2)];
% % % 2. normalize simulation data & experimental data
% % % (normalize max to 1)
% % univariate_results = univariate_results ./ max(univariate_results);
% % simulated_univariate_results = simulated_univariate_results ./ 500;% max(simulated_univariate_results);
% % 
% % % 3. compute error (sum of squares)
% % MSS_error_univariate = sum((univariate_results - simulated_univariate_results).^2);
% % % 4. store error
% % if returnUnivariateError
% %     MSS_error = MSS_error + MSS_error_univariate;
% %     fittedString = [fittedString 'errUniv = ' num2str(MSS_error_univariate), '; '];
% % end
% % 
% % % compute discrepancy between actual and simulated behaviour
% % % 1. collect simulation data
% % simulated_behavioural_results = [mean(behaviour_mismatch_LowClarity), mean(behaviour_match_LowClarity),...
% %     mean(behaviour_mismatch_HighClarity),mean(behaviour_match_HighClarity)];
% % % 2. normalize simulation data & experimental data
% % % (normalization not necessary since both are already percentages)
% % % 3. compute error (sum of squares)
% % MSS_error_behavioural = sum((behavioural_results - simulated_behavioural_results).^2);
% % % 4. store error
% % if returnBehaviouralError
% %     MSS_error = MSS_error + MSS_error_behavioural;
% %     fittedString = [fittedString 'errBehav = ' num2str(MSS_error_behavioural), '; '];
% % end
% % 
% % % compute discrepancy between actual and simulated RSA
% % % 1. collect simulation data
% % % (already collected in simulated_RSA_results)
% % % 2. normalize simulation data & experimental data
% % % (normalize max to 1)
% % RSA_results = RSA_results ./ max(RSA_results);
% % simulated_RSA_results = simulated_RSA_results ./ max(simulated_RSA_results);
% % % 3. compute error (sum of squares)
% % MSS_error_RSA = sum((RSA_results - simulated_RSA_results).^2);
% % % 4. store error
% % if returnRSAError
% %     MSS_error = MSS_error + MSS_error_RSA;
% %     fittedString = [fittedString 'errRSA = ' num2str(MSS_error_RSA), '; '];
% % end
% % 
% % disp(['Ran ' type ' simulationModel. x = ' mat2str(parameters, 4) '; ' fittedString]);
% 
% 
% % %% make specifications
% % type         = 'PC';
% % %type         = 'sharp';
% % 
% % if strcmp(type, 'PC')
% %     loopFunction = @predictiveCodingLoop_Sohoglumodel;
% %     figureTitle  = 'Predictive Coding';
% % else
% %     loopFunction = @sharpeningLoop;
% %     figureTitle  = 'Sharpening';
% % end
% % 
% % rng('default');  % seed the random number generator so that results are replicable...
% % 
