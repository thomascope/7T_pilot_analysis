function output = search_subfunction_bayes_modelling_PINFA(tooptimise)

% An adaptation of Ed Sohoglu's Bayesian Prediction Error account to model the PNFA MEG behavioural data
% Thomas Cope Nov 2016
%
% NEW hierarchical version
% Models two layers of representation
% Layer 1: Acoustic-phonetic features
% Layer 2: Phonemes

global normalised_model_meansarray % Outputs
global optimise_this_subject report_data meansarray t % inputs
%global mu_input_match mu_input_mismatch mu_pred %Bayesian space 

update = [ones(1,9)*0.5]; % determines how much priors are updated by prediction errors

% prior category A: pm = mean [condition 1, 2, 3 etc.] ps = standard deviation, pa = area under curve
% (simulating prior knowledge by changing area parameter)
pm_all{1} = [-2.5 -2.5 -2.5 -2.5 -2.5 -2.5 -2.5 -2.5 -2.5]; 
ps_all{1} = ones(1,9)*tooptimise(1); % Standard deviation of the prior could potentially be a floating variable for optimisation
pa_all{1} = [0 0 0 0 0 0 0.5 0.5 0.5]; %XXX This has changed as otherwise difficult to explain similar responses for mismatch and neutral. Given that there are an infinite number of possible words, I think it is reasonable for mismatch and neutral priors have the same sharpening effect. Possibly in future could optimise this, to reflect 'confusion' of a mismatch, but negatively related to sigma_pred so this is problematic.

% prior category B
pm_all{2} = [2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5]; %At the moment these values are arbitrary, and don't reflect in the distribution of the posterior. This is a candidate for optimisation if we want to model a neutral effect, but not implemented for now.
ps_all{2} = ones(1,9)*tooptimise(1);
pa_all{2} = [0.5 0.5 0.5 0 0 0 0 0 0]; 

% likelihood category A: lm = mean [condition 1, 2, 3 etc.], ls = standard deviation, la = area under curve
% (simulating sensory detail by changing area parameter)
lm_all{1} = [-2.5 -2.5 -2.5 -2.5 -2.5 -2.5 -2.5 -2.5 -2.5]; 
ls_all{1} = [2 2 2 2 2 2 2 2 2]; 

%XXX Define the precision of the sensory input based on vocode report performance (25% chance, 100% perfect)
la_all{1} = (repmat(squeeze(report_data(t,optimise_this_subject,:))',[1 3])-25)./75;

% likelihood category B
lm_all{2} = [2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5]; 
ls_all{2} = [2 2 2 2 2 2 2 2 2];
la_all{2} = ones(1,9)-la_all{1};

condition{1} = 'A heard, Mismatch prior (B>A), Low sensory detail';
condition{2} = 'A heard, Mismatch prior (B>A), Medium sensory detail';
condition{3} = 'A heard, Mismatch prior (B>A), High sensory detail';
condition{4} = 'A heard, Neutral prior (A=B), Low sensory detail';
condition{5} = 'A heard, Neutral prior (A=B), Medium sensory detail';
condition{6} = 'A heard, Neutral prior (A=B), High sensory detail';
condition{7} = 'A heard, Matching prior (A>B), Low sensory detail';
condition{8} = 'A heard, Matching prior (A>B), Medium sensory detail';
condition{9} = 'A heard, Matching prior (A>B), High sensory detail';

x = [-10:1:10];
xdisplay = [-20 20];

for i=1:length(pm_all{1})
        
    % set prior activation levels for phoneme layer
    prior_phonemes{1} = pa_all{1}(i);
    prior_phonemes{2} = pa_all{2}(i);
    
    % set model input (to acoustic-phonetic layer)
    input{1} = la_all{1}(i);
    input{2} = la_all{2}(i);
    
    % compute phoneme-to-feature weights
    phoneme2feature_weights{1} = pdf('norm',x,pm_all{1}(i),ps_all{1}(i)); % category A
    phoneme2feature_weights{2} = pdf('norm',x,pm_all{2}(i),ps_all{2}(i)); % category B
    
    % compute input-to-feature weights
    input2feature_weights{1} = pdf('norm',x,lm_all{1}(i),ls_all{1}(i)); % category A
    input2feature_weights{2} = pdf('norm',x,lm_all{2}(i),ls_all{2}(i)); % category B
    
    % compute activation levels for acoustic-phonetic feature layer
    like_features{1} = input2feature_weights{1}*input{1};
    like_features{2} = input2feature_weights{2}*input{2};
    prior_features{1} = phoneme2feature_weights{1}*prior_phonemes{1};
    prior_features{2} = phoneme2feature_weights{2}*prior_phonemes{2};
    PE_features{1} = like_features{1}-prior_features{1};
    PE_features{2} = like_features{2}-prior_features{2};
    
    % update phoneme predictions with phoneme prediction errors
    PE_phonemes{1} = (PE_features{1}+PE_features{2})*pinv(phoneme2feature_weights{1}); % sum feature prediction errors over categories and map to phoneme layer
    PE_phonemes{2} = (PE_features{1}+PE_features{2})*pinv(phoneme2feature_weights{2});
    prior_phonemes_updated{1} = prior_phonemes{1}+PE_phonemes{1}*update(i); %potentially here account for the lack of distraction in the case of mismatches - assume that participants are able to tell if the prime is incongruent
    prior_phonemes_updated{2} = prior_phonemes{2}+PE_phonemes{2}*update(i);
    prior_phonemes_updated_scaled{1} = prior_phonemes_updated{1} * 1/(prior_phonemes_updated{1}+prior_phonemes_updated{2}); % normalize to 1
    prior_phonemes_updated_scaled{2} = prior_phonemes_updated{2} * 1/(prior_phonemes_updated{1}+prior_phonemes_updated{2}); % normalize to 1
    PE_phonemes_scaled{1} = prior_phonemes_updated_scaled{1}-prior_phonemes{1};  % normalize to 1
    PE_phonemes_scaled{2} = prior_phonemes_updated_scaled{2}-prior_phonemes{2}; % normalize to 1
    
    % store prediction error and identification of A vs B
    error_all(i) = sum(abs(PE_features{1}+PE_features{2}));
    identification_all(i) = log(prior_phonemes_updated{1})-log(prior_phonemes_updated{2});
    raw_above_threshold(i) = (input{1}+(pa_all{1}(i)*(input{1}/(tooptimise(1)^2))))-tooptimise(2); % XXX New - AUC of sensory detail + weighting*AUC of sensory detail*Precision of prior (1/prior variance) - perceptual threshold
    amount_above_threshold(i) = max(raw_above_threshold(i),0); % XXX If below perceptual threshold then zero
                    
    %plot activation levels - not necessary for optimisation and significantly slows it down
%     figure(100+i);
%     subplot(3,3,2); bar(x,like_features{1},'b'); hold on; bar(x,like_features{2},'r'); ylim([-0.2 0.2]); xlim([xdisplay(1) xdisplay(2)]); title('Sensory input','FontSize',12); set(gca,'FontSize',12);
%     subplot(3,3,5); bar(x,prior_features{1},'b'); hold on; bar(x,prior_features{2},'r'); ylim([-0.2 0.2]); xlim([xdisplay(1) xdisplay(2)]);  title('Prediction','FontSize',12); set(gca,'FontSize',12);
%     subplot(3,3,8); bar(x,PE_features{1},'b'); hold on; bar(x,PE_features{2},'r'); ylim([-0.2 0.2]); xlim([xdisplay(1) xdisplay(2)]);  title('Prediction error','FontSize',12); set(gca,'FontSize',12);
%     subplot(3,3,1); bar(1,prior_phonemes{1},'b'); hold on; bar(2,prior_phonemes{2},'r'); ylim([0 1]); title('Prediction','FontSize',12); set(gca,'xtick',[1:2]); set(gca,'xticklabel',{'A' 'B'}); set(gca,'FontSize',12);
%     subplot(3,3,4); bar(1,PE_phonemes_scaled{1},'b'); hold on; bar(2,PE_phonemes_scaled{2},'r'); ylim([-.2 .2]); title('Prediction error','FontSize',12); set(gca,'xtick',[1:2]); set(gca,'xticklabel',{'A' 'B'}); set(gca,'FontSize',12);
%     subplot(3,3,7); bar(1,prior_phonemes_updated_scaled{1},'b'); hold on; bar(2,prior_phonemes_updated_scaled{2},'r'); ylim([0 1]); title('Prediction (updated)','FontSize',12); set(gca,'xtick',[1:2]); set(gca,'xticklabel',{'A' 'B'}); set(gca,'FontSize',12);
%     subplot(3,3,3); bar(1,input{1},'b'); hold on; bar(2,input{2},'r'); ylim([0 1]); title('Model input','FontSize',12); set(gca,'xtick',[1:2]); set(gca,'xticklabel',{'A' 'B'}); set(gca,'FontSize',12);
%     subplot(3,3,6); bar(1,(input{1}+(pa_all{1}(i)*(input{1}/(tooptimise(1)^2)))),'b'); hold on; ylim([0 1]); plot([0 2],[tooptimise(2) tooptimise(2)],'r'); title('Posterior','FontSize',12); set(gca,'xtick',[]); set(gca,'xticklabel',{}); set(gca,'FontSize',12);
%     subplot(3,3,9); bar(1,normalised_model_meansarray(i),'b'); hold on; ylim([0 4]); plot([0 2],tooptimise(2),'r'); title('Clarity Rating','FontSize',12); set(gca,'xtick',[]); set(gca,'xticklabel',{''}); set(gca,'FontSize',12);

    % save figure of activation levels for current condition
    %print(['Activations_condition' num2str(i) '_Sohoglu2012'],'-depsc');
        
end

% figure(200);
% subplot(3,1,1); bar(error_all,'k'); ylim([0 1.5]); title('Prediction errors');
% subplot(3,1,2); bar(identification_all,'k'); ylim([-1 2.5]); title('Identification: A vs B');
% subplot(3,1,3); bar(amount_above_threshold,'k'); ylim([-1 2.5]); title('Modelled Clarity Rating');

% get max and min of meansarray
model_meansarray = amount_above_threshold;
maxVec = max(model_meansarray(:));
minVec = tooptimise(2); %This seems initially counter-intuitive, including the threshold as the minimum value on the rating scale even though it has already been taken off the amount above threshold, but it allows non-linearity in the clarity response that is desirable, and it has face-value validity to have the bottom of the scale as perceptual threshold.
maxVec_R = max(meansarray(:));
minVec_R = min(meansarray(:));
normalised_model_meansarray = ((model_meansarray-minVec)./(maxVec-minVec)).*(maxVec_R-minVec_R)+minVec_R;
meansarray_R = ((meansarray-minVec_R)./(maxVec_R-minVec_R)).*(maxVec_R-minVec_R)+minVec_R;
    
differences = [abs(normalised_model_meansarray(1:3) - meansarray_R(1:3)),abs(normalised_model_meansarray(7:9) - meansarray_R(7:9))];
    
%output = sum(differences(:)); %Unclear whether it makes more face-value sense to minimise residual error or squared residual error, but I have tried it both ways and the statistical outputs don't differ.
output = sum(differences.^2); 