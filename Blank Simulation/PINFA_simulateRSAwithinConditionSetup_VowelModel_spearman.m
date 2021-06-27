function [dist_4ch, dist_12ch, rs_4ch, ps_4ch, rs_12ch, ps_12ch] = ...
    PINFA_simulateRSAwithinConditionSetup_VowelModel_spearman(experimental_setup, wordMatrix, conditionName, modelName, showFig)
% This function calculates the similarity and dissimilarity values
% for the 4 different combinations WITHIN conditions
% Neutral4, Match 4, Neutral 12, Match 12
% and the distance matrix with conditions in this order AVERAGED ACROSS RUNS:
% Neutral4, Match 4, Neutral 12, Match 12

% Input : experimental_setup, wordMatrix, conditionName, modelName, showFig
% Output: similarity and dissimiliarity matrices for the 4 conditions
% and the disctance matrix with the dimensions subjects x n items x n items

rep1 = experimental_setup.rep1;
rep2 = experimental_setup.rep2;
nPhonemes = experimental_setup.nPhonemes;

% size is 24 words x 2 clarity x 24 wordEvidences x n=rep repetitions
% reshape data so that new dimensions are
% 48 words (24 channel 4 + 24 channel 12) x 24 wordEvidences x n=rep repetitions
wordMAT_4ch = reshape(wordMatrix(:,1,:,:), nPhonemes, nPhonemes, rep1, rep2);
wordMAT_12ch = reshape(wordMatrix(:,2,:,:), nPhonemes, nPhonemes, rep1, rep2);

%VowelTripleNanDiag = [NaN 0 0; 0 NaN 0 ; 0 0 NaN ];
shared_segments = zeros(16,16);
shared_segments(1:17:end) = 1;
shared_segments(2:68:end) = 2/3;
shared_segments(3:68:end) = 2/3;
shared_segments(4:68:end) = 1/3;
shared_segments(17:68:end) = 2/3;
shared_segments(19:68:end) = 1/3;
shared_segments(20:68:end) = 2/3;
shared_segments(33:68:end) = 2/3;
shared_segments(34:68:end) = 1/3;
shared_segments(36:68:end) = 2/3;
shared_segments(49:68:end) = 1/3;
shared_segments(50:68:end) = 2/3;
shared_segments(51:68:end) = 2/3;

shared_segments(1,16) = 1/3;
shared_segments(1,14) = 1/3;
shared_segments(16,1) = 1/3;
shared_segments(14,1) = 1/3;
shared_segments(3,16) = 1/3;
shared_segments(3,14) = 1/3;
shared_segments(16,3) = 1/3;
shared_segments(14,3) = 1/3;

shared_segments(5,9) = 1/3;
shared_segments(7,9) = 1/3;
shared_segments(9,5) = 1/3;
shared_segments(9,7) = 1/3;
shared_segments(5,11) = 1/3;
shared_segments(7,11) = 1/3;
shared_segments(11,5) = 1/3;
shared_segments(11,7) = 1/3;

shared_segments(6,10) = 1/3;
shared_segments(8,10) = 1/3;
shared_segments(10,6) = 1/3;
shared_segments(10,8) = 1/3;
shared_segments(6,12) = 1/3;
shared_segments(8,12) = 1/3;
shared_segments(12,6) = 1/3;
shared_segments(12,8) = 1/3;

shared_segments(15,3) = 1/3;
shared_segments(15,4) = 1/3;
shared_segments(16,4) = 1/3;
shared_segments(16,3) = 1/3;
shared_segments(3,16) = 1/3;
shared_segments(4,16) = 1/3;
shared_segments(4,15) = 1/3;
shared_segments(3,15) = 1/3;

shared_segments = 1-shared_segments;

% compute distance between word representations
for s = 1:rep1
    % compute similarity for each subject (s)
    % order of condition in the distance matrix:
    % 4ch all runs, 12 ch all runs
    
    % average across runs
    dist_4ch(s,:,:) = squareform(pdist(mean(wordMAT_4ch(:,:,s,:),4), 'correlation')); % #ok<AGROW>
    dist_12ch(s,:,:) = squareform(pdist(mean(wordMAT_12ch(:,:,s,:),4), 'correlation')); % #ok<AGROW>
    
    searchlightRDMs_4ch = squareform(squeeze(dist_4ch(s,:,:))');
    % spearman correlation
    [rs_4ch(s), ps_4ch(s)] = corr(searchlightRDMs_4ch', squareform(shared_segments)',...
        'type', 'Spearman', 'rows', 'pairwise');
    % [rs, ps] = corr(searchlightRDM', modelRDMs_ltv', 'type', 'Spearman', 'rows', 'pairwise');
    % Fisher transformation
    rs_4ch(s) = fisherTransform(rs_4ch(s));
    
    searchlightRDMs_12ch = squareform(squeeze(dist_12ch(s,:,:))');
    % spearman correlation
    [rs_12ch(s), ps_12ch(s)] = corr(searchlightRDMs_12ch', squareform(shared_segments)',...
        'type', 'Spearman', 'rows', 'pairwise');
    % [rs, ps] = corr(searchlightRDM', modelRDMs_ltv', 'type', 'Spearman', 'rows', 'pairwise');
    % Fisher transformation
    rs_12ch(s) = fisherTransform(rs_12ch(s));
end

%% plot average similarity across subjects
if showFig
    maxFig = max(max(max([dist_4ch; dist_4ch])));
    dist_meanSubj = squeeze(mean(dist_4ch,1));
    figure;
    imagesc(dist_meanSubj);
    caxis([0,maxFig]);
    colorbar;
    title([conditionName ' 4ch ' modelName]);
    
    maxFig = max(max(max([dist_12ch; dist_12ch])));
    dist_meanSubj = squeeze(mean(dist_12ch,1));
    figure;
    imagesc(dist_meanSubj);
    caxis([0,maxFig]);
    colorbar;
    title([conditionName ' 12ch ' modelName]);
end