function [ feature_mat_words ] = PINFA_word2phoneme( phonemelist )
% This function transforms words from DISC transcription into feature
% vectors (words have a consonant - vowel - consonant structure)
% From Blank and Davis 2016:
% Representing the segments of the 24 words in our stimulus set required 13
% consonantal features and 11 vowel features, concatenated into a set of 37
% binary features for the CVC syllables used in the experiment. The 13
% consonantal features were divided into four groups: (1) place of
% articulation (six features: bilabial, labiodental, dental, alveolar,
% palato-alveolar, velar), (2) manner of articulation (three features:
% stop, sibilant, non-sibilant), (3) nasality (three features: nasal,
% oral), and (4) voicing (two features: voiceless, voiced). The 11 vowel
% features were divided into four groups: (1) height (five features: high,
% mid-high, mid, mid-low, low), (2) backness (two features: front, back),
% (3) rounding (two features: rounded, unrounded), and (4) length/diphthong
% (two features: long, short). Based on these position-specific features,
% we constructed a feature-to-word transformation matrix that included
% positive binary values in each row to indicate which phonetic features
% were relevant for each word (see S2 Fig). Each row contained 12 active
% features (four features for each consonant and vowel). This matrix served
% as a set of connection weights to link phonetic features to words in both
% models and thereby encoded long-term knowledge of the form of each spoken
% word.

vowelsDISC     = {'I', 'i', '#', 'D', 'u', '1', '2', '5', 'E'};
consonantsDISC = {'T', 'N', 's', 't', 'd', 'p', 'k', 'b', ...
    'f', 'm', 'S', 'z', 'n', 'l', '_', 'r'};

% 11 vowel features - In Helen's code said 10, but actually 11
% 24.06.21 E Added by TEC
vowels_to_features_all10 = [...
    0	1	0	0	0	1	0	0	1	0  1; ...
    1	0	0	0	0	1	0	0	1	0  1; ...
    0	0	0	0	1	0	1	0	1	0  1; ...
    0	0	0	1	0	0	1	1	0	0  1; ...
    1	0	0	0	0	0	1	1	0	0  1; ...
    0	1	0	0	0	1	0	0	1	1  0; ...
    1	0	0	0	0	1	0	0	1	1  0; ...
    0	0	1	0	0	0	1	1	0	1  0; ...
    0	0	1	0	0   1   0   0   1   0  1]; 

% This matrix was corrected on 05.03.2015 - f = non-sibilant
% 24.06.21 l _ and r Added by TEC
% 13 consonant features
consonants_to_features = [...
    0	0	1	0	0	0	0	0	1	0	1	1	0; ...
    0	0	0	0	0	1	1	0	0	1	0	0	1; ...
    0	0	0	1	0	0	0	1	0	0	1	1	0; ...
    0	0	0	1	0	0	1	0	0	0	1	1	0; ...
    0	0	0	1	0	0	1	0	0	0	1	0	1; ...
    1	0	0	0	0	0	1	0	0	0	1	1	0; ...
    0	0	0	0	0	1	1	0	0	0	1	1	0; ...
    1	0	0	0	0	0	1	0	0	0	1	0	1; ...
    0	1	0	0	0	0	0	0	1	0	1	1	0; ...
    1	0	0	0	0	0	1	0	0	1	0	0	1; ...
    0	0	0	0	1	0	0	1	0	0	1	1	0; ...
    0	0	0	1	0	0	0	1	0	0	1	0	1; ...
    0	0	0	1	0	0	1	0	0	1	0	0	1; ...
    0   0   0   1   0   0   0   0   1   0   1   0   1; ...
    0   0   0   0   1   0   0   1   0   0   1   0   1; ...
    0   0   0   0   1   0   0   0   1   0   1   0   1]; 

% initialize final output matrix
feature_mat_words = zeros(length(phonemelist), 37);
% loop throught phoneme list
for i = 1:numel(phonemelist)
    phoneme = phonemelist{i};
    curr_consonant1 = phoneme(1);
    curr_vowel      = phoneme(2);
    curr_consonant2 = phoneme(3);
    
    feature_consonant1 = ~cellfun('isempty',(regexp(consonantsDISC, curr_consonant1)));
    feature_vowel      = ~cellfun('isempty',(regexp(vowelsDISC, curr_vowel)));
    feature_consonant2 = ~cellfun('isempty',(regexp(consonantsDISC, curr_consonant2)));
    feature_mat_words(i,:) = [consonants_to_features(feature_consonant1, :), ...
        vowels_to_features_all10(feature_vowel, :), ...
        consonants_to_features(feature_consonant2, :)];
end

end