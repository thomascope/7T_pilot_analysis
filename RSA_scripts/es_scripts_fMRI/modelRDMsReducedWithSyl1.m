function Models = modelRDMs()

fldnm = 'Syl2StrongM';
Models.(fldnm) = ones(16,16)/120; % Syl2 Model
Models.(fldnm) = kron([1 0; 0 0],Models.(fldnm)); % Extend across Strong and Weak
Models.(fldnm) = kron([1 0; 0 0],Models.(fldnm)); % Extend across Matching and Mismatching
Models.(fldnm)(65:80,:) = 0; % NoiseInitial
Models.(fldnm)(:,65:80) = 0;
Models.(fldnm)(81:144,:) = 0; % Syl1
Models.(fldnm)(:,81:144) = 0;
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = 0; % Make diagonal 0
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'Syl2WeakM';
Models.(fldnm) = ones(16,16)/120; % Syl2 Model
Models.(fldnm) = kron([0 0; 0 1],Models.(fldnm)); % Extend across Strong and Weak
Models.(fldnm) = kron([1 0; 0 0],Models.(fldnm)); % Extend across Matching and Mismatching
Models.(fldnm)(65:80,:) = 0; % NoiseInitial
Models.(fldnm)(:,65:80) = 0;
Models.(fldnm)(81:144,:) = 0; % Syl1
Models.(fldnm)(:,81:144) = 0;
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = 0; % Make diagonal 0
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'Syl2StrongMM';
Models.(fldnm) = ones(16,16)/120; % Syl2 Model
Models.(fldnm) = kron([1 0; 0 0],Models.(fldnm)); % Extend across Strong and Weak
Models.(fldnm) = kron([0 0; 0 1],Models.(fldnm)); % Extend across Matching and Mismatching
Models.(fldnm)(65:80,:) = 0; % NoiseInitial
Models.(fldnm)(:,65:80) = 0;
Models.(fldnm)(81:144,:) = 0; % Syl1
Models.(fldnm)(:,81:144) = 0;
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = 0; % Make diagonal 0
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'Syl2WeakMM';
Models.(fldnm) = ones(16,16)/120; % Syl2 Model
Models.(fldnm) = kron([0 0; 0 1],Models.(fldnm)); % Extend across Strong and Weak
Models.(fldnm) = kron([0 0; 0 1],Models.(fldnm)); % Extend across Matching and Mismatching
Models.(fldnm)(65:80,:) = 0; % NoiseInitial
Models.(fldnm)(:,65:80) = 0;
Models.(fldnm)(81:144,:) = 0; % Syl1
Models.(fldnm)(:,81:144) = 0;
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = 0; % Make diagonal 0
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'Syl2NoiseSpeech';
Models.(fldnm)(65:80,65:80) = ones(16,16)/120;
Models.(fldnm)(1:64,:) = 0; % % Clear Strong/Weak / Matching/Mismatching
Models.(fldnm)(:,1:64) = 0;
Models.(fldnm)(81:144,:) = 0; % Syl1
Models.(fldnm)(:,81:144) = 0;
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = 0; % Make diagonal 0
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'Syl1Strong';
Models.(fldnm)(81:112,81:112) = ones(32,32)/496;
Models.(fldnm)(113:144,:) = 0; % Syl1 Weak
Models.(fldnm)(:,113:144) = 0;
Models.(fldnm)(1:64,:) = 0; % % Clear Strong/Weak / Matching/Mismatching
Models.(fldnm)(:,1:64) = 0;
Models.(fldnm)(65:80,:) = 0; % NoiseInitial
Models.(fldnm)(:,65:80) = 0;
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = 0; % Make diagonal 0
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'Syl1Weak';
Models.(fldnm)(113:144,113:144) = ones(32,32)/496;
Models.(fldnm)(81:112,:) = 0; % Syl1 Strong
Models.(fldnm)(:,81:112) = 0;
Models.(fldnm)(1:64,:) = 0; % % Clear Strong/Weak / Matching/Mismatching
Models.(fldnm)(:,1:64) = 0;
Models.(fldnm)(65:80,:) = 0; % NoiseInitial
Models.(fldnm)(:,65:80) = 0;
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = 0; % Make diagonal 0
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'Syl2Graded';
load('RDM_phon_syl2.mat');
Models.(fldnm)(1:64,1:64) = RDM(1:2:128,1:2:128);
Models.(fldnm)(65:80,65:80) = Models.(fldnm)(1:16,1:16); % NoiseInitial
Models.(fldnm)(65:80,1:64) = Models.(fldnm)(1:16,1:64);
Models.(fldnm)(1:64,65:80) = Models.(fldnm)(1:64,1:16);
Models.(fldnm)(81:144,:) = NaN; % Syl1
Models.(fldnm)(:,81:144) = NaN;
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'Syl2GradedStrongM';
load('RDM_phon_syl2.mat');
Models.(fldnm)(1:16,1:16) = RDM(1:2:32,1:2:32);
Models.(fldnm) = kron([1 NaN; NaN NaN],Models.(fldnm)); % Extend across Strong and Weak
Models.(fldnm) = kron([1 NaN; NaN NaN],Models.(fldnm)); % Extend across Matching and Mismatching
Models.(fldnm)(65:80,:) = NaN; % NoiseInitial
Models.(fldnm)(:,65:80) = NaN;
Models.(fldnm)(81:144,:) = NaN; % Syl1
Models.(fldnm)(:,81:144) = NaN;
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'Syl2GradedWeakM';
load('RDM_phon_syl2.mat');
Models.(fldnm)(1:16,1:16) = RDM(1:2:32,1:2:32);
Models.(fldnm) = kron([NaN NaN; NaN 1],Models.(fldnm)); % Extend across Strong and Weak
Models.(fldnm) = kron([1 NaN; NaN NaN],Models.(fldnm)); % Extend across Matching and Mismatching
Models.(fldnm)(65:80,:) = NaN; % NoiseInitial
Models.(fldnm)(:,65:80) = NaN;
Models.(fldnm)(81:144,:) = NaN; % Syl1
Models.(fldnm)(:,81:144) = NaN;
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'Syl2GradedStrongMM';
load('RDM_phon_syl2.mat');
Models.(fldnm)(1:16,1:16) = RDM(1:2:32,1:2:32);
Models.(fldnm) = kron([1 NaN; NaN NaN],Models.(fldnm)); % Extend across Strong and Weak
Models.(fldnm) = kron([NaN NaN; NaN 1],Models.(fldnm)); % Extend across Matching and Mismatching
Models.(fldnm)(65:80,:) = NaN; % NoiseInitial
Models.(fldnm)(:,65:80) = NaN;
Models.(fldnm)(81:144,:) = NaN; % Syl1
Models.(fldnm)(:,81:144) = NaN;
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'Syl2GradedWeakMM';
load('RDM_phon_syl2.mat');
Models.(fldnm)(1:16,1:16) = RDM(1:2:32,1:2:32);
Models.(fldnm) = kron([NaN NaN; NaN 1],Models.(fldnm)); % Extend across Strong and Weak
Models.(fldnm) = kron([NaN NaN; NaN 1],Models.(fldnm)); % Extend across Matching and Mismatching
Models.(fldnm)(65:80,:) = NaN; % NoiseInitial
Models.(fldnm)(:,65:80) = NaN;
Models.(fldnm)(81:144,:) = NaN; % Syl1
Models.(fldnm)(:,81:144) = NaN;
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'Syl2GradedNoiseSpeech';
load('RDM_phon_syl2.mat');
Models.(fldnm)(65:80,65:80) = RDM(1:2:32,1:2:32);
Models.(fldnm)(1:64,:) = NaN; % % Clear Strong/Weak / Matching/Mismatching
Models.(fldnm)(:,1:64) = NaN;
Models.(fldnm)(81:144,:) = NaN; % Syl1
Models.(fldnm)(:,81:144) = NaN;
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'Syl1Graded';
load('RDM_phon_syl1.mat');
Models.(fldnm)(81:144,81:144) = RDM(1:64,1:64);
Models.(fldnm)(1:64,:) = NaN; % % Clear Strong/Weak / Matching/Mismatching
Models.(fldnm)(:,1:64) = NaN;
Models.(fldnm)(65:80,:) = NaN; % NoiseInitial
Models.(fldnm)(:,65:80) = NaN;
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'Syl1GradedStrong';
load('RDM_phon_syl1.mat');
Models.(fldnm)(81:112,81:112) = RDM(1:32,1:32);
Models.(fldnm)(113:144,:) = NaN; % Syl1 Weak
Models.(fldnm)(:,113:144) = NaN;
Models.(fldnm)(1:64,:) = NaN; % % Clear Strong/Weak / Matching/Mismatching
Models.(fldnm)(:,1:64) = NaN;
Models.(fldnm)(65:80,:) = NaN; % NoiseInitial
Models.(fldnm)(:,65:80) = NaN;
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'Syl1GradedWeak';
load('RDM_phon_syl1.mat');
Models.(fldnm)(113:144,113:144) = RDM(33:64,33:64);
Models.(fldnm)(81:112,:) = NaN; % Syl1 Strong
Models.(fldnm)(:,81:112) = NaN;
Models.(fldnm)(1:64,:) = NaN; % % Clear Strong/Weak / Matching/Mismatching
Models.(fldnm)(:,1:64) = NaN;
Models.(fldnm)(65:80,:) = NaN; % NoiseInitial
Models.(fldnm)(:,65:80) = NaN;
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);
