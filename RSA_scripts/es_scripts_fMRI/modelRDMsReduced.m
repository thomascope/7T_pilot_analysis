function Models = modelRDMs()

fldnm = 'Syl2';
Models.(fldnm) = [0,NaN,1,1,1,1,1,1,1,1,1,1,1,1,1,1;NaN,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,0,NaN,1,1,1,1,1,1,1,1,1,1,1,1;1,1,NaN,0,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,0,NaN,1,1,1,1,1,1,1,1,1,1;1,1,1,1,NaN,0,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,0,NaN,1,1,1,1,1,1,1,1;1,1,1,1,1,1,NaN,0,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,0,NaN,1,1,1,1,1,1;1,1,1,1,1,1,1,1,NaN,0,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,0,NaN,1,1,1,1;1,1,1,1,1,1,1,1,1,1,NaN,0,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,0,NaN,1,1;1,1,1,1,1,1,1,1,1,1,1,1,NaN,0,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,NaN;1,1,1,1,1,1,1,1,1,1,1,1,1,1,NaN,0]; % Syl2 Model
Models.(fldnm) = kron([1 1; 1 1],Models.(fldnm)); % Extend across Strong and Weak
Models.(fldnm) = kron([1 1; 1 1],Models.(fldnm)); % Extend across Matching and Mismatching
Models.(fldnm)(65:80,:) = NaN; % NoiseFinal Strong
Models.(fldnm)(:,65:80) = NaN;
Models.(fldnm)(81:96,:) = NaN; % NoiseFinal Weak
Models.(fldnm)(:,81:96) = NaN;
Models.(fldnm)(97:112,97:112) = Models.(fldnm)(1:16,1:16); % NoiseInitial
Models.(fldnm)(97:112,1:96) = Models.(fldnm)(1:16,1:96);
Models.(fldnm)(1:96,97:112) = Models.(fldnm)(1:96,1:16);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = NaN;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'Syl2StrongM';
Models.(fldnm) = [0,NaN,1,1,1,1,1,1,1,1,1,1,1,1,1,1;NaN,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,0,NaN,1,1,1,1,1,1,1,1,1,1,1,1;1,1,NaN,0,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,0,NaN,1,1,1,1,1,1,1,1,1,1;1,1,1,1,NaN,0,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,0,NaN,1,1,1,1,1,1,1,1;1,1,1,1,1,1,NaN,0,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,0,NaN,1,1,1,1,1,1;1,1,1,1,1,1,1,1,NaN,0,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,0,NaN,1,1,1,1;1,1,1,1,1,1,1,1,1,1,NaN,0,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,0,NaN,1,1;1,1,1,1,1,1,1,1,1,1,1,1,NaN,0,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,NaN;1,1,1,1,1,1,1,1,1,1,1,1,1,1,NaN,0]; % Syl2 Model
Models.(fldnm) = kron([1 NaN; NaN NaN],Models.(fldnm)); % Extend across Strong and Weak
Models.(fldnm) = kron([1 NaN; NaN NaN],Models.(fldnm)); % Extend across Matching and Mismatching
Models.(fldnm)(65:80,:) = NaN; % NoiseFinal Strong
Models.(fldnm)(:,65:80) = NaN;
Models.(fldnm)(81:96,:) = NaN; % NoiseFinal Weak
Models.(fldnm)(:,81:96) = NaN;
Models.(fldnm)(97:112,:) = NaN; % NoiseInitial
Models.(fldnm)(:,97:112) = NaN;
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = NaN;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'Syl2WeakM';
Models.(fldnm) = [0,NaN,1,1,1,1,1,1,1,1,1,1,1,1,1,1;NaN,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,0,NaN,1,1,1,1,1,1,1,1,1,1,1,1;1,1,NaN,0,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,0,NaN,1,1,1,1,1,1,1,1,1,1;1,1,1,1,NaN,0,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,0,NaN,1,1,1,1,1,1,1,1;1,1,1,1,1,1,NaN,0,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,0,NaN,1,1,1,1,1,1;1,1,1,1,1,1,1,1,NaN,0,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,0,NaN,1,1,1,1;1,1,1,1,1,1,1,1,1,1,NaN,0,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,0,NaN,1,1;1,1,1,1,1,1,1,1,1,1,1,1,NaN,0,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,NaN;1,1,1,1,1,1,1,1,1,1,1,1,1,1,NaN,0]; % Syl2 Model
Models.(fldnm) = kron([NaN NaN; NaN 1],Models.(fldnm)); % Extend across Strong and Weak
Models.(fldnm) = kron([1 NaN; NaN NaN],Models.(fldnm)); % Extend across Matching and Mismatching
Models.(fldnm)(65:80,:) = NaN; % NoiseFinal Strong
Models.(fldnm)(:,65:80) = NaN;
Models.(fldnm)(81:96,:) = NaN; % NoiseFinal Weak
Models.(fldnm)(:,81:96) = NaN;
Models.(fldnm)(97:112,:) = NaN; % NoiseInitial
Models.(fldnm)(:,97:112) = NaN;
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = NaN;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'Syl2StrongMM';
Models.(fldnm) = [0,NaN,1,1,1,1,1,1,1,1,1,1,1,1,1,1;NaN,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,0,NaN,1,1,1,1,1,1,1,1,1,1,1,1;1,1,NaN,0,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,0,NaN,1,1,1,1,1,1,1,1,1,1;1,1,1,1,NaN,0,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,0,NaN,1,1,1,1,1,1,1,1;1,1,1,1,1,1,NaN,0,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,0,NaN,1,1,1,1,1,1;1,1,1,1,1,1,1,1,NaN,0,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,0,NaN,1,1,1,1;1,1,1,1,1,1,1,1,1,1,NaN,0,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,0,NaN,1,1;1,1,1,1,1,1,1,1,1,1,1,1,NaN,0,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,NaN;1,1,1,1,1,1,1,1,1,1,1,1,1,1,NaN,0]; % Syl2 Model
Models.(fldnm) = kron([1 NaN; NaN NaN],Models.(fldnm)); % Extend across Strong and Weak
Models.(fldnm) = kron([NaN NaN; NaN 1],Models.(fldnm)); % Extend across Matching and Mismatching
Models.(fldnm)(65:80,:) = NaN; % NoiseFinal Strong
Models.(fldnm)(:,65:80) = NaN;
Models.(fldnm)(81:96,:) = NaN; % NoiseFinal Weak
Models.(fldnm)(:,81:96) = NaN;
Models.(fldnm)(97:112,:) = NaN; % NoiseInitial
Models.(fldnm)(:,97:112) = NaN;
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = NaN;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'Syl2WeakMM';
Models.(fldnm) = [0,NaN,1,1,1,1,1,1,1,1,1,1,1,1,1,1;NaN,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,0,NaN,1,1,1,1,1,1,1,1,1,1,1,1;1,1,NaN,0,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,0,NaN,1,1,1,1,1,1,1,1,1,1;1,1,1,1,NaN,0,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,0,NaN,1,1,1,1,1,1,1,1;1,1,1,1,1,1,NaN,0,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,0,NaN,1,1,1,1,1,1;1,1,1,1,1,1,1,1,NaN,0,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,0,NaN,1,1,1,1;1,1,1,1,1,1,1,1,1,1,NaN,0,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,0,NaN,1,1;1,1,1,1,1,1,1,1,1,1,1,1,NaN,0,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,NaN;1,1,1,1,1,1,1,1,1,1,1,1,1,1,NaN,0]; % Syl2 Model
Models.(fldnm) = kron([NaN NaN; NaN 1],Models.(fldnm)); % Extend across Strong and Weak
Models.(fldnm) = kron([NaN NaN; NaN 1],Models.(fldnm)); % Extend across Matching and Mismatching
Models.(fldnm)(65:80,:) = NaN; % NoiseFinal Strong
Models.(fldnm)(:,65:80) = NaN;
Models.(fldnm)(81:96,:) = NaN; % NoiseFinal Weak
Models.(fldnm)(:,81:96) = NaN;
Models.(fldnm)(97:112,:) = NaN; % NoiseInitial
Models.(fldnm)(:,97:112) = NaN;
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = NaN;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'Syl2StrongNoise';
Models.(fldnm)(65:80,65:80) = [0,NaN,1,1,1,1,1,1,1,1,1,1,1,1,1,1;NaN,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,0,NaN,1,1,1,1,1,1,1,1,1,1,1,1;1,1,NaN,0,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,0,NaN,1,1,1,1,1,1,1,1,1,1;1,1,1,1,NaN,0,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,0,NaN,1,1,1,1,1,1,1,1;1,1,1,1,1,1,NaN,0,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,0,NaN,1,1,1,1,1,1;1,1,1,1,1,1,1,1,NaN,0,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,0,NaN,1,1,1,1;1,1,1,1,1,1,1,1,1,1,NaN,0,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,0,NaN,1,1;1,1,1,1,1,1,1,1,1,1,1,1,NaN,0,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,NaN;1,1,1,1,1,1,1,1,1,1,1,1,1,1,NaN,0]; % Syl2 Model
Models.(fldnm)(81:96,:) = NaN; % NoiseFinal Weak
Models.(fldnm)(:,81:96) = NaN;
Models.(fldnm)(1:64,:) = NaN; % % Clear Strong/Weak / Matching/Mismatching
Models.(fldnm)(:,1:64) = NaN;
Models.(fldnm)(97:112,:) = NaN; % NoiseInitial
Models.(fldnm)(:,97:112) = NaN;
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = NaN;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'Syl2WeakNoise';
Models.(fldnm)(81:96,81:96) = [0,NaN,1,1,1,1,1,1,1,1,1,1,1,1,1,1;NaN,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,0,NaN,1,1,1,1,1,1,1,1,1,1,1,1;1,1,NaN,0,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,0,NaN,1,1,1,1,1,1,1,1,1,1;1,1,1,1,NaN,0,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,0,NaN,1,1,1,1,1,1,1,1;1,1,1,1,1,1,NaN,0,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,0,NaN,1,1,1,1,1,1;1,1,1,1,1,1,1,1,NaN,0,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,0,NaN,1,1,1,1;1,1,1,1,1,1,1,1,1,1,NaN,0,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,0,NaN,1,1;1,1,1,1,1,1,1,1,1,1,1,1,NaN,0,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,NaN;1,1,1,1,1,1,1,1,1,1,1,1,1,1,NaN,0]; % Syl2 Model
Models.(fldnm)(65:80,:) = NaN; % NoiseFinal Strong
Models.(fldnm)(:,65:80) = NaN;
Models.(fldnm)(1:64,:) = NaN; % % Clear Strong/Weak / Matching/Mismatching
Models.(fldnm)(:,1:64) = NaN;
Models.(fldnm)(97:112,:) = NaN; % NoiseInitial
Models.(fldnm)(:,97:112) = NaN;
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = NaN;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'Syl2NoiseSpeech';
Models.(fldnm)(97:112,97:112) = [0,NaN,1,1,1,1,1,1,1,1,1,1,1,1,1,1;NaN,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1;1,1,0,NaN,1,1,1,1,1,1,1,1,1,1,1,1;1,1,NaN,0,1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,0,NaN,1,1,1,1,1,1,1,1,1,1;1,1,1,1,NaN,0,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,0,NaN,1,1,1,1,1,1,1,1;1,1,1,1,1,1,NaN,0,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,0,NaN,1,1,1,1,1,1;1,1,1,1,1,1,1,1,NaN,0,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,0,NaN,1,1,1,1;1,1,1,1,1,1,1,1,1,1,NaN,0,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1,0,NaN,1,1;1,1,1,1,1,1,1,1,1,1,1,1,NaN,0,1,1;1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,NaN;1,1,1,1,1,1,1,1,1,1,1,1,1,1,NaN,0];
Models.(fldnm)(1:64,:) = NaN; % % Clear Strong/Weak / Matching/Mismatching
Models.(fldnm)(:,1:64) = NaN;
Models.(fldnm)(65:80,:) = NaN; % NoiseFinal Strong
Models.(fldnm)(:,65:80) = NaN;
Models.(fldnm)(81:96,:) = NaN; % NoiseFinal Weak
Models.(fldnm)(:,81:96) = NaN;
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = NaN;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'Syl2Graded';
load('RDM_phon_syl2.mat');
Models.(fldnm)(1:64,1:64) = RDM(1:2:128,1:2:128);
Models.(fldnm)(isnan(Models.Syl2(1:64,1:64))) = NaN;
Models.(fldnm)(65:80,:) = NaN; % NoiseFinal Strong
Models.(fldnm)(:,65:80) = NaN;
Models.(fldnm)(81:96,:) = NaN; % NoiseFinal Weak
Models.(fldnm)(:,81:96) = NaN;
Models.(fldnm)(97:112,97:112) = Models.(fldnm)(1:16,1:16); % NoiseInitial
Models.(fldnm)(97:112,1:96) = Models.(fldnm)(1:16,1:96);
Models.(fldnm)(1:96,97:112) = Models.(fldnm)(1:96,1:16);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'Syl2GradedStrongM';
load('RDM_phon_syl2.mat');
Models.(fldnm)(1:16,1:16) = RDM(1:2:32,1:2:32);
Models.(fldnm)(isnan(Models.Syl2(1:16,1:16))) = NaN;
Models.(fldnm) = kron([1 NaN; NaN NaN],Models.(fldnm)); % Extend across Strong and Weak
Models.(fldnm) = kron([1 NaN; NaN NaN],Models.(fldnm)); % Extend across Matching and Mismatching
Models.(fldnm)(65:80,:) = NaN; % NoiseFinal Strong
Models.(fldnm)(:,65:80) = NaN;
Models.(fldnm)(81:96,:) = NaN; % NoiseFinal Weak
Models.(fldnm)(:,81:96) = NaN;
Models.(fldnm)(97:112,:) = NaN; % NoiseInitial
Models.(fldnm)(:,97:112) = NaN;
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'Syl2GradedWeakM';
load('RDM_phon_syl2.mat');
Models.(fldnm)(1:16,1:16) = RDM(1:2:32,1:2:32);
Models.(fldnm)(isnan(Models.Syl2(1:16,1:16))) = NaN;
Models.(fldnm) = kron([NaN NaN; NaN 1],Models.(fldnm)); % Extend across Strong and Weak
Models.(fldnm) = kron([1 NaN; NaN NaN],Models.(fldnm)); % Extend across Matching and Mismatching
Models.(fldnm)(65:80,:) = NaN; % NoiseFinal Strong
Models.(fldnm)(:,65:80) = NaN;
Models.(fldnm)(81:96,:) = NaN; % NoiseFinal WeakGradedGraded
Models.(fldnm)(:,81:96) = NaN;
Models.(fldnm)(97:112,:) = NaN; % NoiseInitial
Models.(fldnm)(:,97:112) = NaN;
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'Syl2GradedStrongMM';
load('RDM_phon_syl2.mat');
Models.(fldnm)(1:16,1:16) = RDM(1:2:32,1:2:32);
Models.(fldnm)(isnan(Models.Syl2(1:16,1:16))) = NaN;
Models.(fldnm) = kron([1 NaN; NaN NaN],Models.(fldnm)); % Extend across Strong and Weak
Models.(fldnm) = kron([NaN NaN; NaN 1],Models.(fldnm)); % Extend across Matching and Mismatching
Models.(fldnm)(65:80,:) = NaN; % NoiseFinal Strong
Models.(fldnm)(:,65:80) = NaN;
Models.(fldnm)(81:96,:) = NaN; % NoiseFinal Weak
Models.(fldnm)(:,81:96) = NaN;
Models.(fldnm)(97:112,:) = NaN; % NoiseInitial
Models.(fldnm)(:,97:112) = NaN;
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'Syl2GradedWeakMM';
load('RDM_phon_syl2.mat');
Models.(fldnm)(1:16,1:16) = RDM(1:2:32,1:2:32);
Models.(fldnm)(isnan(Models.Syl2(1:16,1:16))) = NaN;
Models.(fldnm) = kron([NaN NaN; NaN 1],Models.(fldnm)); % Extend across Strong and Weak
Models.(fldnm) = kron([NaN NaN; NaN 1],Models.(fldnm)); % Extend across Matching and Mismatching
Models.(fldnm)(65:80,:) = NaN; % NoiseFinal Strong
Models.(fldnm)(:,65:80) = NaN;
Models.(fldnm)(81:96,:) = NaN; % NoiseFinal Weak
Models.(fldnm)(:,81:96) = NaN;
Models.(fldnm)(97:112,:) = NaN; % NoiseInitial
Models.(fldnm)(:,97:112) = NaN;
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'Syl2GradedStrongNoise';
load('RDM_phon_syl2.mat');
Models.(fldnm)(65:80,65:80) = RDM(1:2:32,1:2:32); % Syl2 Model
Models.(fldnm)(isnan(Models.Syl2(65:80,65:80))) = NaN;
Models.(fldnm)(81:96,:) = NaN; % NoiseFinal Weak
Models.(fldnm)(:,81:96) = NaN;
Models.(fldnm)(1:64,:) = NaN; % % Clear Strong/Weak / Matching/Mismatching
Models.(fldnm)(:,1:64) = NaN;
Models.(fldnm)(97:112,:) = NaN; % NoiseInitial
Models.(fldnm)(:,97:112) = NaN;
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'Syl2GradedWeakNoise';
load('RDM_phon_syl2.mat');
Models.(fldnm)(81:96,81:96) = RDM(1:2:32,1:2:32); % Syl2 Model
Models.(fldnm)(isnan(Models.Syl2(81:96,81:96))) = NaN;
Models.(fldnm)(65:80,:) = NaN; % NoiseFinal Strong
Models.(fldnm)(:,65:80) = NaN;
Models.(fldnm)(1:64,:) = NaN; % % Clear Strong/Weak / Matching/Mismatching
Models.(fldnm)(:,1:64) = NaN;
Models.(fldnm)(97:112,:) = NaN; % NoiseInitial
Models.(fldnm)(:,97:112) = NaN;
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'Syl2GradedNoiseSpeech';
load('RDM_phon_syl2.mat');
Models.(fldnm)(97:112,97:112) = RDM(1:2:32,1:2:32); % Syl2 Model
Models.(fldnm)(isnan(Models.Syl2(97:112,97:112))) = NaN;
Models.(fldnm)(1:64,:) = NaN; % % Clear Strong/Weak / Matching/Mismatching
Models.(fldnm)(:,1:64) = NaN;
Models.(fldnm)(65:80,:) = NaN; % NoiseFinal Strong
Models.(fldnm)(:,65:80) = NaN;
Models.(fldnm)(81:96,:) = NaN; % NoiseFinal Weak
Models.(fldnm)(:,81:96) = NaN;
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);