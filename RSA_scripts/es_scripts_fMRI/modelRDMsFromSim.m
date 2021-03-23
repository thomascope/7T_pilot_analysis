function Models = modelRDMs()

%%

fldnm = 'PEabsSyl1';
load('RDM_sim_peAbs_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:128,1:128) = RDM(1:128,1:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PEabsSyl1Within';
load('RDM_sim_peAbs_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:64,1:64) = RDM(1:64,1:64);
Models.(fldnm)(65:128,65:128) = RDM(65:128,65:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PEabsSyl1StrongM';
load('RDM_sim_peAbs_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:32,1:32) = RDM(1:32,1:32);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PEabsSyl1WeakM';
load('RDM_sim_peAbs_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(33:64,33:64) = RDM(33:64,33:64);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PEabsSyl1StrongMM';
load('RDM_sim_peAbs_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(65:96,65:96) = RDM(65:96,65:96);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PEabsSyl1WeakMM';
load('RDM_sim_peAbs_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(97:128,97:128) = RDM(97:128,97:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PEabsSyl1StrongNoise';
load('RDM_sim_peAbs_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(129:160,129:160) = RDM(129:160,129:160);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PEabsSyl1WeakNoise';
load('RDM_sim_peAbs_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(161:192,161:192) = RDM(161:192,161:192);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PEabsSyl2';
load('RDM_sim_peAbs_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:128,1:128) = RDM(1:128,1:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PEabsSyl2Within';
load('RDM_sim_peAbs_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:64,1:64) = RDM(1:64,1:64);
Models.(fldnm)(65:128,65:128) = RDM(65:128,65:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PEabsSyl2StrongM';
load('RDM_sim_peAbs_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:32,1:32) = RDM(1:32,1:32);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PEabsSyl2WeakM';
load('RDM_sim_peAbs_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(33:64,33:64) = RDM(33:64,33:64);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PEabsSyl2StrongMM';
load('RDM_sim_peAbs_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(65:96,65:96) = RDM(65:96,65:96);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PEabsSyl2WeakMM';
load('RDM_sim_peAbs_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(97:128,97:128) = RDM(97:128,97:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PEabsSyl2StrongNoise';
load('RDM_sim_peAbs_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(129:160,129:160) = RDM(129:160,129:160);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PEabsSyl2WeakNoise';
load('RDM_sim_peAbs_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(161:192,161:192) = RDM(161:192,161:192);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PEabsWhole';
load('RDM_sim_peAbs_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:128,1:128) = RDM(1:128,1:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PEabsWholeWithin';
load('RDM_sim_peAbs_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:64,1:64) = RDM(1:64,1:64);
Models.(fldnm)(65:128,65:128) = RDM(65:128,65:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PEabsWholeStrongM';
load('RDM_sim_peAbs_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:32,1:32) = RDM(1:32,1:32);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PEabsWholeWeakM';
load('RDM_sim_peAbs_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(33:64,33:64) = RDM(33:64,33:64);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PEabsWholeStrongMM';
load('RDM_sim_peAbs_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(65:96,65:96) = RDM(65:96,65:96);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PEabsWholeWeakMM';
load('RDM_sim_peAbs_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(97:128,97:128) = RDM(97:128,97:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PEabsWholeStrongNoise';
load('RDM_sim_peAbs_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(129:160,129:160) = RDM(129:160,129:160);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PEabsWholeWeakNoise';
load('RDM_sim_peAbs_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(161:192,161:192) = RDM(161:192,161:192);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

%%

fldnm = 'SensorySyl1';
load('RDM_sim_sensory_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:128,1:128) = RDM(1:128,1:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SensorySyl1Within';
load('RDM_sim_sensory_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:64,1:64) = RDM(1:64,1:64);
Models.(fldnm)(65:128,65:128) = RDM(65:128,65:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SensorySyl1StrongM';
load('RDM_sim_sensory_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:32,1:32) = RDM(1:32,1:32);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SensorySyl1WeakM';
load('RDM_sim_sensory_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(33:64,33:64) = RDM(33:64,33:64);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SensorySyl1StrongMM';
load('RDM_sim_sensory_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(65:96,65:96) = RDM(65:96,65:96);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SensorySyl1WeakMM';
load('RDM_sim_sensory_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(97:128,97:128) = RDM(97:128,97:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SensorySyl1StrongNoise';
load('RDM_sim_sensory_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(129:160,129:160) = RDM(129:160,129:160);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SensorySyl1WeakNoise';
load('RDM_sim_sensory_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(161:192,161:192) = RDM(161:192,161:192);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SensorySyl2';
load('RDM_sim_sensory_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:128,1:128) = RDM(1:128,1:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SensorySyl2Within';
load('RDM_sim_sensory_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:64,1:64) = RDM(1:64,1:64);
Models.(fldnm)(65:128,65:128) = RDM(65:128,65:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SensorySyl2StrongM';
load('RDM_sim_sensory_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:32,1:32) = RDM(1:32,1:32);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SensorySyl2WeakM';
load('RDM_sim_sensory_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(33:64,33:64) = RDM(33:64,33:64);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SensorySyl2StrongMM';
load('RDM_sim_sensory_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(65:96,65:96) = RDM(65:96,65:96);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SensorySyl2WeakMM';
load('RDM_sim_sensory_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(97:128,97:128) = RDM(97:128,97:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SensorySyl2StrongNoise';
load('RDM_sim_sensory_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(129:160,129:160) = RDM(129:160,129:160);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SensorySyl2WeakNoise';
load('RDM_sim_sensory_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(161:192,161:192) = RDM(161:192,161:192);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SensoryWhole';
load('RDM_sim_sensory_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:128,1:128) = RDM(1:128,1:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SensoryWholeWithin';
load('RDM_sim_sensory_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:64,1:64) = RDM(1:64,1:64);
Models.(fldnm)(65:128,65:128) = RDM(65:128,65:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SensoryWholeStrongM';
load('RDM_sim_sensory_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:32,1:32) = RDM(1:32,1:32);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SensoryWholeWeakM';
load('RDM_sim_sensory_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(33:64,33:64) = RDM(33:64,33:64);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SensoryWholeStrongMM';
load('RDM_sim_sensory_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(65:96,65:96) = RDM(65:96,65:96);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SensoryWholeWeakMM';
load('RDM_sim_sensory_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(97:128,97:128) = RDM(97:128,97:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SensoryWholeStrongNoise';
load('RDM_sim_sensory_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(129:160,129:160) = RDM(129:160,129:160);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SensoryWholeWeakNoise';
load('RDM_sim_sensory_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(161:192,161:192) = RDM(161:192,161:192);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

%%

fldnm = 'PredSyl1';
load('RDM_sim_pred_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:128,1:128) = RDM(1:128,1:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PredSyl1Within';
load('RDM_sim_pred_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:64,1:64) = RDM(1:64,1:64);
Models.(fldnm)(65:128,65:128) = RDM(65:128,65:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PredSyl1StrongM';
load('RDM_sim_pred_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:32,1:32) = RDM(1:32,1:32);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PredSyl1WeakM';
load('RDM_sim_pred_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(33:64,33:64) = RDM(33:64,33:64);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PredSyl1StrongMM';
load('RDM_sim_pred_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(65:96,65:96) = RDM(65:96,65:96);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PredSyl1WeakMM';
load('RDM_sim_pred_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(97:128,97:128) = RDM(97:128,97:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PredSyl1StrongNoise';
load('RDM_sim_pred_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(129:160,129:160) = RDM(129:160,129:160);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PredSyl1WeakNoise';
load('RDM_sim_pred_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(161:192,161:192) = RDM(161:192,161:192);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PredSyl2';
load('RDM_sim_pred_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:128,1:128) = RDM(1:128,1:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PredSyl2Within';
load('RDM_sim_pred_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:64,1:64) = RDM(1:64,1:64);
Models.(fldnm)(65:128,65:128) = RDM(65:128,65:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PredSyl2StrongM';
load('RDM_sim_pred_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:32,1:32) = RDM(1:32,1:32);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PredSyl2WeakM';
load('RDM_sim_pred_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(33:64,33:64) = RDM(33:64,33:64);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PredSyl2StrongMM';
load('RDM_sim_pred_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(65:96,65:96) = RDM(65:96,65:96);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PredSyl2WeakMM';
load('RDM_sim_pred_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(97:128,97:128) = RDM(97:128,97:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PredSyl2StrongNoise';
load('RDM_sim_pred_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(129:160,129:160) = RDM(129:160,129:160);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PredSyl2WeakNoise';
load('RDM_sim_pred_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(161:192,161:192) = RDM(161:192,161:192);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PredWhole';
load('RDM_sim_pred_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:128,1:128) = RDM(1:128,1:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PredWholeWithin';
load('RDM_sim_pred_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:64,1:64) = RDM(1:64,1:64);
Models.(fldnm)(65:128,65:128) = RDM(65:128,65:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PredWholeStrongM';
load('RDM_sim_pred_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:32,1:32) = RDM(1:32,1:32);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PredWholeWeakM';
load('RDM_sim_pred_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(33:64,33:64) = RDM(33:64,33:64);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PredWholeStrongMM';
load('RDM_sim_pred_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(65:96,65:96) = RDM(65:96,65:96);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PredWholeWeakMM';
load('RDM_sim_pred_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(97:128,97:128) = RDM(97:128,97:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PredWholeStrongNoise';
load('RDM_sim_pred_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(129:160,129:160) = RDM(129:160,129:160);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PredWholeWeakNoise';
load('RDM_sim_pred_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(161:192,161:192) = RDM(161:192,161:192);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

%%

fldnm = 'SsSyl1';
load('RDM_sim_ss_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:128,1:128) = RDM(1:128,1:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SsSyl1Within';
load('RDM_sim_ss_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:64,1:64) = RDM(1:64,1:64);
Models.(fldnm)(65:128,65:128) = RDM(65:128,65:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SsSyl1StrongM';
load('RDM_sim_ss_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:32,1:32) = RDM(1:32,1:32);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SsSyl1WeakM';
load('RDM_sim_ss_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(33:64,33:64) = RDM(33:64,33:64);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SsSyl1StrongMM';
load('RDM_sim_ss_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(65:96,65:96) = RDM(65:96,65:96);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SsSyl1WeakMM';
load('RDM_sim_ss_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(97:128,97:128) = RDM(97:128,97:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SsSyl1StrongNoise';
load('RDM_sim_ss_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(129:160,129:160) = RDM(129:160,129:160);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SsSyl1WeakNoise';
load('RDM_sim_ss_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(161:192,161:192) = RDM(161:192,161:192);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SsSyl2';
load('RDM_sim_ss_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:128,1:128) = RDM(1:128,1:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SsSyl2Within';
load('RDM_sim_ss_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:64,1:64) = RDM(1:64,1:64);
Models.(fldnm)(65:128,65:128) = RDM(65:128,65:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SsSyl2StrongM';
load('RDM_sim_ss_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:32,1:32) = RDM(1:32,1:32);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SsSyl2WeakM';
load('RDM_sim_ss_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(33:64,33:64) = RDM(33:64,33:64);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SsSyl2StrongMM';
load('RDM_sim_ss_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(65:96,65:96) = RDM(65:96,65:96);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SsSyl2WeakMM';
load('RDM_sim_ss_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(97:128,97:128) = RDM(97:128,97:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SsSyl2StrongNoise';
load('RDM_sim_ss_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(129:160,129:160) = RDM(129:160,129:160);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SsSyl2WeakNoise';
load('RDM_sim_ss_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(161:192,161:192) = RDM(161:192,161:192);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SsWhole';
load('RDM_sim_ss_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:128,1:128) = RDM(1:128,1:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SsWholeWithin';
load('RDM_sim_ss_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:64,1:64) = RDM(1:64,1:64);
Models.(fldnm)(65:128,65:128) = RDM(65:128,65:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SsWholeStrongM';
load('RDM_sim_ss_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:32,1:32) = RDM(1:32,1:32);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SsWholeWeakM';
load('RDM_sim_ss_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(33:64,33:64) = RDM(33:64,33:64);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SsWholeStrongMM';
load('RDM_sim_ss_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(65:96,65:96) = RDM(65:96,65:96);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SsWholeWeakMM';
load('RDM_sim_ss_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(97:128,97:128) = RDM(97:128,97:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SsWholeStrongNoise';
load('RDM_sim_ss_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(129:160,129:160) = RDM(129:160,129:160);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'SsWholeWeakNoise';
load('RDM_sim_ss_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(161:192,161:192) = RDM(161:192,161:192);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

%%

fldnm = 'PsSyl1';
load('RDM_sim_ps_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:128,1:128) = RDM(1:128,1:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PsSyl1Within';
load('RDM_sim_ps_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:64,1:64) = RDM(1:64,1:64);
Models.(fldnm)(65:128,65:128) = RDM(65:128,65:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PsSyl1StrongM';
load('RDM_sim_ps_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:32,1:32) = RDM(1:32,1:32);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PsSyl1WeakM';
load('RDM_sim_ps_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(33:64,33:64) = RDM(33:64,33:64);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PsSyl1StrongMM';
load('RDM_sim_ps_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(65:96,65:96) = RDM(65:96,65:96);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PsSyl1WeakMM';
load('RDM_sim_ps_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(97:128,97:128) = RDM(97:128,97:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PsSyl1StrongNoise';
load('RDM_sim_ps_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(129:160,129:160) = RDM(129:160,129:160);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PsSyl1WeakNoise';
load('RDM_sim_ps_syl1.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(161:192,161:192) = RDM(161:192,161:192);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PsSyl2';
load('RDM_sim_ps_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:128,1:128) = RDM(1:128,1:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PsSyl2Within';
load('RDM_sim_ps_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:64,1:64) = RDM(1:64,1:64);
Models.(fldnm)(65:128,65:128) = RDM(65:128,65:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PsSyl2StrongM';
load('RDM_sim_ps_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:32,1:32) = RDM(1:32,1:32);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PsSyl2WeakM';
load('RDM_sim_ps_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(33:64,33:64) = RDM(33:64,33:64);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PsSyl2StrongMM';
load('RDM_sim_ps_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(65:96,65:96) = RDM(65:96,65:96);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PsSyl2WeakMM';
load('RDM_sim_ps_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(97:128,97:128) = RDM(97:128,97:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PsSyl2StrongNoise';
load('RDM_sim_ps_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(129:160,129:160) = RDM(129:160,129:160);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PsSyl2WeakNoise';
load('RDM_sim_ps_syl2.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(161:192,161:192) = RDM(161:192,161:192);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PsWhole';
load('RDM_sim_ps_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:128,1:128) = RDM(1:128,1:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PsWholeWithin';
load('RDM_sim_ps_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:64,1:64) = RDM(1:64,1:64);
Models.(fldnm)(65:128,65:128) = RDM(65:128,65:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PsWholeStrongM';
load('RDM_sim_ps_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(1:32,1:32) = RDM(1:32,1:32);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PsWholeWeakM';
load('RDM_sim_ps_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(33:64,33:64) = RDM(33:64,33:64);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PsWholeStrongMM';
load('RDM_sim_ps_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(65:96,65:96) = RDM(65:96,65:96);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PsWholeWeakMM';
load('RDM_sim_ps_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(97:128,97:128) = RDM(97:128,97:128);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PsWholeStrongNoise';
load('RDM_sim_ps_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(129:160,129:160) = RDM(129:160,129:160);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

fldnm = 'PsWholeWeakNoise';
load('RDM_sim_ps_whole.mat');
Models.(fldnm) = NaN(224,224);
Models.(fldnm)(161:192,161:192) = RDM(161:192,161:192);
Models.(fldnm)(find(eye(size(Models.(fldnm))))) = NaN; % Make diagonal NaN
imAlpha = ones(size(Models.(fldnm)));
imAlpha(isnan(Models.(fldnm))) = 0;
figure; imagesc(Models.(fldnm),'AlphaData',imAlpha);

%%

close all
r = corr([Models.SensorySyl1(:) Models.PredSyl1(:) Models.PEabsSyl1(:) Models.PsSyl1(:)...
    Models.SensorySyl2(:) Models.PredSyl2(:) Models.PEabsSyl2(:) Models.PsSyl2(:)],'rows','pairwise');
figure; imagesc(r);