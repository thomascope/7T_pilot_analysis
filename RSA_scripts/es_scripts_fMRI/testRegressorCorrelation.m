clearvars

load('SPM.mat');

sim = corr(SPM.xX.X);
sim(find(eye(size(sim,1)))) = NaN;

figure; imagesc(sim);

ind = find(sim>.6);
[indA,indB] = ind2sub(size(sim),ind);

clear sim2display
for i=1:length(ind)
    sim2display{i} = num2str(sim(ind(i)));
end

[SPM.xX.name(indA)' SPM.xX.name(indB)' sim2display']

