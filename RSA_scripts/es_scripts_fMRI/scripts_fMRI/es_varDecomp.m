function [R2c,R2u] = es_varDecomp(Y,X)

nPred = size(X,2);
nRep = size(X,1);

[~,R2full] = es_regress(Y,[X ones(nRep,1)]);
for i=1:nPred
    indCurrent = setdiff(1:nPred,i);
    [~,R2] = es_regress(Y,[X(:,indCurrent) ones(nRep,1)]);
    R2u(i,:) = R2full - R2;
end
R2c = R2full - sum(R2u,1);
