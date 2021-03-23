function conditionVec = es_getSPMconditionVec(SPMdir,conditions)

load(SPMdir);
conditionVec = zeros(length(SPM.xX.name),1);
for c=1:length(conditions)
    ind = find(not(cellfun('isempty',strfind(SPM.xX.name,conditions{c}))))';
    conditionVec(ind) = c;
end