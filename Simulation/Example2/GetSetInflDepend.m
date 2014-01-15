function [S, St] = GetSetInflDepend(spA, theta)
% Get {S_i} the set of points that strongly influence i and {St_i} the set of points that strongly depend on i
% S:        numPoint * 1 cell array, the i-th entry is a vector containg all points
%           that strongly influence i
% St:       numPoint * 1 cell array, the i-th entry is a vector containg all points
%           that strongly depend on i
% spA:      the  numPoint * numPoint LHS matrix A (sparse)
% theta:    the threshold that determines strong influence/dependence

numPoint = size(spA, 1); % number of points in the graph

S = cell(numPoint, 1);
for indexPoint = 1 : numPoint
    nonDiagTemp = [spA(indexPoint, 1 : indexPoint - 1), spA(indexPoint, indexPoint + 1 : numPoint)];
    indexNonDiagTemp = [1 : indexPoint - 1, indexPoint + 1 : numPoint];
    maxNonDiagTemp = max(-nonDiagTemp);
    S{indexPoint} = indexNonDiagTemp(-nonDiagTemp >= theta * maxNonDiagTemp);
end

St = cell(numPoint, 1);
for indexPoint = 1 : numPoint
    SiTemp = S{indexPoint};
    for indexS = 1 : length(SiTemp)
        St{SiTemp(indexS)} = [St{SiTemp(indexS)}, indexPoint];
    end
end
