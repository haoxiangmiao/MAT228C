function [flagC, S, St] = GetFC(spA, theta)
% Get fine point and coarse point
% flagC:     numPoint * 1 logic vector, 1 - Coarse, 0 - Fine
% S:        
% St:       
% spA:      the  numPoint * numPoint LHS matrix A (sparse)
% theta:    the threshold that determines strong influence/dependence

% Initialization
numPoint = size(spA, 1); % number of points in the graph
point = 1 : numPoint;

[S, St] = GetSetInflDepend(spA, theta); % Find the set of points strongly influencing or influenced by each point

flagCata = zeros(numPoint, 1); % The flag value indicating the category of each points: 0 - unassigned, 1 - coarse point, 2 - fine point. 
lambda = zeros(numPoint, 1);  % The value indicating how likely each point is a C point, 
for indexPoint = 1 : numPoint %the i-th entry initialized as the cardinality of St_i
    lambda(indexPoint) = length(St{indexPoint});
end

% The first pass, using H-2 as a guide
flagUnassign = true; 
while flagUnassign % Until all points have been assigned to be C points or F points
    
    [~, pointCTemp] = max(lambda); % Find the first point with the greatest lambda
     
    % Set this point to be a C point
    flagCata(pointCTemp) = 1;
    lambda(pointCTemp) = NaN;
    
    % Set all the unassigned points that are strongly influenced by i to be
    % F point
    StTemp = St{pointCTemp};
    pointFTemp = StTemp(flagCata(StTemp) == 0); % Only find the unassigned points
    flagCata(pointFTemp) = 2;
    lambda(pointFTemp) = NaN;
    
    % For each of the newly found F point j, increase the lambda for each unassigned point in S_j

    for indexPointF = 1 : length(pointFTemp)
        pointIncreTemp = S{pointFTemp(indexPointF)};
        lambda(pointIncreTemp) = lambda(pointIncreTemp) + 1;
    end
    
    % Check whether there are unassigned point
    flagUnassign = (min(flagCata) == 0);
end

% The second pass, rigorously satisfying H-1
pointF = point(flagCata == 2);
numPointF = length(pointF);
for indexPointF = 1 : numPointF - 1 % For each F point i
    if (flagCata(pointF(indexPointF)) == 1) 
        continue;% This point has been modified as a C point
    else       
        Si = S{pointF(indexPointF)}; % Si
        pointFSi = Si(flagCata(Si) == 2); % F points in Si
        pointCSi = Si(flagCata(Si) == 1); % C points in Si
        for indexPointFi = 1 : length(pointFSi) % For each of these F points j,
            if (flagCata(pointFSi(indexPointFi)) == 1)
                continue; % This point has been modified as a C point
            else
                Sj = S{pointFSi(indexPointFi)}; %  test to see whether it has a common C point both in Sj and a C point in Si 
                pointCSj = Sj(flagCata(Sj) == 1);
                if isempty(intersect(pointCSi, pointCSj))
                    flagCata(pointFSi(indexPointFi)) = 1; % If there is no common C point, set F point i to be C point
                    pointCSi = Si(flagCata(Si) == 1); % update C points in Si
                    
%                     flagCata(pointF(indexPointF)) = 1; % If there is no common C point, set F point i to be C point
%                     pointCSi = Si(flagCata(Si) == 1); % update C points in Si
%                     break;
                end
            end
        end
    end
end

flagC = (flagCata == 1);



