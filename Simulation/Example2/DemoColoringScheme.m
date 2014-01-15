clear all;
close all;
clc;

%% 1 Parameter setup
load('boundary.mat');
load('mesh.mat');
load('coeff.mat');

[K, F]=assempde(b,p,e,t,c,a,f);
theta = 0.25; % Threshold to determine strong dependency
[flagC, S, St] = GetFC(K, theta);  % Determine fine points, coarse points as well as influence and dependence

%indexPoint = [57, 624, 1896, 1138, 918, 43, 11, 2]; % The index of the point
% indexPoint = [5, 8, 2, 10]; % The index of the point
% figure;
% pdemesh(p, e, t), hold on;
% plot(p(1, indexPoint), p(2, indexPoint), 'k+', 'linewidth', 2);
% pointS = zeros(1, 0);
% pointSt = zeros(1, 0);
% for indexTemp = 1 : length(indexPoint)
%     pointS = [pointS, S{indexPoint(indexTemp)}];
%     pointSt = [pointSt, St{indexPoint(indexTemp)}];
% end
% plot(p(1, pointS), p(2, pointS), 'ro', 'linewidth', 2);
% plot(p(1, pointSt), p(2, pointSt), 'g+', 'linewidth', 2);
% axis equal;
% set(gca, 'Fontsize', 16);
% legend('Edge', 'Boundary', 'Current Point', 'Points influencing', 'Points influenced');
% text(1.1, 1.6, 'a=1', 'Fontsize', 16, 'Fontweight', 'bold');
% text(4.5, 1.6, 'a=1000', 'Fontsize', 16, 'Fontweight', 'bold');
% text(8.0, 1.6, 'a=1', 'Fontsize', 16, 'Fontweight', 'bold');
% title('-\nabla(a(x,y)\nablau)=f');
% [C, I] =min((p(1, :) - 1.75).^2 + (p(2, :) - 2.75).^2);

numPoint = size(p, 2);
flagCata = zeros(numPoint, 1); % The flag value indicating the category of each points: 0 - unassigned, 1 - coarse point, 2 - fine point
lambda = zeros(numPoint, 1);  % The value indicating how likely each point is a C point, 
for indexPoint = 1 : numPoint %the i-th entry initialized as the cardinality of St_i
    lambda(indexPoint) = length(St{indexPoint});
end
flagUnassign = true; 
figure;
pdemesh(p, e, t), hold on;
axis equal, xlim([0, 10]), ylim([-0.5, 6]);
set(gca, 'Fontsize', 16);
% text(1.1, 1.6, 'a=1', 'Fontsize', 16, 'Fontweight', 'bold');
% text(4.5, 1.6, 'a=1000', 'Fontsize', 16, 'Fontweight', 'bold');
% text(8.0, 1.6, 'a=1', 'Fontsize', 16, 'Fontweight', 'bold');
title('-\nabla(a(x,y)\nablau)=f');
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
    pointIncre = [];
    for indexPointF = 1 : length(pointFTemp)
        pointIncreTemp = S{pointFTemp(indexPointF)};
        pointIncreTemp = pointIncreTemp(~isnan(lambda(pointIncreTemp)));
        pointIncreTemp = pointIncreTemp(flagCata(pointIncreTemp) == 0);
        lambda(pointIncreTemp) = lambda(pointIncreTemp) + 1;
        pointIncre = [pointIncre, pointIncreTemp];
    end
    
    % Check whether there are unassigned point
    flagUnassign = (min(flagCata) == 0);
    
    if (true)
%         point = 1 : numPoint;
%         pointC = point(flagCata == 1);
%         pointF = point(flagCata == 2);
% 
%         plot(p(1, pointC), p(2, pointC), 'ro', 'linewidth', 2), hold on;
%         plot(p(1, pointF), p(2, pointF), 'go', 'linewidth', 2), hold on;
%         plot(p(1, pointIncre), p(2, pointIncre), 'mo', 'linewidth', 2), hold on;
%         legend('Edge', 'Boundary', 'Coarse', 'Fine', 'Increment');        
        %pause;
    end
end

point = 1 : numPoint;
pointC = point(flagCata == 1);
pointF = point(flagCata == 2);

plot(p(1, pointC), p(2, pointC), 'ro', 'linewidth', 2), hold on;
legend('Edge', 'Boundary', 'Coarse');  
text(1.1, 1.6, 'a=1', 'Fontsize', 16, 'Fontweight', 'bold');
text(4.5, 1.6, 'a=1000', 'Fontsize', 16, 'Fontweight', 'bold');
text(8.0, 1.6, 'a=1', 'Fontsize', 16, 'Fontweight', 'bold');

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

% flagC = (flagCata == 1);
% pointC = point(flagCata == 1);
% pointF = point(flagCata == 2);
% 
% plot(p(1, pointC), p(2, pointC), 'ro', 'linewidth', 2), hold on;
% legend('Edge', 'Boundary', 'Coarse');  
% text(1.1, 1.6, 'a=1', 'Fontsize', 16, 'Fontweight', 'bold');
% text(4.5, 1.6, 'a=1000', 'Fontsize', 16, 'Fontweight', 'bold');
% text(8.0, 1.6, 'a=1', 'Fontsize', 16, 'Fontweight', 'bold');



