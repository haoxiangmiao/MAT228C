clear all;
close all;
clc;

%% 1 Parameter setup
load('boundary.mat');
load('mesh.mat');
load('coeff.mat');

[A, F]=assempde(b,p,e,t,c,a,f);
theta = 0.25; % Threshold to determine strong dependency

numPoint = size(p, 2);
point = 1 : numPoint;

maxDepth = 2;
pointC = cell(maxDepth + 1, 1);
pointTemp = point;
ATemp = A;

for depth = 0 : maxDepth   
    [flagC, S, ~] = GetFC(ATemp, theta);  % Determine fine points, coarse points as well as influence and dependence
    pointC{depth + 1} = pointTemp(flagC);
        
    pointTemp = pointTemp(flagC);
    I = GetMatInterp(ATemp, flagC, S); % Get the interpolatory matrix
    ATemp = I' * ATemp * I;
end
    
figure;
pdemesh(p, e, t), hold on;
plot(p(1, pointC{1}), p(2, pointC{1}), 'bo', 'linewidth', 2);
plot(p(1, pointC{2}), p(2, pointC{2}), 'ko', 'linewidth', 2);
plot(p(1, pointC{3}), p(2, pointC{3}), 'ro', 'linewidth', 2);
axis equal;
set(gca, 'Fontsize', 16);
legend('Edge', 'Boundary', ['1st-level: ', num2str(length(pointC{1}))], ['2nd-level: ', num2str(length(pointC{2}))], ['3rd-level: ', num2str(length(pointC{3}))]);
text(1.1, 1.6, 'a=1', 'Fontsize', 16, 'Fontweight', 'bold');
text(4.5, 1.6, 'a=1000', 'Fontsize', 16, 'Fontweight', 'bold');
text(8.0, 1.6, 'a=1', 'Fontsize', 16, 'Fontweight', 'bold');
title('-\nabla(a(x,y)\nablau)=f');
