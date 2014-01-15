clear all;
close all;
clc;

load('boundary.mat');
load('mesh.mat');
load('coeff.mat');

[K, F]=assempde(b,p,e,t,c,a,f);

theta = 0.25; % Threshold to determine strong dependency
[flagC, S, St] = GetFC(K, theta);  % Determine fine points, coarse points as well as influence and dependence

numPoint = length(flagC);
point = 1 : numPoint;
pointC = point(flagC);
pointF = point(~flagC);

maxRow = zeros(numPoint, 1);
for indexPoint = 1 : numPoint
    rowTemp = [K(indexPoint, 1 : indexPoint - 1),  K(indexPoint, indexPoint + 1 : numPoint)];
    maxRow(indexPoint) = max(rowTemp);
end
% load('mesh.mat');
% figure;
% pdemesh(p, e, t), hold on;
% plot(p(1, pointC), p(2, pointC), 'ro', 'linewidth', 2, 'MarkerSize', 10);

uAMG = zeros(numPoint, 1); % Initial guess
w = 2 / 3; % Weighting coefficients for weigted Jacobi relaxation
n1 = 3; % times of relaxation before moving down to coarser grid
n2 = 3; % times of relaxation after moving up from coarser grid
maxDepth = 3; % Maximum depth (times of recursion)ma
theta = 0.25; % Threshold to determine strong dependence

numVCycle = 10;
normErr = zeros(numVCycle, 1);
normRes = zeros(numVCycle, 1);
u = assempde(b,p,e,t,c,a,f); % Solution given by MATLAB solver
% for indexVCycle = 1 : numVCycle
%     disp(['The ', num2str(indexVCycle), ' -th V Cycle']);
%     uAMG = amgUD(K, uAMG, F, w, n1, n2, 0, maxDepth, theta, 1 : numPoint, false); % Solution given by AMG
%     normErr(indexVCycle) = norm(u - uAMG);
%     normRes(indexVCycle) = norm(K * uAMG - F);
% end

% figure;
% semilogy(1 : numVCycle, normErr, 'b+-', 'linewidth', 2);
% grid on, xlim([1, numVCycle]);
% set(gca, 'Fontsize', 16), xlabel('Number of V-cycle'), ylabel('||u - u_{AMG}||');
% 
% figure;
% semilogy(1 : numVCycle, normRes, 'b+-', 'linewidth', 2);
% grid on, xlim([1, numVCycle]);
% set(gca, 'Fontsize', 16), xlabel('Number of V-cycle'), ylabel('||Au_{AMG}-f||');

u = K \ F;
disp(norm(K * u - F));
% figure;
% subplot(1, 2, 1);
% pdeplot(p, [], t, 'xydata', u, 'xystyle', 'interp', 'zdata', u, 'zstyle', 'continuous', 'colormap', jet), colorbar, caxis([0, 1]);
% grid on, set(gca, 'Fontsize', 16), title('Actual solution')
% subplot(1, 2, 2);
% pdeplot(p, [], t, 'xydata', uAMG, 'xystyle', 'interp', 'zdata', uAMG, 'zstyle', 'continuous', 'colormap', jet), colorbar, caxis([0, 1]);
% grid on, set(gca, 'Fontsize', 16), title(['AMG solution, ', num2str(numVCycle), ' V-cycle']);
% figure;
% pdeplot(p, [], t, 'xydata', u - uAMG, 'xystyle', 'interp', 'zdata', u - uAMG, 'zstyle', 'continuous', 'colormap', jet);
% grid on, set(gca, 'Fontsize', 16), title(['Difference between AMG and actual solution, ', num2str(numVCycle), ' V-cycle']);
