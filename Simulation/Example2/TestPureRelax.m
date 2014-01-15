clear all;
close all;
clc;

load('boundary.mat');
load('mesh.mat');
load('coeff.mat');

[K, F]=assempde(b,p,e,t,c,a,f);

u = zeros(size(K, 1), 1); % Initial guess
w = 2 / 3; % Weighting coefficient for weighted Jacobi method
numRelax = 1000; % Maximum times of relaxation
normRes = zeros(numRelax, 1);
for indexRelax = 1 : numRelax
    u = RelaxWJ(K, u, F, w, 1); % Relax n2 times
    normRes(indexRelax) = norm(K * u - F);
    clc;
    disp([num2str(indexRelax), '/', num2str(numRelax), ' completed'])
end

figure;
semilogy(1 : numRelax, normRes, 'b-', 'linewidth', 2);
grid on, xlim([1, numRelax]);
set(gca, 'Fontsize', 16), xlabel('Number of relaxation'), ylabel('||Au-f||');