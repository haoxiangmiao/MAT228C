function v = amgUD(A, v, f, w, n1, n2, depth, maxDepth, theta, indexPoint, flagPlot)
% Recursive AMG solving Av = f with one V-cycle
% v: solution
% A: the LHS matrix
% f: the RHS vector
% w: weighting coefficients for weighted Jacobian method
% n1, n2: at each level, relax n1 times before going down to coarser grid, relax n2 times after returning from coarser grid 
% depth: current depth in the V-cycle
% maxDepth: when the V-cycle reaches the maximum depth solve directly
% theta: threshold to determine strong dependence for AMG coarsening
% indexPoint: The index of the points on current level in original points 
% flagPlot: if true, plot the residual and interpolated residual

if (depth < maxDepth)
    v = RelaxWJ(A, v, f, w, n1); % Relax n1 times
    r = f - A * v; % Residual
    
    [flagC, S, ~] = GetFC(A, theta);  % Determine fine points, coarse points as well as influence and dependence
    I = GetMatInterp(A, flagC, S); % Get the interpolatory matrix

    dimC = nnz(flagC); % Number of coarse points 
    rC = I' * r; % restrict the residual to coarser grid
    
    if flagPlot
        load('mesh.mat');
        figure;
        subplot(1, 2, 1);
        tri = delaunay(p(1, indexPoint), p(2, indexPoint));
        trisurf(tri, p(1, indexPoint), p(2, indexPoint), r), colorbar;
        set(gca, 'Fontsize', 16), title(['Residual after pre-relaxation, ',num2str(depth), '-th iteration']);
        subplot(1, 2, 2);
        tri = delaunay(p(1, indexPoint(flagC)), p(2, indexPoint(flagC)));
        trisurf(tri, p(1,indexPoint(flagC)), p(2, indexPoint(flagC)), rC), colorbar;
        set(gca, 'Fontsize', 16), title('Interpolated residual');
    end

    
    AC = I' * A * I; % restrict the LHS matrix to coarser grid
    errC = zeros(dimC, 1); % initial guess of the error
    errC = amgUD(AC, errC, rC, w, n1, n2, depth + 1, maxDepth, theta, indexPoint(flagC), flagPlot); % Go to the coarser grid
    
    err = I * errC; % interpolate to get the error vector from the coarser grid
    vTemp = v + err; % correct v
    v = RelaxWJ(A, vTemp, f, w, n2); % Relax n2 times
    
    if flagPlot

        load('mesh.mat');
        figure;
        subplot(1, 2, 1);
        tri = delaunay(p(1, indexPoint(flagC)), p(2, indexPoint(flagC)));
        trisurf(tri, p(1,indexPoint(flagC)), p(2, indexPoint(flagC)), errC), colorbar;
        set(gca, 'Fontsize', 16), title(['Error, ', num2str(depth), '-th iteration']);
        
        subplot(1, 2, 2);
        tri = delaunay(p(1, indexPoint), p(2, indexPoint));
        trisurf(tri, p(1, indexPoint), p(2, indexPoint), err), colorbar;
        set(gca, 'Fontsize', 16), title('Prolonged error');
        
        rTemp = f - A * vTemp;
        r = f - A * v; % Residual
        figure;
        subplot(1, 2, 1);
        tri = delaunay(p(1, indexPoint), p(2, indexPoint));
        trisurf(tri, p(1,indexPoint), p(2, indexPoint), rTemp), colorbar;
        set(gca, 'Fontsize', 16), title(['Corrected residual, ', num2str(depth), '-th iteration']);
        
        subplot(1, 2, 2);
        tri = delaunay(p(1, indexPoint), p(2, indexPoint));
        trisurf(tri, p(1, indexPoint), p(2, indexPoint), r), colorbar;
        set(gca, 'Fontsize', 16), title('Residual after pro-relaxation');
    end

else % Upon reaching the bottom of the V cycle, solve the linear equation directly
    v = A \ f;
end

disp(['The ', num2str(depth), '/', num2str(maxDepth) ' iteration completed.']);