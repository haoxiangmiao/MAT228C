function I = GetMatInterp(spA, flagC, S)

%% Initialization
numPoint = length(flagC);
point = 1 : numPoint;
pointC = point(flagC);
pointF = point(~flagC);
numPointF = length(pointF);
numPointC = length(pointC);

I = zeros(numPoint, numPointC); % Initialize the interpolatory matrix

%% Interpolatory weight corresponding to the coarse point
for indexPointC = 1 : numPointC
    I(pointC(indexPointC), indexPointC) = 1;
end

for indexPointF = 1 : numPointF
    pointi = pointF(indexPointF); % Current point
    Ni = find(spA(pointi, :)); % The set of all j where aij ~= 0
    Si = S{pointi}; % The set of points that strongly influence i
    Ci = intersect(Ni, Si(flagC(Si))); % The set of coarse-grid points that strongly influence i
    Dsi = intersect(Ni,Si(~flagC(Si))); % The set of fine grid points that strongly influence i
    Dwi = setdiff(Ni, [Si, pointi]); % The set of points that do not strongly influence i

    denomTemp = spA(pointi, pointi) + sum(spA(pointi, Dwi)); % Denominator, same for all wij
    for indexPointCi = 1 : length(Ci) % For each coarse point, compute the numerator and the weight to interolate point i
        pointj = Ci(indexPointCi); 
        numTemp = spA(pointi, pointj);
        for indexPointDsi = 1 : length(Dsi)
            pointm = Dsi(indexPointDsi);
            numTemp = numTemp + spA(pointi, pointm) * spA(pointm, pointj) / sum(spA(pointm, Ci));
        end
        I(pointi, (pointj == pointC)) = - numTemp / denomTemp;
    end
end