function v = RelaxWJ(A, v, f, w, n)
% Weighted Jacobi relaxation n times for Av = f

D = diag(diag(A));
LU = D - A;

for indexTimes = 1 : n
    v = (1 - w) * v + w * sparse(diag(1./diag(D))) * (LU * v + f);
end

