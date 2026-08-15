%% =========================================================
% CONVERGENCE OF TRANSFORMED-VARIABLE SHOOTING WITH EPSILON
% =========================================================
clear
clc

cases = [ ...
    0.50,  6.0000;
    0.90, 29.2404];

epsValues = [1e-2 3e-3 1e-3 3e-4 1e-4 3e-5 1e-5 3e-6 1e-6];

opts.RelTol     = 1e-11;
opts.AbsTol     = 1e-13;
opts.tolA       = 1e-11;
opts.maxIter    = 120;
opts.maxBracket = 100;

for k = 1:size(cases,1)

    n   = cases(k,1);
    phi = cases(k,2);

    phiStar = sqrt(2*(1+n))/(1-n);
    xdzRef  = 1 - phiStar/phi;

    fprintf('\n========================================\n');
    fprintf('n = %.2f, phi = %.4f\n',n,phi);
    fprintf('reference x_dz = %.12f\n',xdzRef);
    fprintf('----------------------------------------\n');
    fprintf('epsilon       x_dz(shoot)       error\n');

    for j = 1:length(epsValues)

        opts.eps = epsValues(j);

        result = solve_supercritical_slab_shooting_y(phi,n,opts);

        fprintf('%.1e      %.12f      %.6e\n', ...
            opts.eps,result.xdz,abs(result.xdz-xdzRef));
    end
end
