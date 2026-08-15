%% =========================================================
% SUPERCRITICAL SLAB:
% TRANSFORMED-VARIABLE CLASSICAL SHOOTING
% =========================================================
clear
clc

cases = [ ...
    0.50,  6.0000;
    0.90, 29.2404];

opts.RelTol     = 1e-10;
opts.AbsTol     = 1e-12;
opts.tolA       = 1e-10;
opts.eps        = 1e-5;
opts.maxIter    = 100;
opts.maxBracket = 100;

nRepeat = 10;

fprintf('\n');
fprintf('Supercritical slab: transformed-variable shooting\n');
fprintf('===================================================\n');

for k = 1:size(cases,1)

    n   = cases(k,1);
    phi = cases(k,2);

    % Analytical value used ONLY for validation/error reporting.
    phiStar = sqrt(2*(1+n))/(1-n);
    xdzRef  = 1 - phiStar/phi;

    result = solve_supercritical_slab_shooting_y(phi,n,opts);

    xdzShoot = result.xdz;
    errXdz   = abs(xdzShoot-xdzRef);

    % Timing.
    times = zeros(nRepeat,1);

    for j = 1:nRepeat
        tic
        solve_supercritical_slab_shooting_y(phi,n,opts);
        times(j) = toc;
    end

    fprintf('\nCase: n = %.2f, phi = %.4f\n',n,phi);
    fprintf('q                       = %.6f\n',result.q);
    fprintf('phi* (validation)       = %.12f\n',phiStar);
    fprintf('x_dz shooting           = %.12f\n',xdzShoot);
    fprintf('x_dz reference          = %.12f\n',xdzRef);
    fprintf('|Delta x_dz|            = %.6e\n',errXdz);
    fprintf('iterations              = %d\n',result.iterations);
    fprintf('epsilon                 = %.1e\n',result.eps);
    fprintf('mean time (%d runs)     = %.6f s\n',nRepeat,mean(times));
    fprintf('time std. dev.          = %.6f s\n',std(times));
end

fprintf('\n');
