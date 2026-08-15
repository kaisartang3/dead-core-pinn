%% =========================================================
% FINAL SHOOTING TIMING FOR TABLE
% =========================================================
clear
clc

cases = [ ...
    0.50,  6.0000;
    0.90, 29.2404];

% Set this after checking the epsilon convergence study.
opts.RelTol     = 1e-10;
opts.AbsTol     = 1e-12;
opts.tolA       = 1e-10;
opts.eps        = 1e-4;
opts.maxIter    = 100;
opts.maxBracket = 100;

nRepeat = 20;

fprintf('\n');
fprintf('FINAL SHOOTING RESULTS FOR TABLE\n');
fprintf('=================================\n');
fprintf('n       phi       x_dz            error            time (s)\n');

for k = 1:size(cases,1)

    n   = cases(k,1);
    phi = cases(k,2);

    phiStar = sqrt(2*(1+n))/(1-n);
    xdzRef  = 1 - phiStar/phi;

    % Warm-up.
    solve_supercritical_slab_shooting_y(phi,n,opts);

    times = zeros(nRepeat,1);

    for j = 1:nRepeat
        tic
        result = solve_supercritical_slab_shooting_y(phi,n,opts);
        times(j) = toc;
    end

    fprintf('%.2f   %.4f   %.12f   %.6e   %.6f +/- %.6f\n', ...
        n,phi,result.xdz,abs(result.xdz-xdzRef), ...
        mean(times),std(times));
end
