%% =========================================================
% SUPERCRITICAL SLAB: SHOOTING COMPARISON
% =========================================================
clear
clc

cases = [ ...
    0.50,  6.0000;
    0.90, 29.2404];

opts.RelTol     = 1e-9;
opts.AbsTol     = 1e-11;
opts.tolS       = 1e-12;
opts.maxIter    = 100;
opts.maxBracket = 100;

% Number of repetitions used only for timing.
% Increase this if the individual run is very short.
nRepeat = 10;

fprintf('\nSupercritical slab: classical shooting\n');
fprintf('-----------------------------------------------\n');

for k = 1:size(cases,1)

    n   = cases(k,1);
    phi = cases(k,2);

    % Analytical reference used ONLY for error reporting.
    phiStarRef = sqrt(2*(1+n))/(1-n);
    xdzRef     = 1 - phiStarRef/phi;

    % One run for the numerical result.
    result = solve_supercritical_slab_by_shooting(phi,n,opts);

    xdzShoot = result.xdz;
    errXdz   = abs(xdzShoot-xdzRef);

    % Average timing over repeated runs.
    times = zeros(nRepeat,1);

    for j = 1:nRepeat
        tic
        tmp = solve_supercritical_slab_by_shooting(phi,n,opts); %#ok<NASGU>
        times(j) = toc;
    end

    meanTime = mean(times);
    stdTime  = std(times);

    fprintf('\nCase: n = %.2f, phi = %.4f\n',n,phi);
    fprintf('phi* (reference)       = %.12f\n',phiStarRef);
    fprintf('x_dz (shooting)        = %.12f\n',xdzShoot);
    fprintf('x_dz (reference)       = %.12f\n',xdzRef);
    fprintf('|Delta x_dz|           = %.6e\n',errXdz);
    fprintf('surface slope          = %.12f\n',result.shootingSlope);
    fprintf('shooting iterations    = %d\n',result.iterations);
    fprintf('mean time (%d runs)    = %.6f s\n',nRepeat,meanTime);
    fprintf('time standard deviation= %.6f s\n',stdTime);
end
