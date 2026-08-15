%% =========================================================
% RUN EVERYTHING:
% TRANSFORMED SHOOTING SOLVER, CONVERGENCE, TIMING, AND FIGURES
% =========================================================
clear
clc
close all

fprintf('\n');
fprintf('============================================================\n');
fprintf('  SUPERCRITICAL SLAB: COMPLETE SHOOTING STUDY\n');
fprintf('============================================================\n');

%% ---------------------------------------------------------
% Cases
%% ---------------------------------------------------------
cases = [ ...
    0.50,  6.0000;
    0.90, 29.2404];

%% ---------------------------------------------------------
% Numerical options
%% ---------------------------------------------------------
opts.RelTol     = 1e-10;
opts.AbsTol     = 1e-12;
opts.tolA       = 1e-10;
opts.eps        = 1e-4;
opts.maxIter    = 100;
opts.maxBracket = 100;

nRepeat = 20;

%% ---------------------------------------------------------
% Check that solver is available
%% ---------------------------------------------------------
if exist('solve_supercritical_slab_shooting_y','file') ~= 2
    error(['Cannot find solve_supercritical_slab_shooting_y.m. ', ...
           'Place this script in the same folder as the solver.']);
end

%% =========================================================
% 1. EPSILON CONVERGENCE
% =========================================================
fprintf('\n');
fprintf('1. EPSILON CONVERGENCE\n');
fprintf('============================================================\n');

epsValues = [1e-2 3e-3 1e-3 3e-4 1e-4 3e-5 1e-5 3e-6 1e-6];

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

%% =========================================================
% 2. FINAL SHOOTING RESULTS AND TIMING
% =========================================================
fprintf('\n\n');
fprintf('2. FINAL SHOOTING RESULTS AND TIMING\n');
fprintf('============================================================\n');
fprintf('epsilon = %.1e\n\n',opts.eps);

fprintf('n       phi       x_dz            error            time (s)\n');

finalResults = zeros(size(cases,1),5);

for k = 1:size(cases,1)

    n   = cases(k,1);
    phi = cases(k,2);

    phiStar = sqrt(2*(1+n))/(1-n);
    xdzRef  = 1 - phiStar/phi;

    % One warm-up call, excluded from timing.
    solve_supercritical_slab_shooting_y(phi,n,opts);

    times = zeros(nRepeat,1);

    for j = 1:nRepeat
        tic
        result = solve_supercritical_slab_shooting_y(phi,n,opts);
        times(j) = toc;
    end

    meanTime = mean(times);
    stdTime  = std(times);
    errXdz   = abs(result.xdz-xdzRef);

    finalResults(k,:) = [n,phi,result.xdz,errXdz,meanTime];

    fprintf('%.2f   %.4f   %.12f   %.6e   %.6f +/- %.6f\n', ...
        n,phi,result.xdz,errXdz,meanTime,stdTime);
end

%% =========================================================
% 3. SOLUTION FIGURES
% =========================================================
fprintf('\n\n');
fprintf('3. GENERATING SOLUTION FIGURES\n');
fprintf('============================================================\n');

% Reset epsilon to final value.
opts.eps = 1e-4;

for k = 1:size(cases,1)

    n   = cases(k,1);
    phi = cases(k,2);

    result = solve_supercritical_slab_shooting_y(phi,n,opts);

    phiStar = sqrt(2*(1+n))/(1-n);
    xdzRef  = 1 - phiStar/phi;

    xFine = linspace(0,1,2000);

    % Exact solution.
    uRef = zeros(size(xFine));
    active = xFine >= xdzRef;
    uRef(active) = ...
        ((xFine(active)-xdzRef)/(1-xdzRef)).^(2/(1-n));

    % Shooting solution including dead-core endpoint.
    xPlot = [result.xdz, result.x(:)'];
    uPlot = [0, result.u(:)'];

    [xPlot,idx] = unique(xPlot,'stable');
    uPlot = uPlot(idx);

    % -----------------------------------------------------
    % Solution figure
    % -----------------------------------------------------
    fig = figure('Color','w');
    hold on

    plot(xFine,uRef,'LineWidth',2);
    plot(xPlot,uPlot,'--','LineWidth',1.8);
    xline(xdzRef,':','LineWidth',1.2);

    xlabel('$x$','Interpreter','latex');
    ylabel('$u(x)$','Interpreter','latex');

    title(sprintf('$n=%.2f,\\ \\phi=%.4f$',n,phi), ...
        'Interpreter','latex');

    legend({'Exact solution','Shooting solution', ...
            '$x_{dz}^{\rm ref}$'}, ...
            'Interpreter','latex', ...
            'Location','southeast');

    xlim([0 1]);
    ylim([0 1.02]);
    grid on
    box on
    set(gca,'FontSize',11);

    filename = sprintf('shooting_solution_n%.2f_phi%.4f',n,phi);

    exportgraphics(fig,[filename '.png'],'Resolution',300);
    exportgraphics(fig,[filename '.pdf'],'ContentType','vector');

    close(fig);

    % -----------------------------------------------------
    % Pointwise error figure
    % -----------------------------------------------------
    uShootFine = zeros(size(xFine));
    activeShoot = xFine >= result.xdz;

    uShootFine(activeShoot) = interp1( ...
        xPlot,uPlot,xFine(activeShoot),'linear','extrap');

    absError = abs(uShootFine-uRef);

    fig2 = figure('Color','w');

    semilogy(xFine,max(absError,1e-16),'LineWidth',1.5);
    hold on
    xline(xdzRef,':','LineWidth',1.2);

    xlabel('$x$','Interpreter','latex');
    ylabel('$|u_{\rm shoot}-u_{\rm ref}|$','Interpreter','latex');

    title(sprintf('Absolute error: $n=%.2f$, $\\phi=%.4f$', ...
        n,phi),'Interpreter','latex');

    legend({'Absolute error','$x_{dz}^{\rm ref}$'}, ...
        'Interpreter','latex','Location','best');

    xlim([0 1]);
    grid on
    box on
    set(gca,'FontSize',11);

    filename2 = sprintf('shooting_error_n%.2f_phi%.4f',n,phi);

    exportgraphics(fig2,[filename2 '.png'],'Resolution',300);
    exportgraphics(fig2,[filename2 '.pdf'],'ContentType','vector');

    close(fig2);

    fprintf('Created %s.[png,pdf]\n',filename);
    fprintf('Created %s.[png,pdf]\n',filename2);
end

%% =========================================================
% 4. COMBINED TWO-PANEL SOLUTION FIGURE
% =========================================================
fprintf('\n');
fprintf('Creating combined solution figure...\n');

fig = figure('Color','w');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

for k = 1:size(cases,1)

    n   = cases(k,1);
    phi = cases(k,2);

    result = solve_supercritical_slab_shooting_y(phi,n,opts);

    phiStar = sqrt(2*(1+n))/(1-n);
    xdzRef  = 1 - phiStar/phi;

    xFine = linspace(0,1,2000);

    uRef = zeros(size(xFine));
    active = xFine >= xdzRef;
    uRef(active) = ...
        ((xFine(active)-xdzRef)/(1-xdzRef)).^(2/(1-n));

    xPlot = [result.xdz, result.x(:)'];
    uPlot = [0, result.u(:)'];

    [xPlot,idx] = unique(xPlot,'stable');
    uPlot = uPlot(idx);

    nexttile
    hold on

    plot(xFine,uRef,'LineWidth',2);
    plot(xPlot,uPlot,'--','LineWidth',1.8);
    xline(xdzRef,':','LineWidth',1.2);

    xlabel('$x$','Interpreter','latex');
    ylabel('$u(x)$','Interpreter','latex');

    title(sprintf('$n=%.2f,\\ \\phi=%.4f$',n,phi), ...
        'Interpreter','latex');

    legend({'Exact','Shooting','$x_{dz}$'}, ...
        'Interpreter','latex','Location','southeast');

    xlim([0 1]);
    ylim([0 1.02]);
    grid on
    box on
    set(gca,'FontSize',10);
end

exportgraphics(fig,'shooting_solution_comparison.pdf', ...
    'ContentType','vector');
exportgraphics(fig,'shooting_solution_comparison.png', ...
    'Resolution',300);

close(fig);

%% =========================================================
% 5. SAVE NUMERICAL RESULTS
% =========================================================
resultsTable = array2table(finalResults, ...
    'VariableNames',{'n','phi','xdz_shoot','interface_error','time_s'});

writetable(resultsTable,'shooting_final_results.csv');

save('shooting_final_results.mat','finalResults','opts','cases');

fprintf('\n');
fprintf('============================================================\n');
fprintf('COMPLETE SHOOTING STUDY FINISHED\n');
fprintf('============================================================\n');
fprintf('Generated:\n');
fprintf('  - individual solution PNG/PDF files\n');
fprintf('  - individual error PNG/PDF files\n');
fprintf('  - shooting_solution_comparison.pdf\n');
fprintf('  - shooting_solution_comparison.png\n');
fprintf('  - shooting_final_results.csv\n');
fprintf('  - shooting_final_results.mat\n');
fprintf('============================================================\n');
