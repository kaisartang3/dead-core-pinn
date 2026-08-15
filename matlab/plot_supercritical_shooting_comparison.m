%% =========================================================
% COMBINED FIGURE:
% SHOOTING VS EXACT SOLUTIONS FOR BOTH CASES
% =========================================================
clear
clc
close all

cases = [ ...
    0.50,  6.0000;
    0.90, 29.2404];

opts.RelTol     = 1e-10;
opts.AbsTol     = 1e-12;
opts.tolA       = 1e-10;
opts.eps        = 1e-4;
opts.maxIter    = 100;
opts.maxBracket = 100;

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

    xPlot = [xdzRef, result.x(:)'];
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

exportgraphics(fig,'shooting_solution_comparison.pdf','ContentType','vector');
exportgraphics(fig,'shooting_solution_comparison.png','Resolution',300);

fprintf('\nSaved combined figures:\n');
fprintf('  shooting_solution_comparison.pdf\n');
fprintf('  shooting_solution_comparison.png\n');
