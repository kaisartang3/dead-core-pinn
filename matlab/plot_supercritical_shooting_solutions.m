%% =========================================================
% SUPERCRITICAL SLAB:
% TRANSFORMED SHOOTING + SOLUTION PLOTS
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

for k = 1:size(cases,1)

    n   = cases(k,1);
    phi = cases(k,2);

    % -----------------------------------------------------
    % Shooting solution
    % -----------------------------------------------------
    result = solve_supercritical_slab_shooting_y(phi,n,opts);

    xdzShoot = result.xdz;

    % -----------------------------------------------------
    % Exact reference solution
    % -----------------------------------------------------
    phiStar = sqrt(2*(1+n))/(1-n);
    xdzRef  = 1 - phiStar/phi;

    % Fine grid for plotting.
    xFine = linspace(0,1,2000);

    uRef = zeros(size(xFine));
    active = xFine >= xdzRef;
    uRef(active) = ...
        ((xFine(active)-xdzRef)/(1-xdzRef)).^(2/(1-n));

    % Reconstruct the numerical shooting solution.
    % result.x contains only the active interval.
    xShoot = result.x;
    uShoot = result.u;

    % The numerical solution starts at x_dz+epsilon.
    % Add the asymptotic dead-core portion for plotting.
    xPlot = [xdzShoot, xShoot(:)'];
    uPlot = [0, uShoot(:)'];

    % Remove possible duplicate first point.
    [xPlot,idx] = unique(xPlot,'stable');
    uPlot = uPlot(idx);

    % -----------------------------------------------------
    % Plot solution
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

    grid on
    box on
    xlim([0 1]);
    ylim([0 1.02]);

    set(gca,'FontSize',11);

    % Save PNG and PDF.
    filename = sprintf('shooting_solution_n%.2f_phi%.4f',n,phi);
    exportgraphics(fig,[filename '.png'],'Resolution',300);
    exportgraphics(fig,[filename '.pdf'],'ContentType','vector');

    % -----------------------------------------------------
    % Error plot
    % -----------------------------------------------------
    % Interpolate the shooting solution onto the fine grid.
    uShootFine = zeros(size(xFine));
    activeShoot = xFine >= xdzShoot;

    uShootFine(activeShoot) = interp1( ...
        xPlot,uPlot,xFine(activeShoot),'linear','extrap');

    absError = abs(uShootFine-uRef);

    fig2 = figure('Color','w');
    semilogy(xFine,max(absError,1e-16),'LineWidth',1.5);
    hold on
    xline(xdzRef,':','LineWidth',1.2);

    xlabel('$x$','Interpreter','latex');
    ylabel('$|u_{\rm shoot}-u_{\rm ref}|$','Interpreter','latex');

    title(sprintf('Absolute error: $n=%.2f$, $\\phi=%.4f$',n,phi), ...
        'Interpreter','latex');

    legend({'Absolute error','$x_{dz}^{\rm ref}$'}, ...
        'Interpreter','latex', ...
        'Location','best');

    grid on
    box on
    xlim([0 1]);

    set(gca,'FontSize',11);

    filename2 = sprintf('shooting_error_n%.2f_phi%.4f',n,phi);
    exportgraphics(fig2,[filename2 '.png'],'Resolution',300);
    exportgraphics(fig2,[filename2 '.pdf'],'ContentType','vector');

    % -----------------------------------------------------
    % Console output
    % -----------------------------------------------------
    fprintf('\nCase n = %.2f, phi = %.4f\n',n,phi);
    fprintf('Reference x_dz = %.12f\n',xdzRef);
    fprintf('Shooting  x_dz = %.12f\n',xdzShoot);
    fprintf('|Delta x_dz|    = %.6e\n',abs(xdzShoot-xdzRef));
    fprintf('Figures saved: %s.[png,pdf]\n',filename);
    fprintf('Error figures saved: %s.[png,pdf]\n',filename2);
end
