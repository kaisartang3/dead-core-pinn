clear; clc; close all;

% =========================
% GLOBAL GRAPHICS SETTINGS (ENHANCED)
% =========================
set(groot,'defaultFigureColor','w')
set(groot,'defaultAxesColor','w')

set(groot,'defaultTextInterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')

set(groot,'defaultAxesFontName','Times')
set(groot,'defaultTextFontName','Times')

% ---- Larger fonts ----
fntsize=20;
set(groot,'defaultAxesFontSize',fntsize)
set(groot,'defaultTextFontSize',fntsize)

set(groot,'defaultAxesBox','on')
set(groot,'defaultAxesLineWidth',1.2)

set(groot,'defaultLineLineWidth',2)

% =========================
% Parameters
% =========================
n = 0.5;
mu_star = sqrt(2*(n+1)) / (1-n);

phi_values = [2, mu_star, 6];

fprintf('phi* = %.6f\n', mu_star);

% =========================
% Figure setup
% =========================
figure('Renderer','painters');
set(gcf,'Position',[100 100 1400 420]);

tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

for k = 1:3

    phi = phi_values(k);
    xdz = max(0, 1 - mu_star/phi); % ensure non-negative

    % =========================
    % Title
    % =========================
    if abs(phi - mu_star) < 1e-8
        title_str = sprintf('$\\phi = \\phi^* = %.3f$ (critical)', mu_star);
    elseif phi < mu_star
        title_str = sprintf('$\\phi = %.2f$ (subcritical)', phi);
    else
        title_str = sprintf('$\\phi = %.2f$ (supercritical)', phi);
    end

    % =========================
    % ODE (no artificial regularization)
    % =========================
    odefun = @(x,y) [y(2); phi^2 * (max(y(1),0)).^n];

    % =========================
    % Solve BVP
    % =========================
    if phi < mu_star
        % ---- Subcritical ----
        bcfun = @(ya,yb) [ya(1)-1; yb(1)-1];

        sol = bvp4c(odefun, bcfun, ...
            bvpinit(linspace(-1,1,100), @(x)[max(1-x.^2,0); -2*x]));

        x_full = linspace(-1,1,400);
        y = deval(sol,x_full);
        u_full = y(1,:);

    elseif abs(phi - mu_star) < 1e-8
        % ---- Critical ----
        bcfun = @(ya,yb) [ya(1)-1; yb(1)-1];

        sol = bvp4c(odefun, bcfun, ...
            bvpinit(linspace(-1,1,100), @(x)[max(1-x.^2,0); -2*x]));

        x_full = linspace(-1,1,400);
        y = deval(sol,x_full);
        u_full = y(1,:);

    else
        % ---- Supercritical ----
        sol = bvp4c(odefun, @(ya,yb)[ya(1); yb(1)-1], ...
            bvpinit(linspace(xdz,1,100), @(x)[(x-xdz)/(1-xdz);1]));

        x_pos = linspace(xdz,1,250);
        y_pos = deval(sol,x_pos);
        u_pos = y_pos(1,:);

        % ---- Full domain reconstruction (vectorized) ----
        x_full = linspace(-1,1,400);
        u_full = zeros(size(x_full));

        xi = abs(x_full);
        mask = xi > xdz;

        u_full(mask) = interp1(x_pos,u_pos,xi(mask),'pchip');
        u_full(~mask) = 0;
    end

    % =========================
    % Restrict to [0,1]
    % =========================
    mask_plot = x_full >= 0;
    x = x_full(mask_plot);
    u = u_full(mask_plot);

    % =========================
    % Plot
    % =========================
    nexttile;
    plot(x,u,'k'); hold on;

    % Dead zone
    if phi > mu_star
        xline(xdz,'k--','LineWidth',1.5);
        text(xdz/2,0.1, sprintf('$x_{dz}=%.3f$',xdz),'FontSize',fntsize-2);
    end

    xlabel('$x$');
    if k == 1
        ylabel('$u(x)$');
    end

    title(title_str);

    ofs=0.03;
    xlim([0-ofs 1+ofs]);
    ylim([0-ofs 1+ofs]);

    grid on;
end

% =========================
% EXPORT
% =========================
exportgraphics(gcf,'proceeding_comparison.pdf','ContentType','vector');

disp('Saved as proceeding_comparison.pdf');
