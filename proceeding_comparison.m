clear; clc; close all;

% =========================
% GLOBAL GRAPHICS SETTINGS (CRITICAL)
% =========================
set(groot,'defaultFigureColor','w')
set(groot,'defaultAxesColor','w')

set(groot,'defaultTextInterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')

set(groot,'defaultAxesFontName','Times')
set(groot,'defaultTextFontName','Times')
set(groot,'defaultAxesFontSize',11)

set(groot,'defaultAxesBox','on')
set(groot,'defaultAxesLineWidth',0.8)

set(groot,'defaultLineLineWidth',1.5)

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
set(gcf,'Position',[100 100 1200 350]);

tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

for k = 1:3

    phi = phi_values(k);

    xdz = 1 - mu_star/phi;

    % =========================
    % Title (correct LaTeX)
    % =========================
    if abs(phi - mu_star) < 1e-6
        title_str = sprintf('$\\phi = \\phi^* = %.3f$ (critical)', mu_star);
    elseif phi < mu_star
        title_str = sprintf('$\\phi = %.0f$ (subcritical)', phi);
    else
        title_str = sprintf('$\\phi = %.0f$ (supercritical)', phi);
    end

    % =========================
    % ODE
    % =========================
    odefun = @(x,y) [y(2); phi^2 * max(y(1),1e-12)^n];

    if phi <= mu_star
        % ---- Subcritical / Critical ----
        bcfun = @(ya,yb) [ya(1)-1; yb(1)-1];

        sol = bvp4c(odefun, bcfun, ...
            bvpinit(linspace(-1,1,100), @(x)[1-(1-x.^2).^2;0]));

        x_full = linspace(-1,1,400);
        y = deval(sol,x_full);
        u_full = y(1,:);
        %xdz = 0;

    else
        % ---- Supercritical ----
        
        sol = bvp4c(odefun, @(ya,yb)[ya(1); yb(1)-1], ...
            bvpinit(linspace(xdz,1,100), @(x)[(x-xdz)/(1-xdz);1]));

        x_pos = linspace(xdz,1,200);
        y_pos = deval(sol,x_pos);
        u_pos = y_pos(1,:);

        x_full = linspace(-1,1,400);
        u_full = zeros(size(x_full));

        for i = 1:length(x_full)
            xi = abs(x_full(i));
            if xi <= xdz
                u_full(i) = 0;
            else
                u_full(i) = interp1(x_pos,u_pos,xi,'pchip');
            end
        end
    end

    % =========================
    % Restrict to [0,1]
    % =========================
    mask = x_full >= 0;
    x = x_full(mask);
    u = u_full(mask);

    % =========================
    % Plot
    % =========================
    nexttile;
    plot(x,u,'k'); hold on;

    % Dead zone annotation (FIXED LATEX)
    if xdz >= 0
        xline(xdz,'k--','LineWidth',1);
        text(0.1,0.2,sprintf('$x_{dz}=%.3f$',xdz));
    end

    xlabel('$x$');
    ylabel('$u(x)$');
    title(title_str);

    xlim([-0.1 1.1]);
    ylim([-0.1 1.1]);

    grid on;
end

% =========================
% EXPORT (MODERN + CORRECT)
% =========================
exportgraphics(gcf,'proceeding_comparison.pdf','ContentType','vector');

disp('Saved as proceeding_comparison.pdf');