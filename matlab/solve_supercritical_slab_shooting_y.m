function result = solve_supercritical_slab_shooting_y(phi,n,opts)
%SOLVE_SUPERCRITICAL_SLAB_SHOOTING_Y
% Robust shooting formulation using
%
%       y = u^((1-n)/2).
%
% The original active problem is
%
%       u'' = phi^2 u^n,   u(x_dz)=u'(x_dz)=0,   u(1)=1.
%
% Set q=2/(1-n), so u=y^q. Then
%
%       y*y'' + (q-1)*(y')^2 = phi^2/q.
%
% We shoot on the unknown interface a=x_dz.
%
% Near the interface,
%
%       y(x) ~ C (x-a),
%
% where
%
%       C = phi/sqrt(q*(q-1)).
%
% Thus the numerical initialization is O(eps), rather than
% O(eps^q). This avoids the extremely small concentrations that
% occur in direct-u shooting when n is close to one.
%
% IMPORTANT:
% The analytical x_dz is NOT used by this solver.

if nargin < 3
    opts = struct();
end

if ~isfield(opts,'RelTol'),     opts.RelTol = 1e-10; end
if ~isfield(opts,'AbsTol'),     opts.AbsTol = 1e-12; end
if ~isfield(opts,'tolA'),       opts.tolA = 1e-10; end
if ~isfield(opts,'eps'),        opts.eps = 1e-5; end
if ~isfield(opts,'maxIter'),    opts.maxIter = 100; end
if ~isfield(opts,'maxBracket'), opts.maxBracket = 100; end

if n <= 0 || n >= 1
    error('This transformed solver requires 0 < n < 1.');
end

q = 2/(1-n);
C = phi/sqrt(q*(q-1));

%% =========================================================
% Residual for a trial interface a
%
% Start at a+eps with the leading-order interface behavior.
%% =========================================================

    function [F,sol] = residual(a)

        epsLocal = min(opts.eps,0.1*(1-a));
        x0 = a + epsLocal;

        % Leading-order expansion y ~ C(x-a).
        y0  = C*epsLocal;
        yp0 = C;

        odeOpts = odeset( ...
            'RelTol',opts.RelTol, ...
            'AbsTol',opts.AbsTol);

        [x,Y] = ode45(@odefun,[x0 1],[y0;yp0],odeOpts);

        F = Y(end,1)-1;

        sol.x = x;
        sol.y = Y(:,1);
        sol.yp = Y(:,2);
        sol.x0 = x0;
        sol.y0 = y0;
        sol.yp0 = yp0;
    end

%% =========================================================
% Bracket a=x_dz
%% =========================================================

aLow = 0;
[FLow,~] = residual(aLow);

aHigh = 1 - 10*opts.eps;
[FHigh,~] = residual(aHigh);

kBracket = 0;

while FLow*FHigh > 0 && kBracket < opts.maxBracket

    aHigh = 1 - 0.5*(1-aHigh);

    if 1-aHigh <= 10*opts.eps
        break
    end

    [FHigh,~] = residual(aHigh);
    kBracket = kBracket + 1;
end

if FLow*FHigh > 0
    error(['Unable to bracket x_dz. ', ...
        'F(0)=%.6e, F(aHigh)=%.6e. ', ...
        'Try changing eps.'],FLow,FHigh);
end

%% =========================================================
% Bisection
%% =========================================================

for iter = 1:opts.maxIter

    aMid = 0.5*(aLow+aHigh);
    [FMid,~] = residual(aMid);

    if FLow*FMid <= 0
        aHigh = aMid;
        FHigh = FMid;
    else
        aLow = aMid;
        FLow = FMid;
    end

    if abs(aHigh-aLow) <= opts.tolA
        break
    end
end

aStar = 0.5*(aLow+aHigh);
[Ffinal,sol] = residual(aStar);

% Recover u from y for output.
u = sol.y.^q;

result.phi = phi;
result.n = n;
result.q = q;
result.C = C;
result.xdz = aStar;
result.residual = Ffinal;
result.iterations = iter;
result.x = sol.x;
result.y = sol.y;
result.yp = sol.yp;
result.u = u;
result.x0 = sol.x0;
result.y0 = sol.y0;
result.yp0 = sol.yp0;
result.eps = opts.eps;
result.success = true;

%% =========================================================
% Transformed ODE
%
% y*y'' + (q-1)*(y')^2 = phi^2/q
%
% Hence
%
% y'' = [phi^2/q - (q-1)*(y')^2]/y.
%% =========================================================

    function dydx = odefun(~,z)

        y  = max(z(1),realmin);
        yp = z(2);

        ypp = (phi^2/q - (q-1)*yp^2)/y;

        dydx = [yp;ypp];
    end

end
