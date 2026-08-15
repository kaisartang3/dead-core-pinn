function result = solve_supercritical_slab_shooting(phi,n,opts)
%SOLVE_SUPERCRITICAL_SLAB_SHOOTING
% Free-boundary shooting solver for the supercritical slab problem
%
%   u'' = phi^2 u^n,       x in (0,1),
%   u(x)=0 on [0,x_dz],
%   u(1)=1.
%
% The unknown shooting parameter is the dead-core interface a=x_dz.
% The active solution is initialized at x=a+eps using the exact local
% leading-order asymptotics
%
%   u ~ A (x-a)^q,   q=2/(1-n),
%
% where
%
%   A = [phi^2/(q(q-1))]^(1/(1-n)).
%
% For a trial a, the active solution is integrated from a+eps to x=1.
% The shooting residual is F(a)=u(1)-1.  Bisection is used to determine
% a such that F(a)=0.
%
% The analytical formula for x_dz is NOT used in determining a.
% It is used only by the comparison driver for the final error.

if nargin < 3
    opts = struct();
end

if ~isfield(opts,'RelTol'),       opts.RelTol = 1e-10; end
if ~isfield(opts,'AbsTol'),       opts.AbsTol = 1e-12; end
if ~isfield(opts,'tolA'),         opts.tolA = 1e-11; end
if ~isfield(opts,'eps'),          opts.eps = 1e-6; end
if ~isfield(opts,'maxIter'),      opts.maxIter = 100; end
if ~isfield(opts,'maxBracket'),   opts.maxBracket = 100; end

if n <= 0 || n >= 1
    error('This solver requires 0 < n < 1.');
end

q = 2/(1-n);

% Leading-order coefficient from
% A*q*(q-1) = phi^2*A^n.
A = (phi^2/(q*(q-1)))^(1/(1-n));

odeOpts = odeset( ...
    'RelTol',opts.RelTol, ...
    'AbsTol',opts.AbsTol);

%% =========================================================
%  Residual for a trial dead-core interface a
%% =========================================================

    function [F,sol] = residual(a)

        epsLocal = min(opts.eps,0.1*(1-a));

        x0 = a + epsLocal;
        z  = epsLocal;

        u0  = A*z^q;
        up0 = A*q*z^(q-1);

        [x,Y] = ode45(@odefun,[x0 1],[u0;up0],odeOpts);

        F = Y(end,1)-1;

        sol.x = x;
        sol.u = Y(:,1);
        sol.up = Y(:,2);
        sol.x0 = x0;
        sol.u0 = u0;
        sol.up0 = up0;
    end

%% =========================================================
%  Bracket the interface
%
%  If a is too small, the active region is too long and u(1)>1.
%  If a is too large, the active region is too short and u(1)<1.
%% =========================================================

aLow = 0;
[FLow,~] = residual(aLow);

aHigh = 1 - 10*opts.eps;
[FHigh,~] = residual(aHigh);

kBracket = 0;

while FLow*FHigh > 0 && kBracket < opts.maxBracket

    % Move the upper interface closer to 1.
    aHigh = 1 - 0.5*(1-aHigh);

    if 1-aHigh <= 10*opts.eps
        break
    end

    [FHigh,~] = residual(aHigh);
    kBracket = kBracket + 1;
end

if FLow*FHigh > 0
    error(['Unable to bracket x_dz. ', ...
           'Check eps and the selected supercritical case. ', ...
           'F(0)=%.6e, F(aHigh)=%.6e.'],FLow,FHigh);
end

%% =========================================================
%  Bisection in a=x_dz
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

result.phi = phi;
result.n = n;
result.xdz = aStar;
result.residual = Ffinal;
result.iterations = iter;
result.x = sol.x;
result.u = sol.u;
result.up = sol.up;
result.x0 = sol.x0;
result.u0 = sol.u0;
result.up0 = sol.up0;
result.q = q;
result.A = A;
result.eps = opts.eps;
result.success = true;

%% =========================================================
% Nested ODE
%% =========================================================

    function dydx = odefun(~,y)
        uSafe = max(y(1),0);
        dydx = [y(2);
                phi^2*uSafe^n];
    end

end
