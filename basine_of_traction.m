%% BASIN OF ATTRACTION (2D) for an energy-transition ODE
% This script:
% 1) Sweeps a grid of initial conditions in (N,R) (or (x,y))
% 2) Integrates with ode45
% 3) Classifies the long-time outcome as:
%       - Equilibrium #1 / #2 / ...  (closest equilibrium)
%       - Limit cycle (if not close to any equilibrium, but oscillatory)
%       - Divergent / invalid (if it blows up or leaves domain)
% 4) Plots a basin map
%
% NOTE: You MUST implement your ODE in the function f_rhs below
% and list the equilibria you want to classify against.

clear; %clc; %close all;
% base=[pwd(),filesep(),'..',filesep(),'DDE_Biftool_Nov2025',filesep()];
% ddebiftool_path(base,'sym');
  base=[pwd(),filesep(),'..',filesep(),'DDE_Biftool_Nov2025',filesep()];
addpath([base,filesep,'ddebiftool'],...
    [base,filesep,'ddebiftool_extra_psol'],...
    [base,filesep,'ddebiftool_utilities'],...
    [base,filesep,'ddebiftool_extra_rotsym'],...
    [base,filesep,'ddebiftool_extra_nmfm'],...
    [base,filesep,'ddebiftool_extra_symbolic'],...
    [base,filesep,'ddebiftool_coco']);
format compact
format short 
load("steady_state_br.mat")
%% ------------------- USER SETTINGS -------------------
% Choose integration horizon (long enough to settle)
Tend   = 50;%300;       % increase if convergence is slow
Ttrans = 0.7*Tend;  % transient cutoff when detecting periodicity

% Grid of initial conditions
Nmin = 0;  Nmax = 1.2;   % adjust domain to your variables (scaled vs unscaled)
Rmin = 0;  Rmax = 1.2;
nN = 30;             % grid resolution
nR = 30;

% Classification tolerances
eqTol      = 1e-3;  % distance to equilibrium threshold
ampTol     = 1e-3;  % oscillation amplitude threshold to call it a cycle
divBound   = 50;    % if state norm exceeds this, call divergent

% ODE solver options (tighten if needed)
opts = odeset('RelTol',1e-4,'AbsTol',1e-4,'MaxStep',0.5);

%% ------------------- PARAMETERS -------------------
% Put your parameter values here.
% Example placeholders (edit to your model):
% par.sigma = 0.25;
% par.K     = 1.0;
% par.a     = 0.8;
% par.b     = 1.0;
% par.muN   = 0.2;
% par.muR   = 0.4;

% If you're using the "scaled version" of your other model:
% par.alpha = ...; par.gamma = ...; par.eps = ...; par.d = ...; par.kappa = ...;

%% ------------------- EQUILIBRIA LIST -------------------
% You need equilibria to classify basins.
% Provide them as columns: EQ = [N1 N2 ...; R1 R2 ...]
%
% If you don't know them analytically, you can:
% - compute them separately with fsolve, or
% - use MATCONT equilibrium points.
%
% Example: for your original (N+R)/K model in (x,y) you might have:
%   fossil-only: (x*,0), renewable-only: (0,y*), coexistence: (x~,y~)
%
% Put YOUR equilibria here:
[~,ind1]=min(abs(xpar_br1-2.9));
[~,ind2]=min(abs(xpar_br2-2.9));
pp2=eq_branch2.point(ind2);
pp3=eq_branch1.point(ind1);
pp2.parameter(in.sigma)=pp3.parameter(in.sigma);
%%
[p_br,suc0]=SetupStst(funcs,...
    'parameter',pp2.parameter,'x',pp2.x,...
    'contpar',in.sigma,'step',0.05,parbd{:});
%[pp2,success]=p_correc(funcs,pp2.x,[],[],method.point)
%%
EQ = [ ...
    p_br.point(1).x'% (stable)0.8,  0.0;  % eq1: (N*,R*)
    pp3.x'    %(unstable) eq2: (N*,R*)
]';

% If you have 3 equilibria, add a 3rd column, etc.

%% ------------------- GRID SETUP -------------------
Ngrid = linspace(Nmin,Nmax,nN);
Rgrid = linspace(Rmin,Rmax,nR);
basin = zeros(nR,nN);  % integer label for attractor
% Labels convention:
% 1..Neq  => equilibrium index
% Neq+1   => limit cycle
% Neq+2   => divergent/invalid

Neq = size(EQ,2);
LABEL_CYCLE = Neq + 1;
LABEL_DIV   = Neq + 2;
LABEL_UNKNOWN = Neq + 3;
%% ------------------- MAIN LOOP -------------------
fprintf('Computing basins on %dx%d grid...\n', nR, nN);

for i = 1:nN
    for j = 1:nR
        x0 = [Ngrid(i); Rgrid(j)];

        % Skip negative or invalid initial conditions if your model needs positivity
        if any(x0 < 0)
            basin(j,i) = LABDEL_DIV;
            continue;
        end

        try
            %[t,x] = ode45(@(t,x) f_rhs(t,x,par), [0 Tend], x0, opts);
            Fode45=@(t,X)funcs.wrap_rhs(X,eq_branch2.point(ind2).parameter);
            [t,x]=ode45(Fode45,[0,300],x0,opts);
%
        catch
            basin(j,i) = LABEL_DIV;
            continue;
        end

        if any(~isfinite(x(:))) || max(vecnorm(x,2,2)) > divBound
            basin(j,i) = LABEL_DIV;
            continue;
        end

        % Use tail segment for classification
        tailIdx = find(t >= Ttrans, 1, 'first');
        if isempty(tailIdx), tailIdx = floor(numel(t)/2); end
        xtail = x(tailIdx:end,:);

        xend = xtail(end,:)';

        % 1) Check if close to any equilibrium
        dists = vecnorm(xend - EQ, 2, 1);
        [dmin, kmin] = min(dists);

        if dmin < eqTol
            basin(j,i) = kmin; % converged to equilibrium kmin
            continue;
        end

        % 2) If not near equilibrium, test for sustained oscillation (limit cycle)
        % Simple robust test: amplitude over tail window
        ampN = max(xtail(:,1)) - min(xtail(:,1));
        ampR = max(xtail(:,2)) - min(xtail(:,2));

        % Also check whether it keeps moving (not settling)
        drift = norm(xtail(end,:)-xtail(1,:));

        if (max(ampN, ampR) > ampTol) && (drift < 0.1*max(ampN,ampR) + 1e-6)
            % Sustained oscillation → limit cycle
            basin(j,i) = LABEL_CYCLE;
        else
            % Neither equilibrium nor clear cycle → unknown / transient
            basin(j,i) = LABEL_UNKNOWN;
        end

    end
    if mod(i,10)==0
        fprintf('  %d / %d columns done\n', i, nN);
    end
end

%% ------------------- PLOT BASIN MAP -------------------
figure('Color','w');
imagesc(Ngrid, Rgrid, basin);
set(gca,'YDir','normal');
xlabel('N(0)'); ylabel('R(0)');
title(sprintf('Basin map (sigma=%.4g)', eq_branch2.point(ind2).parameter(in.sigma)));
colorbar;
colormap(parula(max(basin(:))));

% Make a legend-like text (basic)
hold on;
text(Nmin + 0.02*(Nmax-Nmin), Rmax - 0.05*(Rmax-Rmin), ...
    sprintf('1..%d = equilibria, %d = limit cycle, %d = divergent', ...
    Neq, LABEL_CYCLE, LABEL_DIV), 'Color','k','FontSize',10, 'FontWeight','bold');

%% ------------------- OPTIONAL: OVERLAY EQUILIBRIA -------------------
plot(EQ(1,:), EQ(2,:), 'ko', 'MarkerFaceColor','w','MarkerSize',7);
hold off;

%% ================================================================
%% RHS FUNCTION: PUT YOUR MODEL HERE
% function dx = f_rhs(~, x, par)
%     N = x(1);
%     R = x(2);
% 
%     % ------------------ CHOOSE ONE MODEL ------------------
%     % (A) Earlier model with shared carrying capacity (N+R)/K and substitution sigma
%     %     dR/dt = a R (1-(N+R)/K) + sigma N R - muR R
%     %     dN/dt = b N (1-(N+R)/K) - sigma N R - muN N
%     %
%     % EDIT par.* to match your notation:
%     a    = par.a;
%     b    = par.b;
%     muR  = par.muR;
%     muN  = par.muN;
%     K    = par.K;
%     sig  = par.sigma;
% 
%     dR = a*R*(1 - (N+R)/K) + sig*N*R - muR*R;
%     dN = b*N*(1 - (N+R)/K) - sig*N*R - muN*N;
% 
%     % (B) If instead you want the scaled "alpha,b,gamma,d,kappa,eps" model, comment (A)
%     % and uncomment this block:
%     % alpha = par.alpha;  b = par.b;
%     % gamma = par.gamma;  d = par.d;
%     % kappa = par.kappa;  eps = par.eps;
%     % dN = alpha - b*N - gamma*N*R;
%     % dR = d*R*(1 - R/kappa) + eps*N*R;
% 
%     dx = [dN; dR];
% end
