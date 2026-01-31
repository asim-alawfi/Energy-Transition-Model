%% BASIN MAP: "Converges to an equilibrium" vs "Does NOT converge"
% Classification rule:
%   - If the final state is within eqTol of ANY equilibrium in EQ  -> label = 1 (equilibrium basin)
%   - Else                                                       -> label = 0 (non-equilibrium: cycle/transient/escape)
%
% You must provide:
%   1) The ODE RHS in f_rhs()
%   2) The equilibrium list EQ (columns are equilibria)
%
% This is intentionally minimal and matches your criterion exactly.

%%clear; clc; close all;

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
% Tend   = 50;%300;       % increase if convergence is slow
% Ttrans = 0.7*Tend;  % transient cutoff when detecting periodicity
% 
% % Grid of initial conditions
% Nmin = 0;  Nmax = 1.2;   % adjust domain to your variables (scaled vs unscaled)
% Rmin = 0;  Rmax = 1.2;
% nN = 30;             % grid resolution
% nR = 30;
% 
% % Classification tolerances
% eqTol      = 1e-3;  % distance to equilibrium threshold
% ampTol     = 1e-3;  % oscillation amplitude threshold to call it a cycle
% divBound   = 50;    % if state norm exceeds this, call divergent
% 
% % ODE solver options (tighten if needed)
% opts = odeset('RelTol',1e-4,'AbsTol',1e-4,'MaxStep',0.5);

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
%% ------------------- USER SETTINGS -------------------
Tend   = 100;          % integration time (increase if slow convergence)
Ttrans = 0.7*Tend;     % not used here, but kept if you later add extra checks

% Grid of initial conditions (scaled variables)
Nmin = 0; Nmax = 1.2;
Rmin = 0; Rmax = 1.2;
nN = 50;
nR = 50;

eqTol    = 1e-4;       % distance-to-equilibrium threshold
divBound = 50;         % if ||x|| exceeds this -> treat as "not converged"

opts = odeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',0.5);

%% ------------------- EQUILIBRIA LIST -------------------
% Columns are equilibria: EQ = [N1 N2 ...; R1 R2 ...]
% Replace with equilibria at your chosen sigma (from MATCONT or analytic formulas).
EQ = [ ...
    p_br.point(1).x'% (stable)0.8,  0.0;  % eq1: (N*,R*)
    pp3.x'    %(unstable) eq2: (N*,R*)
]';

Neq = size(EQ,2);

%% ------------------- GRID SETUP -------------------
Ngrid = linspace(Nmin,Nmax,nN);
Rgrid = linspace(Rmin,Rmax,nR);

% basin(j,i) = 1 means converged to an equilibrium; 0 means not
basin = zeros(nR,nN,'uint8');

%fprintf('Computing basin map %dx%d (sigma=%.4g)...\n', nR, nN, par.sigma);

for i = 1:nN
    for j = 1:nR
        x0 = [Ngrid(i); Rgrid(j)];

        % optional: enforce nonnegativity
        if any(x0 < 0)
            basin(j,i) = 0;
            continue;
        end

        try
             Fode45=@(t,X)funcs.wrap_rhs(X,eq_branch2.point(ind2).parameter);
            [~,x]=ode45(Fode45,[0,300],x0,opts);
        catch
            basin(j,i) = 0;
            continue;
        end

        if any(~isfinite(x(:))) || max(vecnorm(x,2,2)) > divBound
            basin(j,i) = 0;
            continue;
        end

        xend = x(end,:)';

        % Distance to nearest equilibrium
        dmin = min(vecnorm(xend - EQ, 2, 1));

        % Classification: converged to equilibrium or not
        if dmin < eqTol
            basin(j,i) = 1;
        else
            basin(j,i) = 0;
        end
    end
    if mod(i,10)==0
        fprintf('  %d/%d columns done\n', i, nN);
    end
end

%% ------------------- PLOT -------------------
figure('Color','w');
imagesc(Ngrid, Rgrid, double(basin));
set(gca,'YDir','normal');
xlabel('N(0)'); ylabel('R(0)');
title(sprintf('Basin: equilibrium convergence (1=yes, 0=no),  \\sigma=%.4g', par.sigma));

% Two-color colormap: 0=white, 1=black (you can invert if you like)
colormap([1 1 1; 0 0 0]);
c = colorbar;
c.Ticks = [0 1];
c.TickLabels = {'not equilibrium','equilibrium'};

hold on;
plot(EQ(1,:), EQ(2,:), 'ko', 'MarkerFaceColor','y', 'MarkerSize',7);
hold off;

%% ================================================================
