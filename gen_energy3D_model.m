clear;
base=[pwd(),filesep(),'..',filesep(),'DDE_Biftool_Nov2025',filesep()];
ddebiftool_path(base,'sym');
%%
format compact
format short g
%% ------------------------------------------------------------------------
%  Three-state energy-carbon substitution model (paper system (1))
%
%   Ndot = r*N*(1 - N/k) - beta*N*R/(1+h*N) - mu*N*C
%   Rdot = R*(-d + gamma*N/(1+h*N) + sigma*C) + aR*R*(1 - R/KR)
%   Cdot = ep*N - delta*C - phi*R*C
%
%  N : nonrenewable generation,  R : renewable generation,
%  C : accumulated-emissions / carbon-policy pressure.
%  The exogenous term rho*N has been absorbed into effective (r,k)
%  (see Remark 2.1), so it does not appear as a free parameter here.
% ------------------------------------------------------------------------
%% Parameters
%  N-equation:  r, k, beta, h, mu
%  R-equation:  d, gamma, sigma, aR, KR
%  C-equation:  ep (=epsilon, emission coeff.), delta, phi
%  tau : delay slot, currently unused (ntau=0); kept for the delayed
%        extension discussed in the paper.
%  NOTE: the emission coefficient epsilon is named 'ep' to avoid shadowing
%        MATLAB's built-in eps (machine epsilon).
parnames = {'r','k','beta','h','mu',...
            'd','gamma','sigma','aR','KR',...
            'ep','delta','phi',...
            'tau'};
cind = [parnames; num2cell(1:length(parnames))];
in   = struct(cind{:});
syms(parnames{:});
par = sym(parnames);
%% Symbolic state setup (3 states, single strain, no delays)
[dim,nstrains,ntau]=deal(3,1,1);
N = sym('N', [nstrains, ntau+1]);   % non-renewable energy
R = sym('R', [nstrains, ntau+1]);   % renewable energy
C = sym('C', [nstrains, ntau+1]);   % accumulated emissions / carbon policy
N1=N(1); R1=R(1); C1=C(1); N2=N(2);
%% Right-hand side
dNdt = r.*N1.*(1 - N1./k) ...
       - (beta.*N1.*R1)./(1 + h.*N1) ...
       - mu.*N1.*C1;

dRdt = R1.*(-d + (gamma.*N2)./(1 + h.*N2) + sigma.*C1) ...
       + aR.*R1.*(1 - R1./KR);

dCdt = ep.*N1- delta.*C1 - phi.*R1.*C1;
%% Generate the symbolic right-hand-side functions
%  maxorder=5 so that generalized-Hopf (Bautin) normal-form coefficients
%  are available; use maxorder=3 if only Hopf criticality (first Lyapunov
%  coefficient) and zero-Hopf are needed (smaller/faster generated file).
[fstr,derivs]=dde_sym2funcs([dNdt; dRdt; dCdt], [N; R; C], par, ...
    'filename','sym_energy_carbon_delay', ...
    'maxorder',5, ...
    'directional_derivative',true);