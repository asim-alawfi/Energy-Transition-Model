clear;
base=[pwd(),filesep(),'..',filesep(),'DDE_Biftool_Nov2025',filesep()];
ddebiftool_path(base,'sym');
%%
format compact
format short g
%% ------------------------------------------------------------------------
%  Three-state energy-carbon substitution model 

% ------------------------------------------------------------------------
%% Parameters

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
%% Generate the symbolic right-hand-side functions.
[fstr,derivs]=dde_sym2funcs([dNdt; dRdt; dCdt], [N; R; C], par, ...
    'filename','sym_energy_carbon_delay', ...
    'maxorder',5, 'directional_derivative',true);