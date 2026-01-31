clear;
base=[pwd(),filesep(),'..',filesep(),'DDE_Biftool_Nov2025',filesep()];
ddebiftool_path(base,'sym');%[getenv('HOME'),'/sourceforge-ddebiftool/releases/git-2023-06-27'],'sym');
%%
format compact
format short g
%%
%% Define parameters
parnames = {'a','b','sigma' ,'mu_n','mu_r','tau_delay'};
cind = [parnames; num2cell(1:length(parnames))];
in = struct(cind{:});
syms(parnames{:});
par = sym(parnames);
%% Symbolic setup
[dim,nstrains,ntau]=deal(2,1,0);
N = sym('N', [nstrains, ntau+1]);  % non-renewable enrgy
R = sym('R', [nstrains, ntau+1]);  % renewable enrgy
%% 
dNdt= N.*( b- b.*N -b.*R - mu_n) - sigma.*N.*R;%Lambda.*((1-w)+(1-gamma).*(w-R))-N; %% need to check the term gamma*w or just gamma
dRdt= R.*(a- a.*N -a.*R - mu_r)  + sigma.*N.*R;%Lambda.*(1-R)-e*R;
%%
[fstr,erives]=dde_sym2funcs([dNdt;dRdt],[N;R],par,'filename','sym_energy_E',...
    'maxorder',2,'directional_derivative',true);
%%


