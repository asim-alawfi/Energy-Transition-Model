clear;
base=[pwd(),filesep(),'..',filesep(),'DDE_Biftool_Nov2025',filesep()];
ddebiftool_path(base,'sym');%[getenv('HOME'),'/sourceforge-ddebiftool/releases/git-2023-06-27'],'sym');
%%
format compact
format short g
%%
%% Define parameters
parnames = {'alpha','b','gamma','d','ku','beta','lambda','tau_delay'};
cind = [parnames; num2cell(1:length(parnames))];
in = struct(cind{:});
syms(parnames{:});
par = sym(parnames);
%% Symbolic setup
[dim,nstrains,ntau]=deal(2,1,0);
N = sym('N', [nstrains, ntau+1]);  % non-renewable enrgy
R = sym('R', [nstrains, ntau+1]);  % renewable enrgy
%% 
dNdt= alpha - b.*N-lambda.*N.^2- gamma.*N.*R;%Lambda.*((1-w)+(1-gamma).*(w-R))-N; %% need to check the term gamma*w or just gamma
dRdt= d*R.*(1 - R./ku)+ gamma.*N.*R;%Lambda.*(1-R)-e*R;
%%
[fstr,erives]=dde_sym2funcs([dNdt;dRdt],[N;R],par,'filename','sym_energy',...
    'maxorder',2,'directional_derivative',true);
%%


