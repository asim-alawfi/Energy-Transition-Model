clear;
base=[pwd(),filesep(),'..',filesep(),'DDE_Biftool_Nov2025',filesep()];
ddebiftool_path(base,'sym');
%   base=[pwd(),filesep(),'..',filesep(),'DDE_Biftool_Nov2025',filesep()];
% addpath([base,filesep,'ddebiftool'],...
%     [base,filesep,'ddebiftool_extra_psol'],...
%     [base,filesep,'ddebiftool_utilities'],...
%     [base,filesep,'ddebiftool_extra_rotsym'],...
%     [base,filesep,'ddebiftool_extra_nmfm'],...
%     [base,filesep,'ddebiftool_extra_symbolic'],...
%     [base,filesep,'ddebiftool_coco']);
format compact
format short 
%%
parnames = {'a','lambda','sigma','mu_n','mu_r','k','tau_delay'};
cind = [parnames; num2cell(1:length(parnames))];
in = struct(cind{:});
par([in.a,      in.lambda,     in.sigma,          in.mu_n,           in.mu_r,   in.k   in.tau_delay])=... 
    [0.2,         0.15,            2,              0.3,               0.6,            0.9       0 ];
ndim=2;
stst=dde_stst_create('parameter',par,...
    'x',[0;par(in.k)*(1-(par(in.mu_r)/par(in.a)))]);
%x0=[0;0.5]%zeros(ndim,1);
    parbd={'min_bound',[in.sigma, 0;  in.tau_delay,0.0; in.a,0],...
    'max_bound',   [in.sigma, 5;     in.tau_delay, 5; in.a,1],...
    'max_step',    [in.sigma, 0.1;   in.tau_delay, 0.1; in.a, 0.1; 0,0.1]};
%%
funcs=set_symfuncs(@sym_energy_E_version2,'sys_tau',@ in.tau_delay);
%%
[eq_branch1,suc0]=SetupStst(funcs,...
    'parameter',par,'x',stst.x,...
    'contpar',in.tau_delay,'step',0.1,parbd{:});
figure(1)
clf; hold on
eq_branch1=br_contn(funcs,eq_branch1,400);
%eq_1234_0=br_rvers(eq_1234_0);
%eq_1234_0=br_contn(funcs,eq_1234_0,400);
[eq_branch1,nunst_st]=br_stabl(funcs,eq_branch1,0,0);
%% plot stability
xpar_br1=arrayfun(@(x)x.parameter(in.tau_delay),eq_branch1.point);
ypar_br1=arrayfun(@(x)x.x(1),eq_branch1.point);
biflocation=find(diff(nunst_st));%find((nunst_st>0),1,'last');
%p1=eq_1234_0.point(biflocation);
%
% Extract stability information and plot results
figure(1);
clf
hold on 
plot(xpar_br1(nunst_st==0), ypar_br1(nunst_st==0), 'k.','Linewidth',2,'MarkerSize',10);
plot(xpar_br1(nunst_st==1), ypar_br1(nunst_st==1), 'r.','Linewidth',2,'MarkerSize',10);
plot(xpar_br1(nunst_st>=2), ypar_br1(nunst_st>=2), 'g.','Linewidth',2,'MarkerSize',10);
plot(xpar_br1(biflocation), ypar_br1(biflocation), 'r.','Linewidth',2,'MarkerSize',40);
xlabel('sigma');
ylabel('y_{1}');
title('Stability Analysis of Equilibria');
%%
pt=eq_branch1.point(biflocation+0);
xx=[pt.x]
% F_dde23=@(t,X,Ch)funcs.wrap_rhs(cat(2,X,Ch),pt.parameter)%,pt.parameter);
F_dde23=@(t,X,Ch)funcs.sys_rhs(cat(2,X,Ch),pt.parameter)%,pt.parameter);

rng(1)
simuldde23=dde23(F_dde23,pt.parameter(in.tau_delay),xx(:,1),[0,10])%,ddeset('RelTol',1e-6,'AbsTol',1e-6));
%

figure(88)
clf;
hold on
plot(simuldde23.x,simuldde23.y,'LineWidth',2)
%%

% [ffuncs,sts_fold,suc]=SetupFold(funcs,eq_1234_0,biflocation,...
%     'contpar',[in.sigma,in.a],...
%     'dir',in.sigma,'step',0.01,'print_residual_info',1,'outputfuncs',true)%'norm',false)
% figure(5)
% clf; hold on
% sts_fold=br_contn(ffuncs,sts_fold,80);

%%
  parbd={'min_bound',[in.sigma, 0;  in.mu_n,0.0; in.a,0],...
    'max_bound',   [in.sigma, 3;     in.mu_n, 2.5; in.a,1],...
    'max_step',    [in.sigma, 0.05;   in.mu_n, 0.02; in.a, 0.01; 0,0.05]};
[eq_branch2,suc0]=SetupStst(funcs,...
    'parameter',pt.parameter,'x',simuldde23.y(:,end),...
    'contpar',in.sigma,'step',0.05,parbd{:});
figure(2)
clf; hold on
eq_branch2=br_contn(funcs,eq_branch2,1000);
eq_branch2=br_rvers(eq_branch2);
eq_branch2=br_contn(funcs,eq_branch2,40);
[eq_branch2,nunst_st2]=br_stabl(funcs,eq_branch2,0,0);
%% plot stability
xpar_br2=arrayfun(@(x)x.parameter(in.sigma),eq_branch2.point);
ypar_br2=arrayfun(@(x)x.x(1),eq_branch2.point);
biflocation2=find(diff(nunst_st2));%find((nunst_st>0),1,'last');
%p1=eq_1234_0.point(biflocation);
%%
% Extract stability information and plot results
figure(1);
%clf
hold on 
plot(xpar_br2(nunst_st2==0), ypar_br2(nunst_st2==0), 'k.','Linewidth',2,'MarkerSize',10);
plot(xpar_br2(nunst_st2==1), ypar_br2(nunst_st2==1), 'r.','Linewidth',2,'MarkerSize',10);
plot(xpar_br2(biflocation2(1)), ypar_br2(biflocation2(1)), 'b.','Linewidth',2,'MarkerSize',30);
plot(xpar_br2(biflocation2(2)), ypar_br2(biflocation2(2)), 'g.','Linewidth',2,'MarkerSize',30);
%%%
plot(xpar_br1(nunst_st==0), ypar_br1(nunst_st==0), 'k.','Linewidth',2,'MarkerSize',10);
plot(xpar_br1(nunst_st==1), ypar_br1(nunst_st==1), 'r.','Linewidth',2,'MarkerSize',10);
plot(xpar_br1(biflocation), ypar_br1(biflocation), 'r.','Linewidth',2,'MarkerSize',20);
xlabel('sigma');
ylabel('y_{1}');
title('Stability Analysis of Equilibria');
%%

%plot(xpar_1234(biflocation), ypar_1234(biflocation), 'r.','Linewidth',2,'MarkerSize',20);
%%
%save("steady_state_br.mat")
%%
tspan = [0 30000];
%pt2=eq_branch1.point(biflocation2(2));
F_dde23=@(t,X,Ch)funcs.wrap_rhs(cat(2,X,Ch),pt.parameter);
pt=eq_branch1.point(biflocation(1)+3);
inc=[pt.x]
lags = par(in.tau_delay);
rng(2)
history = @(t) inc(:,1)+0.01*rand(length(pt.x),1);
simuldde23=dde23(F_dde23,lags,history,tspan);%,odeset('RelTol',1e-7,'AbsTol',1e-7));
%
figure(88)
clf;
hold on
plot(simuldde23.y(1,:),simuldde23.y(2,:),'LineWidth',2)