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
parnames = {'a','r','gamma','d','k','beta','mu_n','h','tau_delay'};
cind = [parnames; num2cell(1:length(parnames))];
in = struct(cind{:});
par([in.a,      in.r,        in.gamma,      in.d,        in.k,      in.beta,   in.mu_n    in.h,   in.tau_delay])=... 
    [0.01,        1,            1.2,        0.6,        0.3,           1,      0.01,       0.5,       0 ];
ndim=2;
stst=dde_stst_create('parameter',par,...
    'x',[0;0])%par(in.k)*(1-(par(in.mu_n)/par(in.a)))]);
%x0=[0;0.5]%zeros(ndim,1);
    parbd={'min_bound',[in.gamma, 0;  in.k,0.0; in.a,0],...
    'max_bound',   [in.gamma, 5;     in.k, 5; in.a,1],...
    'max_step',    [in.gamma, 0.1;   in.k, 0.05; in.a, 0.1; 0,0.05]};
%%
funcs=set_symfuncs(@sym_energy_hopf,'sys_tau',@ () []);
%%
[eq_branch1,suc1]=SetupStst(funcs,...
    'parameter',par,'x',stst.x,...
    'contpar',in.k,'step',0.05,parbd{:});
figure(1)
clf; hold on
eq_branch1=br_contn(funcs,eq_branch1,400);
eq_branch1=br_rvers(eq_branch1);
eq_branch1=br_contn(funcs,eq_branch1,400);
[eq_branch1,nunst_st]=br_stabl(funcs,eq_branch1,0,0);
%% plot stability
xpar_br1=arrayfun(@(x)x.parameter(in.k),eq_branch1.point);
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
pt=eq_branch1.point(50);
Fode45=@(t,X)funcs.wrap_rhs(X,pt.parameter);
rng(1)
simulode45=ode45(Fode45,[0,20000],pt.x+0.1*rand(length(pt.x),1),odeset('RelTol',1e-6,'AbsTol',1e-6));
%

figure(88)
clf;
hold on
plot(simulode45.x,simulode45.y,'LineWidth',2)
% pt=eq_branch1.point(50+0);
% xx=[pt.x]
% % F_dde23=@(t,X,Ch)funcs.wrap_rhs(cat(2,X,Ch),pt.parameter)%,pt.parameter);
% F_dde23=@(t,X,Ch)funcs.sys_rhs(cat(2,X,Ch),pt.parameter)%,pt.parameter);
% 
% rng(1)
% simuldde23=dde23(F_dde23,pt.parameter(in.tau_delay),xx(:,1),[0,10])%,ddeset('RelTol',1e-6,'AbsTol',1e-6));
% %
% 
% figure(88)
% clf;
% hold on
% plot(simuldde23.x,simuldde23.y,'LineWidth',2)
%%

% [ffuncs,sts_fold,suc]=SetupFold(funcs,eq_1234_0,biflocation,...
%     'contpar',[in.sigma,in.a],...
%     'dir',in.sigma,'step',0.01,'print_residual_info',1,'outputfuncs',true)%'norm',false)
% figure(5)
% clf; hold on
% sts_fold=br_contn(ffuncs,sts_fold,80);

%%
  % parbd={'min_bound',[in.sigma, 0;  in.mu_n,0.0; in.a,0],...
  %   'max_bound',   [in.sigma, 3;     in.mu_n, 2.5; in.a,1],...
  %   'max_step',    [in.sigma, 0.05;   in.mu_n, 0.02; in.a, 0.01; 0,0.05]};
[eq_branch2,suc0]=SetupStst(funcs,...
    'parameter',pt.parameter,'x',simulode45.y(:,end),...
    'contpar',in.k,'step',0.05,parbd{:});
figure(1)
 hold on
eq_branch2=br_contn(funcs,eq_branch2,1000);
eq_branch2=br_rvers(eq_branch2);
%%
eq_branch2=br_contn(funcs,eq_branch2,120);
[eq_branch2,nunst_st2]=br_stabl(funcs,eq_branch2,0,0);
%% plot stability
[eq_branch2,bif1testfuncs,biflocation2,bif1types]=LocateSpecialPoints(funcs,eq_branch2,'debug',true);
xpar_br2=arrayfun(@(x)x.parameter(in.k),eq_branch2.point);
ypar_br2=arrayfun(@(x)x.x(1),eq_branch2.point);
%biflocation2=find(diff(nunst_st2));%find((nunst_st>0),1,'last');
%p1=eq_1234_0.point(biflocation);
[eq_branch2,nunst_st2]=br_stabl(funcs,eq_branch2,0,0);
%%
biflocation2=find(diff(nunst_st2));

%%
% Extract stability information and plot results
figure(1);
clf
hold on 
plot(xpar_br2(nunst_st2==0), ypar_br2(nunst_st2==0), 'k.','Linewidth',2,'MarkerSize',10);
plot(xpar_br2(nunst_st2>=1), ypar_br2(nunst_st2>=1), 'r.','Linewidth',2,'MarkerSize',10);
plot(xpar_br2(biflocation2(1)), ypar_br2(biflocation2(1)), 'b.','Linewidth',2,'MarkerSize',30);
plot(xpar_br2(biflocation2(2)), ypar_br2(biflocation2(2)), 'g.','Linewidth',2,'MarkerSize',30);
%
% plot(xpar_br1(nunst_st==0), ypar_br1(nunst_st==0), 'k.','Linewidth',2,'MarkerSize',10);
% plot(xpar_br1(nunst_st==1), ypar_br1(nunst_st==1), 'r.','Linewidth',2,'MarkerSize',10);
% plot(xpar_br1(biflocation), ypar_br1(biflocation), 'r.','Linewidth',2,'MarkerSize',20);
xlabel('k');
ylabel('y_{1}');
title('Stability Analysis of Equilibria');
%%
[psol_branch,suc]=SetupPsol(funcs,eq_branch2,biflocation2(1),'degree',6,'intervals',120,...
    'max_step',[in.k,0.1],'max_bound',[in.k,9]);
figure(19);
hold on%clf;ax2=gca;
[psol_branch,s,f,r]=br_contn(funcs,psol_branch,70);
%%

%%
% look at the period along the branch:
figure(4); clf;
Plot2dBranch(psol_branch,'y',@(p)p.period,'funcs',funcs);
xlabel('a21');
ylabel('period');

%plot(xpar_1234(biflocation), ypar_1234(biflocation), 'r.','Linewidth',2,'MarkerSize',20);
%%
%save("steady_state_br.mat")
%%
pt2=eq_branch2.point(end);
Fode45=@(t,X)funcs.wrap_rhs(X,pt2.parameter);
rng(1)
simulode45=ode45(Fode45,[0,2000],pt2.x+0.01*rand(length(pt2.x),1),odeset('RelTol',1e-6,'AbsTol',1e-6));
%

figure(88)
clf;
hold on
plot(simulode45.x,simulode45.y,'LineWidth',2)
%%
[eq_branch3,suc0]=SetupStst(funcs,...
    'parameter',pt.parameter,'x',simulode45.y(:,end),...
    'contpar',in.k,'step',0.1,parbd{:});
figure(1)
 hold on
eq_branch3=br_contn(funcs,eq_branch3,1000);
eq_branch3=br_rvers(eq_branch3);
eq_branch3=br_contn(funcs,eq_branch3,400);
[eq_branch3,nunst_st3]=br_stabl(funcs,eq_branch3,0,0);

%%
  parbd2={'min_bound',[in.gamma, 0;  in.k,0.0; in.d,0],...
    'max_bound',   [in.gamma, 5;     in.k, 5; in.d,1],...
    'max_step',    [in.gamma, 0.1;   in.k, 0.05; in.d, 0.02; 0,0.05]};
%indhopf2=indbifhopf(find(strcmp('hoho',bifhopftypes),1,'first'));
[fhopf,hopf2_branch,suc]=SetupHopf(funcs,eq_branch2,biflocation2(1),'contpar',[in.k,in.d],...
    'dir',in.k,'step',0.05,'outputfuncs',true,'max_bound',[in.k,9],'use_tangent',true);
%hopf2_branch.method.point.print_residual_info=0;
figure(5)
clf;
hold on
[hopf2_branch,s,f,r]=br_contn(fhopf,hopf2_branch,70);
hopf2_branch=br_rvers(hopf2_branch);
[hopf2_branch,s,f,r]=br_contn(fhopf,hopf2_branch,70);
xlabel('a21');ylabel('\tau_s');
[hopf2_branch,hf_sb3tests,hf_sb3bifs,hf_sb3bifind]=MonitorChange(fhopf,hopf2_branch,...
    'printlevel',1,'min_iterations',5,'exclude_trivial',true)

%%
[ffold,fold2_branch,suc]=SetupFold(funcs,eq_branch2,biflocation2(2),'contpar',[in.k,in.d],...
    'dir',in.k,'step',0.05,'outputfuncs',true,'max_bound',[in.k,9])%,'use_tangent',true);
%hopf2_branch.method.point.print_residual_info=0;
figure(5)
[fold2_branch,s,f,r]=br_contn(ffold,fold2_branch,100);
fold2_branch=br_rvers(fold2_branch);
[fold2_branch,s,f,r]=br_contn(ffold,fold2_branch,100);
xlabel('a21');ylabel('\tau_s');
%%
[fold2_branch,ff_tests,ff_bifs,ff_bifind]=MonitorChange(ffold,fold2_branch,...
    'printlevel',1,'min_iterations',5,'exclude_trivial',true)
%%
%[hopf2_branch,nunst_hf]=br_stabl(fhopf,hopf2_branch,0,0);
po=psol_branch.point(30)%.profile(:,1)%.point(50);
Fode45=@(t,X)funcs.wrap_rhs(X,po.parameter);
rng(1)
simulode45=ode45(Fode45,[0,200],po.profile(:,1)+0.0*rand(length(pt.x),1),odeset('RelTol',1e-6,'AbsTol',1e-6));
%

figure(88)
clf;
hold on
plot(simulode45.x,simulode45.y,'LineWidth',2)