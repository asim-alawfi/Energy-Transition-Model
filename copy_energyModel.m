clear;
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
%%
parnames = {'a','r','gamma','d','k','beta','mu_n','h','tau_delay'};
cind = [parnames; num2cell(1:length(parnames))];
in = struct(cind{:});
par([in.a,      in.r,        in.gamma,      in.d,        in.k,      in.beta,   in.mu_n    in.h,   in.tau_delay])=... 
    [0.00,        0.04,          0.13,        0.04,        1.6,           0.065,      0.00,       2,       0 ];
ndim=2;
stst=dde_stst_create('parameter',par,...
    'x',[0.8;0.8])%par(in.k)*(1-(par(in.mu_n)/par(in.a)))]);
%x0=[0;0.5]%zeros(ndim,1);
    parbd={'min_bound',[in.gamma, 0;  in.k,0.01; in.d,0],...
    'max_bound',   [in.gamma, 1;     in.k, 5; in.d,1],...
    'max_step',    [in.gamma, 0.001;   in.k, 0.05; in.d, 0.001; 0,0.05]};
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
eq_branch1=br_contn(funcs,eq_branch1,300);
%%
[eq_branch1ref,nunst_st,eq_brbifs,biflocation]=...
    MonitorChange(funcs,eq_branch1,'min_iterations',5);
% fprintf('1st Symmetry breaking at a=%g\neigenvalues, eigenvectors:',eq_brbifs(1).parameter(in.k));
% disp(eq_brbifs(1).stability.l0)
% disp(eq_brbifs(1).stability.v)
% [eq_branch1,bif1testfuncs,biflocation,biftypes]=LocateSpecialPoints(funcs,eq_branch,'debug',true);

%[eq_branch1,nunst_st]=br_stabl(funcs,eq_branch1,0,0);
%% plot stability
xpar_br1=arrayfun(@(x)x.parameter(in.k),eq_branch1ref.point);
ypar_br1=arrayfun(@(x)x.x(1),eq_branch1ref.point);
%biflocation=find(diff(nunst_st));%find((nunst_st>0),1,'last');
%p1=eq_1234_0.point(biflocation);
%%
% Extract stability information and plot results
%%
[psol_branch,sucp]=SetupPsol(funcs,eq_branch1ref,biflocation(1),'degree',6,'intervals',120,...
    'max_step',[in.k,0.1],'max_bound',[in.k,5]);
figure(1);
hold on%clf;ax2=gca;
[psol_branch,s,f,r]=br_contn(funcs,psol_branch,70);
%%
[psol_branchref,nunst_po,bif,p_bif]=MonitorChange(funcs,psol_branch,...
    'range',2:length(psol_branch.point),'printlevel',1,'print_residual_info',0,'exclude_trival',true,...
    'min_iterations',5);
nunst_po(1)=0;
kp_x=arrayfun(@(x)x.parameter(in.k),psol_branchref.point);
max_p=arrayfun(@(x)max(x.profile(1,:)),psol_branchref.point);
min_p=arrayfun(@(x)min(x.profile(1,:)),psol_branchref.point);
%%

figure(1)
clf
hold on 
plot(xpar_br1(nunst_st==0), ypar_br1(nunst_st==0), 'k-','Linewidth',3,'MarkerSize',10);
plot(xpar_br1(nunst_st==1), ypar_br1(nunst_st==1), 'k--','Linewidth',3,'MarkerSize',10);
plot(xpar_br1(nunst_st>=2), ypar_br1(nunst_st>=2), 'k--','Linewidth',3,'MarkerSize',10);
plot(xpar_br1(biflocation(1)), ypar_br1(biflocation(1)), 'k.','Linewidth',2,'MarkerSize',40);
plot(xpar_br1(biflocation(2)), ypar_br1(biflocation(2)), 'kx','Linewidth',2,'MarkerSize',15);

xlabel('$k$','Interpreter','latex','FontSize',18);
ylabel('$N$','Interpreter','latex','FontSize',18);
title('Stability Analysis');
plot(kp_x,max_p,'k.',kp_x,min_p,'k.','MarkerSize',10)
%% 
[fhopf,hopf_branchd,such1]=SetupHopf(funcs,eq_branch1ref,biflocation(1),'contpar',[in.k,in.d],...
    'dir',in.k,'step',0.05,'outputfuncs',true,parbd{:},'use_tangent',true);
%hopf2_branch.method.point.print_residual_info=0;
figure(5)
clf;
hold on
[hopf_branchd,s1,f1,r1]=br_contn(fhopf,hopf_branchd,70);
hopf_branchd=br_rvers(hopf_branchd);
[hopf_branchd,s11,f11,r11]=br_contn(fhopf,hopf_branchd,70);
xlabel('k');ylabel('d');
[hopf_branchdref,hf_tests,hf_bifs,hf_bifind]=MonitorChange(fhopf,hopf_branchd,...
    'printlevel',1,'min_iterations',5,'exclude_trivial',true)
%%
[ffold_d,foldd_branch,suc]=SetupFold(funcs,eq_branch1ref,biflocation(2),'contpar',[in.k,in.d],...
    'dir',in.k,'step',0.05,'outputfuncs',true,parbd{:})%,'use_tangent',true);
%hopf2_branch.method.point.print_residual_info=0;
figure(5)
[foldd_branch,s,f,r]=br_contn(ffold_d,foldd_branch,120);
foldd_branch=br_rvers(foldd_branch);
[foldd_branch,s,f,r]=br_contn(ffold_d,foldd_branch,100);
xlabel('k');ylabel('d');
[foldd_branch,ff_tests,ff_bifs,ff_bifind]=MonitorChange(ffold_d,foldd_branch,...
    'printlevel',1,'min_iterations',5,'exclude_trivial',true)
%%
kh1_x=arrayfun(@(x)x.parameter(in.k),hopf_branchdref.point);
dh1_y=arrayfun(@(x)x.parameter(in.d),hopf_branchdref.point);
kf1_x=arrayfun(@(x)x.parameter(in.k),foldd_branch.point);
df1_y=arrayfun(@(x)x.parameter(in.d),foldd_branch.point);
figure(3)
clf
hold on
plot(kh1_x,dh1_y,'k-',...
    kf1_x,df1_y,'k--','LineWidth',2)
xlabel('$k$','FontSize',18,'Interpreter','latex')
ylabel('$d$','FontSize',18,'Interpreter','latex')
legend('Hopf','Saddle-node','FontSize',18,'Location','best')
%title('(k,d)-paramter','FontSize',18,'Interpreter','latex')
%%
[fhopf2,hopf_branchg,such2]=SetupHopf(funcs,eq_branch1ref,biflocation(1),'contpar',[in.k,in.gamma],...
    'dir',in.k,'step',0.05,'outputfuncs',true,parbd{:});%,'use_tangent',true);
%hopf2_branch.method.point.print_residual_info=0;
figure(6)
clf;
hold on
[hopf_branchg,s2,f2,r2]=br_contn(fhopf2,hopf_branchg,70);
hopf_branchg=br_rvers(hopf_branchg);
[hopf_branchg,s22,f22,r22]=br_contn(fhopf2,hopf_branchg,70);
xlabel('k');ylabel('\gamma');
[hopf_branchgref,hfg_tests,hfg_bifs,hf_bifing]=MonitorChange(fhopf2,hopf_branchg,...
    'printlevel',1,'min_iterations',5,'exclude_trivial',true);
%%
[ffold_g,foldg_branch,suc]=SetupFold(funcs,eq_branch1ref,biflocation(2),'contpar',[in.k,in.gamma],...
    'dir',in.k,'step',0.05,'outputfuncs',true,parbd{:})%,'use_tangent',true);
%hopf2_branch.method.point.print_residual_info=0;
figure(5)
[foldg_branch,s,f,r]=br_contn(ffold_g,foldg_branch,120);
foldg_branch=br_rvers(foldg_branch);
[foldg_branch,s,f,r]=br_contn(ffold_g,foldg_branch,100);
xlabel('k');ylabel('d');
[foldg_branch,ffg_tests,ff_bifsg,ff_bifind_g]=MonitorChange(ffold_g,foldg_branch,...
    'printlevel',1,'min_iterations',5,'exclude_trivial',true)
kh2_x=arrayfun(@(x)x.parameter(in.k),hopf_branchgref.point);
gh2_y=arrayfun(@(x)x.parameter(in.gamma),hopf_branchgref.point);
kf2_x=arrayfun(@(x)x.parameter(in.k),foldg_branch.point);
gf2_y=arrayfun(@(x)x.parameter(in.gamma),foldg_branch.point);
%%
figure(2)
clf
hold on
plot(kh2_x,gh2_y,'k-', ...
    kf2_x,gf2_y,'k--','LineWidth',2)
xlabel('$k$','FontSize',18,'Interpreter','latex')
ylabel('$\gamma$','FontSize',18,'Interpreter','latex')
legend('Hopf','saddle-node','FontSize',18)
title('$(k,\gamma)$-paramter','FontSize',18,'Interpreter','latex')
ylim([0.08,0.2])
%%
save('Energy_model.mat')
%%
clrs=lines();
figure(99)
clf;
tiledlayout(4,5)
nexttile([4 3])
hold on
plot(xpar_br1(nunst_st==0), ypar_br1(nunst_st==0), 'k-','Linewidth',3,'MarkerSize',10);
plot(xpar_br1(nunst_st==1), ypar_br1(nunst_st==1), 'k--','Linewidth',2,'MarkerSize',10);
plot(xpar_br1(nunst_st>=2), ypar_br1(nunst_st>=2), 'k--','Linewidth',2,'MarkerSize',10);
plot(xpar_br1(biflocation(1)), ypar_br1(biflocation(1)), 'k.','Linewidth',2,'MarkerSize',40);
plot(xpar_br1(biflocation(2)), ypar_br1(biflocation(2)), 'kx','Linewidth',2,'MarkerSize',15);
plot(kp_x,max_p,'k.',kp_x,min_p,'k.','MarkerSize',10)
plot(kp_x(20),max_p(20),'.','MarkerSize',30,'Color',clrs(6,:))
plot(kp_x(40),max_p(40),'.','MarkerSize',30,'Color',clrs(5,:))
plot(kp_x(60),max_p(60),'.','MarkerSize',30,'Color',clrs(4,:))
set(gca,"Box","on",'FontSize',12,'FontWeight','bold','LineWidth',2)
legend('stable Eqs','unstable Eqs','','HB','SN','stable POs','FontSize',14,'FontWeight','normal')
xlabel('$k$','Interpreter','latex','FontSize',14,'FontWeight','bold');
ylabel('$N$','Interpreter','latex','FontSize',14,'FontWeight','bold');
title('(a)','FontSize',12,'FontWeight','normal');
nexttile([1 2])
p_x=arrayfun(@(x)x.period,psol_branchref.point);
plot(kp_x,p_x,'LineWidth',2)
set(gca,"Box","on",'FontSize',9,'FontWeight','bold','LineWidth',2)
xlabel('$k$','Interpreter','latex','FontSize',14,'FontWeight','bold');
ylabel('$period$','Interpreter','latex','FontSize',12,'FontWeight','bold');

title('(b)','FontSize',12,'FontWeight','normal');
nexttile([1,2])
hold on
po1=psol_branchref.point(20);
plot(po1.mesh*po1.period,po1.profile,'LineWidth',2)
plot(40,0.9,'.','MarkerSize',30,'Color',clrs(6,:))
xlim([0,po1.period])
set(gca,"Box","on",'FontSize',9,'FontWeight','bold','LineWidth',2)
title('(c)','FontSize',12,'FontWeight','normal');
nexttile([1,2])
hold on
po2=psol_branchref.point(40);
plot(po2.mesh*po2.period,po2.profile,'LineWidth',2)
xx=[50;0.9];
plot(40,0.9,'.','MarkerSize',30,'Color',clrs(5,:))
xlim([0,po2.period])
set(gca,"Box","on",'FontSize',9,'FontWeight','bold','LineWidth',2)
title('(d)','FontSize',12,'FontWeight','normal');
nexttile([1,2])
hold on
po3=psol_branchref.point(60);
plot(po3.mesh*po3.period,po3.profile,'LineWidth',2)
plot(100,3,'.','MarkerSize',30,'Color',clrs(4,:))
xlim([0,po3.period])
set(gca,"Box","on",'FontSize',9,'FontWeight','bold','LineWidth',2)
xlabel('$period$','Interpreter','latex','FontSize',12,'FontWeight','bold');
title('(e)','FontSize',12,'FontWeight','normal');
legend('$N$','$R$','Interpreter','latex','FontSize',14)
%%
figure(88)
clf;
tiledlayout(2,4)
nexttile([2,2])
hold on
plot(kh1_x,dh1_y,'k-',...
    kf1_x,df1_y,'k--','LineWidth',2)
set(gca,"Box","on",'FontSize',9,'FontWeight','bold','LineWidth',2)
xlabel('$k$','FontSize',14,'Interpreter','latex')
ylabel('$d$','FontSize',14,'Interpreter','latex')
%legend('HB','feasability boundary','FontSize',12,'Location','best')
title('(a)','FontSize',12,'FontWeight','normal');
ylim([0,0.08])
nexttile([2,2])
hold on
plot(kh2_x,gh2_y,'k-', ...
    kf2_x,gf2_y,'k--','LineWidth',2)
set(gca,"Box","on",'FontSize',9,'FontWeight','bold','LineWidth',2)
xlabel('$k$','FontSize',14,'Interpreter','latex')
ylabel('$\gamma$','FontSize',14,'Interpreter','latex')
legend('HB','feasability boundary','FontSize',12,'Location','best')
title('(b)','FontSize',12,'FontWeight','normal');
ylim([0.08,0.2])
%%
pt=eq_branch1ref.point(biflocation(1)-10);
Fode45=@(t,X)funcs.wrap_rhs(X,pt.parameter);
rng(1)
simulode45=ode45(Fode45,[0,250000],pt.x+0.05*rand(length(pt.x),1),odeset('RelTol',1e-6,'AbsTol',1e-6));
%

figure(88)
clf;
hold on
plot(simulode45.x,simulode45.y,'LineWidth',2)
pt=eq_branch1.point(50+0);
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
% [eq_branch2,suc0]=SetupStst(funcs,...
%     'parameter',pt.parameter,'x',simulode45.y(:,end),...
%     'contpar',in.k,'step',0.05,parbd{:});
% figure(1)
%  hold on
% eq_branch2=br_contn(funcs,eq_branch2,1000);
% eq_branch2=br_rvers(eq_branch2);
% %%
% eq_branch2=br_contn(funcs,eq_branch2,120);
% [eq_branch2,nunst_st2]=br_stabl(funcs,eq_branch2,0,0);
% %% plot stability
% [eq_branch2,bif1testfuncs,biflocation2,bif1types]=LocateSpecialPoints(funcs,eq_branch2,'debug',true);
% xpar_br2=arrayfun(@(x)x.parameter(in.k),eq_branch2.point);
% ypar_br2=arrayfun(@(x)x.x(1),eq_branch2.point);
%biflocation2=find(diff(nunst_st2));%find((nunst_st>0),1,'last');
%p1=eq_1234_0.point(biflocation);
%[eq_branch2,nunst_st2]=br_stabl(funcs,eq_branch2,0,0);
%%
% biflocation2=find(diff(nunst_st2));
% 
% %%
% % Extract stability information and plot results
% figure(1);
% clf
% hold on 
% plot(xpar_br2(nunst_st2==0), ypar_br2(nunst_st2==0), 'k.','Linewidth',2,'MarkerSize',10);
% plot(xpar_br2(nunst_st2>=1), ypar_br2(nunst_st2>=1), 'r.','Linewidth',2,'MarkerSize',10);
% plot(xpar_br2(biflocation2(1)), ypar_br2(biflocation2(1)), 'b.','Linewidth',2,'MarkerSize',30);
% plot(xpar_br2(biflocation2(2)), ypar_br2(biflocation2(2)), 'g.','Linewidth',2,'MarkerSize',30);
%
% plot(xpar_br1(nunst_st==0), ypar_br1(nunst_st==0), 'k.','Linewidth',2,'MarkerSize',10);
% plot(xpar_br1(nunst_st==1), ypar_br1(nunst_st==1), 'r.','Linewidth',2,'MarkerSize',10);
% plot(xpar_br1(biflocation), ypar_br1(biflocation), 'r.','Linewidth',2,'MarkerSize',20);
% xlabel('k');
% ylabel('y_{1}');
% title('Stability Analysis of Equilibria');
% %%
% [psol_branch,suc]=SetupPsol(funcs,eq_branch2,biflocation2(1),'degree',6,'intervals',120,...
%     'max_step',[in.k,0.1],'max_bound',[in.k,9]);
% figure(19);
% hold on%clf;ax2=gca;
% [psol_branch,s,f,r]=br_contn(funcs,psol_branch,70);
% %%
% 
% %%
% % look at the period along the branch:
% figure(4); clf;
% Plot2dBranch(psol_branch,'y',@(p)p.period,'funcs',funcs);
% xlabel('a21');
% ylabel('period');

%plot(xpar_1234(biflocation), ypar_1234(biflocation), 'r.','Linewidth',2,'MarkerSize',20);
%%
%save("steady_state_br.mat")
%%
% pt2=eq_branch2.point(end);
% Fode45=@(t,X)funcs.wrap_rhs(X,pt2.parameter);
% rng(1)
% simulode45=ode45(Fode45,[0,2000],pt2.x+0.01*rand(length(pt2.x),1),odeset('RelTol',1e-6,'AbsTol',1e-6));
% %
% 
% figure(88)
% clf;
% hold on
% plot(simulode45.x,simulode45.y,'LineWidth',2)
% %%
% [eq_branch3,suc0]=SetupStst(funcs,...
%     'parameter',pt.parameter,'x',simulode45.y(:,end),...
%     'contpar',in.k,'step',0.1,parbd{:});
% figure(1)
%  hold on
% eq_branch3=br_contn(funcs,eq_branch3,1000);
% eq_branch3=br_rvers(eq_branch3);
% eq_branch3=br_contn(funcs,eq_branch3,400);
% [eq_branch3,nunst_st3]=br_stabl(funcs,eq_branch3,0,0);
% 
% %%
%   parbd2={'min_bound',[in.gamma, 0;  in.k,0.0; in.d,0],...
%     'max_bound',   [in.gamma, 5;     in.k, 5; in.d,1],...
%     'max_step',    [in.gamma, 0.1;   in.k, 0.05; in.d, 0.02; 0,0.05]};
% %indhopf2=indbifhopf(find(strcmp('hoho',bifhopftypes),1,'first'));
% [fhopf,hopf_branchd,suc]=SetupHopf(funcs,eq_branch2,biflocation2(1),'contpar',[in.k,in.d],...
%     'dir',in.k,'step',0.05,'outputfuncs',true,'max_bound',[in.k,9],'use_tangent',true);
% %hopf2_branch.method.point.print_residual_info=0;
% figure(5)
% clf;
% hold on
% [hopf_branchd,s,f,r]=br_contn(fhopf,hopf_branchd,70);
% hopf_branchd=br_rvers(hopf_branchd);
% [hopf_branchd,s,f,r]=br_contn(fhopf,hopf_branchd,70);
% xlabel('a21');ylabel('\tau_s');
% [hopf_branchd,hf_sb3tests,hf_sb3bifs,hf_sb3bifind]=MonitorChange(fhopf,hopf_branchd,...
%     'printlevel',1,'min_iterations',5,'exclude_trivial',true)

%%
% [ffold,fold2_branch,suc]=SetupFold(funcs,eq_branch2,biflocation2(2),'contpar',[in.k,in.d],...
%     'dir',in.k,'step',0.05,'outputfuncs',true,'max_bound',[in.k,9])%,'use_tangent',true);
% %hopf2_branch.method.point.print_residual_info=0;
% figure(5)
% [fold2_branch,s,f,r]=br_contn(ffold,fold2_branch,100);
% fold2_branch=br_rvers(fold2_branch);
% [fold2_branch,s,f,r]=br_contn(ffold,fold2_branch,100);
% xlabel('a21');ylabel('\tau_s');
% %%
% [fold2_branch,ff_tests,ff_bifs,ff_bifind]=MonitorChange(ffold,fold2_branch,...
%     'printlevel',1,'min_iterations',5,'exclude_trivial',true)
%%
%[hopf2_branch,nunst_hf]=br_stabl(fhopf,hopf2_branch,0,0);
% po=psol_branch.point(30)%.profile(:,1)%.point(50);
% Fode45=@(t,X)funcs.wrap_rhs(X,po.parameter);
% rng(1)
% simulode45=ode45(Fode45,[0,200],po.profile(:,1)+0.0*rand(length(pt.x),1),odeset('RelTol',1e-6,'AbsTol',1e-6));
% %
% 
% figure(88)
% clf;
% hold on
% plot(simulode45.x,simulode45.y,'LineWidth',2)

