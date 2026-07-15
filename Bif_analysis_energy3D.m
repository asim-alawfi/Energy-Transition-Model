
clear;
base=[pwd(),filesep(),'..',filesep(),'DDE_Biftool_Nov2025',filesep()];
ddebiftool_path(base,'sym');
format compact
format short 
%%
parnames = {'r','k','beta','h','mu',...
            'd','gamma','sigma','aR','KR',...
            'ep','delta','phi',...
            'tau'};
cind = [parnames; num2cell(1:length(parnames))];
in   = struct(cind{:});
syms(parnames{:});


par([in.r,  in.k, in.beta, in.h, in.mu, in.d,  in.gamma, in.sigma, in.aR, in.KR, in.ep, in.delta, in.phi, in.tau])=...
    [0.7,  1.6,   1.5,   2,    0.2,   0.4,  2.0,      0.1,      0.4,   0.6,    2.0,   1.0,      0.5,    0 ];

ndim=3;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create initial point for continuation
p = par;                        
rhs = @(u) [ ...
    p(in.r)*u(1)*(1-u(1)/p(in.k)) - p(in.beta)*u(1)*u(2)/(1+p(in.h)*u(1)) - p(in.mu)*u(1)*u(3); ...
    u(2)*(-p(in.d) + p(in.gamma)*u(1)/(1+p(in.h)*u(1)) + p(in.sigma)*u(3)) + p(in.aR)*u(2)*(1-u(2)/p(in.KR)); ...
    p(in.ep)*u(1) - p(in.delta)*u(3) - p(in.phi)*u(2)*u(3) ];

% initial guess from the planar interior equilibrium + a carbon estimate
N0 = p(in.d)/(p(in.gamma) - p(in.d)*p(in.h));           % = 0.8  (planar N*)
R0 = p(in.r)*(1+p(in.h)*N0)/p(in.beta)*(1 - N0/p(in.k));% planar R* guess
C0 = p(in.ep)*N0/p(in.delta);                            % C from Cdot=0, R small

u0 = [N0; R0; C0];
opts = optimoptions('fsolve','Display','off','TolFun',1e-12,'TolX',1e-12);
ustar = fsolve(rhs, u0, opts);
fprintf('E* = (N,R,C) = (%.6f, %.6f, %.6f)\n', ustar);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u0=[0.220867; 0.512335; 0.351653];
stst=dde_stst_create('parameter',par,...
    'x',u0);
    parbd={'min_bound',[in.gamma, 0;  in.k,0.01; in.tau,0],...
    'max_bound',   [in.gamma, 3;     in.k, 5; in.tau,6],...
    'max_step',    [in.gamma, 0.001;   in.k, 0.05; in.tau, 0.1; 0,0.1]};
%%
funcs=set_symfuncs(@sym_energy_carbon_delay,'sys_tau',@ () in.tau);
%%
[eq_branch1,suc1]=SetupStst(funcs,...
    'parameter',par,'x',stst.x,...
    'contpar',in.tau,'step',0.1,'dir',[],parbd{:},'print_residual_info',true);
figure(1)
clf; hold on
eq_branch1=br_contn(funcs,eq_branch1,400);
eq_branch1=br_rvers(eq_branch1);
eq_branch1=br_contn(funcs,eq_branch1,300);
%%
[eq_branch1_withbif,testfuncs,indices,biftype]=LocateSpecialPoints(funcs,eq_branch1)
[eq_branch1_with_stab,nunst_st,eq_brbifs,biflocation]=...
    MonitorChange(funcs,eq_branch1_withbif,'min_iterations',5);

%%
xpar_br1=arrayfun(@(x)x.parameter(in.tau),eq_branch1_with_stab.point);
ypar_br1=arrayfun(@(x)x.x(1),eq_branch1_with_stab.point);
figure(1)
clf
hold on 
plot(xpar_br1(nunst_st==0), ypar_br1(nunst_st==0), 'k-','Linewidth',3,'MarkerSize',10);
plot(xpar_br1(nunst_st==1), ypar_br1(nunst_st==1), 'g--','Linewidth',3,'MarkerSize',10);
plot(xpar_br1(nunst_st>=2), ypar_br1(nunst_st>=2), 'k--','Linewidth',3,'MarkerSize',10);
plot(xpar_br1(biflocation(1)), ypar_br1(biflocation(1)), 'k.','Linewidth',2,'MarkerSize',40);
%plot(xpar_br1(biflocation(2)), ypar_br1(biflocation(2)), 'kx','Linewidth',2,'MarkerSize',15);
%%
[psol_branch,sucp]=SetupPsol(funcs,eq_branch1_with_stab,biflocation(1),'degree',6,'intervals',120,...
    'max_step',[in.tau,0.1],'max_bound',[in.tau,6],'print_residual_info',true);
figure(1);
hold on%clf;ax2=gca;
[psol_branch,s,f,r]=br_contn(funcs,psol_branch,70);
%%
[psol_branchref,nunst_po,bif,p_bif]=MonitorChange(funcs,psol_branch,...
    'range',2:length(psol_branch.point),'printlevel',1,'print_residual_info',0,'exclude_trival',true,...
    'min_iterations',5);
%nunst_po(1)=0;
%%
tau_x=arrayfun(@(x)x.parameter(in.tau),psol_branchref.point);
max_p=arrayfun(@(x)max(x.profile(1,:)),psol_branchref.point);
min_p=arrayfun(@(x)min(x.profile(1,:)),psol_branchref.point);

%%
nunst_po(1)=0;
figure(1)
clf
hold on 
plot(xpar_br1(nunst_st==0), ypar_br1(nunst_st==0), 'k-','Linewidth',3,'MarkerSize',10);
plot(xpar_br1(nunst_st>=1), ypar_br1(nunst_st>=1), 'k--','Linewidth',3,'MarkerSize',10);
plot(xpar_br1(biflocation(1)), ypar_br1(biflocation(1)), 'k.','Linewidth',2,'MarkerSize',40);

xlabel('$\tau$','Interpreter','latex','FontSize',18);
ylabel('$N$','Interpreter','latex','FontSize',18);
title('Stability Analysis');
plot(tau_x,max_p,'k.',tau_x,min_p,'k.','MarkerSize',10)

%%
    parbd={'min_bound',[in.gamma, 0;  in.k,0.01; in.tau,0],...
    'max_bound',   [in.gamma, 3;     in.k, 20; in.tau,15],...
    'max_step',    [in.gamma, 0.001;   in.k, 0.1; in.tau, 0.1; 0,0.1]};
[fhopf,hopf_branchd,such1]=SetupHopf(funcs,eq_branch1_with_stab,biflocation(1),'contpar',[in.k,in.tau],...
    'dir',in.k,'step',0.05,'outputfuncs',true,parbd{:},'use_tangent',true);
figure(5)
clf;
hold on
[hopf_branchd,s1,f1,r1]=br_contn(fhopf,hopf_branchd,1000);
hopf_branchd=br_rvers(hopf_branchd);
[hopf_branchd,s11,f11,r11]=br_contn(fhopf,hopf_branchd,300);
xlabel('k');ylabel('d');
[hopf_branchdref,hf_tests,hf_bifs,hf_bifind]=MonitorChange(fhopf,hopf_branchd,...
    'printlevel',1,'min_iterations',5,'exclude_trivial',true)
%%
kh1_x=arrayfun(@(x)x.parameter(in.k),hopf_branchdref.point);
tauh1_y=arrayfun(@(x)x.parameter(in.tau),hopf_branchdref.point);
figure(3)
clf
hold on
plot(kh1_x,tauh1_y,'k-','LineWidth',2)
xlabel('$k$','FontSize',18,'Interpreter','latex')
ylabel('$\tau$','FontSize',18,'Interpreter','latex')
legend('Hopf','Saddle-node','FontSize',18,'Location','best')

%%
figure(3)
clf
hold on
teal = [0 0.43 0.43];   navy = [0.08 0.16 0.35];   gold = [0.80 0.55 0.10];

tau_inf = 2.01;                              % large-k asymptote (Prop. 6.2)
yTop    = 7;

% shaded oscillatory region ABOVE the Hopf curve
xs = kh1_x(:).';  ys = tauh1_y(:).';
h_fill = fill([xs, fliplr(xs)], [ys, yTop*ones(size(ys))], teal, ...
     'FaceAlpha',0.07,'EdgeColor','none');

% asymptote tau_infty (dashed gold)
h_as = yline(tau_inf,'--','Color',gold,'LineWidth',1.6);

% the Hopf curve
h_hopf = plot(kh1_x, tauh1_y,'k-','LineWidth',2.5);

% baseline point (k, tau_c) = (1.6, 3.40)
tau_c_base = tauh1_y(find(kh1_x>=1.6,1));     % or hard-code 3.40

text(9.5, 1.2,{'steady transition','(below curve)'},...
     'Color',navy,'FontSize',12)

% oscillatory label -> move UP into the shaded region (above the curve)
text(6, 4.0,'oscillatory transition (above)',...
     'Color',teal,'FontSize',12)
text(0.8, 6.2,'$\tau_c\!\to\!\infty$ as $k\!\downarrow\! k_f$',...
     'Color','k','FontSize',11,'Interpreter','latex')

xlabel('effective nonrenewable capacity $k$','FontSize',16,'Interpreter','latex')
ylabel('critical delay $\tau_c$ (years)','FontSize',16,'Interpreter','latex')
xlim([0 20]); ylim([0 yTop])
set(gca,'FontSize',13,'LineWidth',0.8); box on

legend([h_hopf h_as],{'Hopf curve $\tau_c(k)$','asymptote $\tau_\infty\approx 2.01$'},...
       'FontSize',12,'Location','northeast','Interpreter','latex')

exportgraphics(gcf,'fig_hopf_curve2.pdf','ContentType','vector')

%%
figure(1)
clf
hold on
teal = [0 0.43 0.43];   navy = [0.08 0.16 0.35];
tau_c = xpar_br1(biflocation(1));

% equilibrium branch
h1 = plot(xpar_br1(nunst_st==0), ypar_br1(nunst_st==0),'k-','LineWidth',3);
h2 = plot(xpar_br1(nunst_st>=1), ypar_br1(nunst_st>=1),'k--','LineWidth',3);

% PO extrema — one handle for the legend, rest hidden
h3 = plot(tau_x(1), max_p(1),'k.','MarkerSize',10);
plot(tau_x, max_p,'k.', tau_x, min_p,'k.','MarkerSize',10,...
     'HandleVisibility','off')

% HB point
plot(tau_c, ypar_br1(biflocation(1)),'k.','MarkerSize',40,...
     'HandleVisibility','off')
text(tau_c-0.15, ypar_br1(biflocation(1))+0.10,'HB',...
     'FontSize',14,'FontWeight','bold','HorizontalAlignment','right')

% region label (only the steady one — POs now in legend)
% text(0.9, 0.62,{'stable steady','transition'},'Color',navy,...
%      'FontSize',12,'HorizontalAlignment','center')

xlabel('delay $\tau$ (years)','Interpreter','latex','FontSize',16)
ylabel('$N$ (nonrenewable, extrema)','Interpreter','latex','FontSize',16)
xlim([0 6]); ylim([0 0.7])
set(gca,'FontSize',13,'LineWidth',0.8); box on

legend([h1 h2 h3],{'stable equilibrium','unstable equilibrium','stable POs'},...
       'FontSize',12,'Location','northwest')

exportgraphics(gcf,'fig_hopf_1param.pdf','ContentType','vector')

%%
[~,it]=min(abs(xpar_br1- 4))
pi=eq_branch1_with_stab.point(it)
fun_sim=@(t,X,Ch)funcs.wrap_rhs(cat(2,X,Ch),pi.parameter);
his=@(t)dde_coll_eva(pi.profile,pi.mesh,1+t/pi.period,pi.degree);
lags=pi.parameter(in.tau);
rng(1)
sol23_test=dde23(fun_sim,lags,pi.x+0.03*rand(3,1),[0,10000], ddeset('RelTol',1e-5,'AbsTol',1e-5));
%%%%%
[~,it2]=min(abs(xpar_br1- 1.4))
pi2=eq_branch1_with_stab.point(it2)
fun_sim2=@(t,X,Ch)funcs.wrap_rhs(cat(2,X,Ch),pi2.parameter);
lags2=pi2.parameter(in.tau);
rng(2)
sol23_test2=dde23(fun_sim2,lags2,pi.x+0.03*rand(3,1),[0,10000], ddeset('RelTol',1e-5,'AbsTol',1e-5));
%%
figure
plot(sol23_test2.x,sol23_test2.y,'LineWidth',2);
%%
%% ---- simulations  --------------------------------------
[~,it]  = min(abs(xpar_br1 - 4.51));
pt      = eq_branch1_with_stab.point(it);
fun_sim = @(t,X,Ch) funcs.wrap_rhs(cat(2,X,Ch), pt.parameter);
lags    = pt.parameter(in.tau);
rng(1)
sol_osc = dde23(fun_sim, lags, pt.x + 0.03*rand(3,1), [0,10000], ...
                ddeset('RelTol',1e-5,'AbsTol',1e-5));

[~,it2]  = min(abs(xpar_br1 - 2.5));
pt2      = eq_branch1_with_stab.point(it2);
fun_sim2 = @(t,X,Ch) funcs.wrap_rhs(cat(2,X,Ch), pt2.parameter);
lags2    = pt2.parameter(in.tau);
rng(2)
sol_std  = dde23(fun_sim2, lags2, pt2.x + 0.03*rand(3,1), [0,10000], ...
                 ddeset('RelTol',1e-5,'AbsTol',1e-5));

%% ---- two-panel figure in the paper's style -----------------------------
teal = [0 0.43 0.43];  navy = [0.08 0.16 0.35];

figure('Units','centimeters','Position',[2 2 26 9])
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% (a) steady transition, tau = 1.4 < tau_c
nexttile
hold on
plot(sol_std.x, sol_std.y(1,:), 'k-', 'LineWidth',1.6)
plot(sol_std.x, sol_std.y(2,:), '-',  'Color',teal, 'LineWidth',1.6)
plot(sol_std.x, sol_std.y(3,:), '--', 'Color',navy, 'LineWidth',1.3)
yline(pt2.x(1), 'k:',            'LineWidth',0.7, 'Alpha',0.7)  % E* guides
yline(pt2.x(2), ':', 'Color',teal,'LineWidth',0.7, 'Alpha',0.7)
xlabel('$t$ (years)','Interpreter','latex','FontSize',13)
ylabel('normalized state','FontSize',13)
title('(a) $\tau=2.5<\tau_c$: steady transition', ...
      'Interpreter','latex','FontSize',12,'FontWeight','normal')
xlim([0 10000]); set(gca,'FontSize',11); box on
legend({'$N$','$R$','$C$'},'Interpreter','latex', ...
       'FontSize',11,'Location','northeast')

% (b) oscillatory transition, tau = 4.0 > tau_c
nexttile
hold on
plot(sol_osc.x, sol_osc.y(1,:), 'k-', 'LineWidth',1.6)
plot(sol_osc.x, sol_osc.y(2,:), '-',  'Color',teal, 'LineWidth',1.6)
plot(sol_osc.x, sol_osc.y(3,:), '--', 'Color',navy, 'LineWidth',1.3)
xlabel('$t$ (years)','Interpreter','latex','FontSize',13)
title('(b) $\tau=4.5>\tau_c$: oscillatory transition', ...
      'Interpreter','latex','FontSize',12,'FontWeight','normal')
xlim([0 10000]); ylim([0 0.6]); set(gca,'FontSize',11); box on

exportgraphics(gcf,'fig_timeseries2.pdf','ContentType','vector')

%%
tau_x=arrayfun(@(x)x.parameter(in.tau),psol_branchref.point);
period_p=arrayfun(@(x)x.period,psol_branchref.point);

figure(999)
plot(tau_x, period_p,'LineWidth',3)
xlabel('$\tau$','Interpreter','latex','FontSize',16)
ylabel('period','FontSize',16)

%%
save('computations_3d.mat')

%%

figure(1)
clf
tiledlayout(1,2,'TileSpacing','compact','Padding','compact')
teal = [0 0.43 0.43];   navy = [0.08 0.16 0.35];

% ---------- (a) bifurcation diagram ----------
nexttile
hold on
tau_c = xpar_br1(biflocation(1));
h1 = plot(xpar_br1(nunst_st==0), ypar_br1(nunst_st==0),'k-','LineWidth',3);
h2 = plot(xpar_br1(nunst_st>=1), ypar_br1(nunst_st>=1),'k--','LineWidth',3);
h3 = plot(tau_x(1), max_p(1),'k.','MarkerSize',10);
plot(tau_x, max_p,'k.', tau_x, min_p,'k.','MarkerSize',10,'HandleVisibility','off')
plot(tau_c, ypar_br1(biflocation(1)),'k.','MarkerSize',40,'HandleVisibility','off')
text(tau_c-0.15, ypar_br1(biflocation(1))+0.10,'HB',...
     'FontSize',14,'FontWeight','bold','HorizontalAlignment','right')
xlabel('delay $\tau$ (years)','Interpreter','latex','FontSize',16)
ylabel('$N$ (nonrenewable, extrema)','Interpreter','latex','FontSize',16)
title('(a)','FontSize',13,'FontWeight','normal')
xlim([0 6]); ylim([0 0.7])
set(gca,'FontSize',13,'LineWidth',0.8); box on
legend([h1 h2 h3],{'stable equilibrium','unstable equilibrium','stable POs'},...
       'FontSize',12,'Location','northwest')

% ---------- (b) period along the PO branch ----------
nexttile
hold on
tau_x2    = arrayfun(@(x)x.parameter(in.tau), psol_branchref.point);
period_p  = arrayfun(@(x)x.period,            psol_branchref.point);
plot(tau_x2, period_p,'-','Color',teal,'LineWidth',3)
xlabel('delay $\tau$ (years)','Interpreter','latex','FontSize',16)
ylabel('cycle period (years)','FontSize',16)
title('(b)','FontSize',13,'FontWeight','normal')
xlim([tau_c 6]); set(gca,'FontSize',13,'LineWidth',0.8); box on

exportgraphics(gcf,'fig_hopf_1param2.pdf','ContentType','vector')