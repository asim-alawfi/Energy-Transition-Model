clear;
rhs_enrgy_model

  base=[pwd(),filesep(),'..',filesep(),'MyBifurcation_functions',filesep()];
 addpath([base,filesep,'MyFinalwork_forCW2_testing'])
%     [base,filesep,'ddebiftool_extra_psol'],...
%     [base,filesep,'ddebiftool_utilities'],...
%     [base,filesep,'ddebiftool_extra_rotsym'],...
%     [base,filesep,'ddebiftool_extra_nmfm'],...
%     [base,filesep,'ddebiftool_extra_symbolic'],...
%     [base,filesep,'ddebiftool_coco']);
%hj=1e-6; tol=1e-8; maxit=10; 
%%
hj=1e-6; tol=1e-7; maxit=10; nmax=100;
T=2*pi/pi; N=300;
%%         plotting fixed points when a=0 
% r=linspace(0,2.25,200);
% g=@(r) pi^2.*r.^4-Jbar*r.^3-eta*r.^2-((Delta+Gamma.*r).^2)/(4*pi^2);
% figure(10)
% plot(r,g(r));
% grid on;
% % the plot shows that there are zeros near r=[1.05,0.45,0.095]
% r0=[1.05,0.45,0.095];
% % to find the valuse of v we subtitute r0 in g2 ( q(3) in CW2).
% g2=@(r) -(Delta+Gamma*r)/(2*pi*r);
% v0=[g2(r0(1)),g2(r0(2)),g2(r0(3))];
%% the initial guesess for fixed point are
uini=[0;0.5]%[r0;v0];
%%  Find fixed point when a=0 
M=@(u,sigma)MyIVRk(@(t,u)rhs(u,sigma,t),u,[0,T],N);  % defined flow map M with parameter value a
dMa=@(u,sigma)MyJacobian(@(u)M(u,sigma),u,hj);
sigma=0;
Fa=@(u)M(u,sigma)-u; 
dFa=@(u)MyJacobian(Fa,u,hj);
n=2;
x_star=zeros(n,1)%length(uini)); % fixed points of F when a=0
conv=false(1,length(uini));
for i=1:1%length(uini)
[x_star(:,i),conv(:,i),~]=MySolve(Fa,uini(:,i),dFa,tol,maxit) ;
eigv(:,i)=eig(MyJacobian(@(u)M(u,sigma),x_star(:,i),hj)); % check stability
end

% ud=x_star(:,1); % stabel
% uc=x_star(:,2);% saddle
% ub=x_star(:,3); % stabel
%% Stability of x_star
ev=NaN(2,1);
for i=1:1%length(x_star)
ev(:,i)=eig(MyJacobian(@(u)Fa(u),x_star(:,i),hj));
end
is_sink=all(abs(ev)<=1,1)
%% TrackCurve along ud
F=@(u)M(u(1:2),u(3))-u(1:2);
dF=@(u)MyJacobian(F,u,hj);
us=[x_star(:,1);sigma]; 
u_tg=[0;0;1];
s=0.01;
nm=100;
Mlist=MyTrackCurve(F,dF,us,u_tg,s,nm);
figure(1)
clf;
hold on
plot(Mlist(3,:),Mlist(1,:),'.-')

set(gca,'fontweight','bold');
grid on 
xlabel('a');  
ylabel('r');
title('Track Curve for Neural Model with periodic Foring')
%%  Stability  
%(i) find the eignvalues correspoiding to the equilibria ylist
plist=[Mlist(1,:);Mlist(2,:)] ; %Mlist where the third rwo is removed
eginlist=NaN(n,length(Mlist));
par=Mlist(3,:); % needed to evaluate the Jacobian at of Mlist(i) with a varies (%plist)
eignMlist=NaN(2,length(Mlist)); %eigvalues for Mlist(''plit'')
for i=1:100%length(Mlist(3,:))
    a=par(i); % vary a
    f=@(y)M(y,a); 
    J=@(y)MyJacobian(f,y,hj);
    eignMlist(:,i)=eig(J(Mlist(1:2,i)));
end
sink=all(abs(eignMlist)<=1,1);
saddle=any(abs(eignMlist)<=1,1) & any(abs(eignMlist)>1,1); 
source=all(abs(eignMlist)>1,1);
%% 
apd_loc=find(diff(sink)); % identify the location of the value of 'a' where the stability changes
uf=Mlist(1:2,apd_loc);  %  point when a=apd1 (meets conditions mentioned in the cw2).
upd1=[uf(1);Mlist(3,apd_loc)];  % 
%% TrackCurve along ub (x_3)
a=0;
us2=[x_star(:,3);a]; 
Mlist2=MyTrackCurve(F,dF,us2,u_tg,s,nm);
figure(22)
clf;
hold on
plot(Mlist2(3,:),Mlist2(1,:),'.-')
set(gca,'fontweight','bold');
grid on 
xlabel('a');  
ylabel('r');
title('Track Curve for Neural Model with periodic Foring')
%% plot stability
figure(202)
clf;
hold on
plot(par,plist(1,:));
%
clrs=colormap(lines);
ldeco={'linewidth',1};
% plot branch along ud
plsink=plot(par(sink),plist(1,sink),'o','color',clrs(1,:),ldeco{:});
plsaddle=plot(par(saddle),plist(1,saddle),'+','color',clrs(2,:),ldeco{:});
plsource=plot(par(source),plist(1,source),'s','color',clrs(3,:),ldeco{:});
%%  Stability  
%(i) find the eignvalues correspoiding to the equilibria ylist
plist2=[Mlist2(1,:);Mlist2(2,:)] ; %Mlist where the third rwo is removed
eginlist2=NaN(length(ub),length(plist2));
par2=Mlist2(3,:); % needed to evaluate the Jacobian at of Mlist(i) with a varies (%plist)
eignMlist2=NaN(length(ub),length(plist2)); %eigvalues for Mlist(''plit'')
uns_source2=NaN(length(ub),length(plist2)); % consists of unstable source the points(vectors)   
stable2=NaN(length(ub),length(plist2));     % ""       "" stable points
uns_saddle2=NaN(length(ub),length(plist2)); %  "" "" saddle points
%sink=false(1,length(plist)); % need later this for plotting 
%source=false(1,length(plist));% need latter for plotting
%saddle=false(1,length(plist));
for i=1:length(plist2)
    a=par2(i); % vary a
    f2=@(y)M(y,a); 
    J2=@(y)MyJacobian(f2,y,hj);
    eignMlist2(:,i)=eig(J2(plist2(:,i)));
end
sink2=all(abs(eignMlist2)<=1,1);
saddle2=any(abs(eignMlist2)<=1,1) & any(abs(eignMlist2)>1,1) ;
source2=all(abs(eignMlist2)>1,1);

                                    %% TrackCurve along uc (x_2)a=0;
us3=[x_star(:,2);0]; %a=0
Mlist3=MyTrackCurve(F,dF,us3,u_tg,s,nm);
figure(23)
clf;
hold on
plot(Mlist3(3,:),Mlist3(1,:),'.-')
set(gca,'fontweight','bold');
grid on 
xlabel('a');  
ylabel('r');
title('Track Curve for Neural Model with periodic Foring')
%%  Stability  
%(i) find the eignvalues correspoiding to the equilibria ylist
plist3=[Mlist3(1,:);Mlist3(2,:)] ; %Mlist where the third rwo is removed
eginlist3=NaN(length(ub),length(plist3));
par3=Mlist3(3,:); % needed to evaluate the Jacobian at of Mlist(i) with a varies (%plist)
eignMlist3=NaN(length(ub),length(plist3)); %eigvalues for Mlist(''plit'')
uns_source3=NaN(length(ub),length(plist3)); % consists of unstable source the points(vectors)   
stable3=NaN(length(ub),length(plist3));     % ""       "" stable points
uns_saddle3=NaN(length(ub),length(plist3)); %  "" "" saddle points
%sink=false(1,length(plist)); % need later this for plotting 
%source=false(1,length(plist));% need latter for plotting
%saddle=false(1,length(plist));
for i=1:length(plist3)
    a=par3(i); % vary a
    f3=@(y)M(y,a); 
    J3=@(y)MyJacobian(f3,y,hj);
    eignMlist3(:,i)=eig(J3(plist3(:,i)));
end
sink3=all(abs(eignMlist3)<=1,1);
saddle3=any(abs(eignMlist3)<=1,1) & any(abs(eignMlist3)>1,1) ;
source3=all(abs(eignMlist3)>1,1);

                          %% Add Stability Information (a,r)-plane
% plot stability in (beta,I)-plane
figure(202)
clf;
hold on
plot(par,plist(1,:));
hold on
plot(par2,plist2(1,:));
hold on
plot(par3,plist3(1,:));
%
clrs=colormap(lines);
ldeco={'linewidth',1};
% plot branch along ud
plsink=plot(par(sink),plist(1,sink),'o','color',clrs(1,:),ldeco{:});
plsaddle=plot(par(saddle),plist(1,saddle),'+','color',clrs(2,:),ldeco{:});
plsource=plot(par(source),plist(1,source),'s','color',clrs(3,:),ldeco{:});
%%
% plot branch along ub
plsink2=plot(par(sink2),plist2(1,sink2),'o','color',clrs(1,:),ldeco{:});
plsaddle2=plot(par(saddle2),plist2(1,saddle2),'+','color',clrs(2,:),ldeco{:});
plsource2=plot(par(source2),plist2(1,source2),'s','color',clrs(3,:),ldeco{:});
%%
% plot branch along uc

plsink3=plot(par3(sink3),plist3(1,sink3),'o','color',clrs(1,:),ldeco{:});
plsaddle3=plot(par3(saddle3),plist3(1,saddle3),'+','color',clrs(2,:),ldeco{:});
plsource3=plot(par3(source3),plist3(1,source3),'s','color',clrs(3,:),ldeco{:});
%
eqlabels={'sink','saddle','source'};
pleq=[plsink;plsaddle;plsource;plsink2;plsaddle2;plsource2;plsink3;plsaddle3;plsource3];
legend(pleq,eqlabels,'location','northwest')
set(gca,'fontweight','bold');
grid on 
xlabel('a');  
ylabel('r');
title('Stability along three branches (a,r)-plane')
                      %% Add Stability Information (a,v)-plane
figure(303)
clf;
hold on
plot(par(1:100),plist(2,100));
%%
hold on
plot(par2,plist2(2,:));
hold on
plot(par3,plist3(2,:));

clrs=colormap(lines);
ldeco={'linewidth',1};
%%
plsink=plot(par(sink),plist(2,sink),'o','color','r');
%%
plsaddle=plot(par(saddle),plist(2,saddle),'+','color','r',ldeco{:});
plsource=plot(par(source),plist(2,source),'s','color','b',ldeco{:});
%%
%
plsink2=plot(par2(sink2),plist2(2,sink2),'o','color',clrs(1,:),ldeco{:});
plsaddle2=plot(par2(saddle2),plist2(2,saddle2),'+','color',clrs(2,:),ldeco{:});
plsource2=plot(par2(source2),plist2(2,source2),'s','color',clrs(3,:),ldeco{:});
%
%
plsink3=plot(par3(sink3),plist3(2,sink3),'o','color',clrs(1,:),ldeco{:});
plsaddle3=plot(par3(saddle3),plist3(2,saddle3),'+','color',clrs(2,:),ldeco{:});
plsource3=plot(par3(source3),plist3(2,source3),'s','color',clrs(3,:),ldeco{:});
%
eqlabels={'sink','saddle','source'};
pleq=[plsink;plsaddle;plsource;plsink2;plsaddle2;plsource2;plsink3;plsaddle3;plsource3];
legend(pleq,eqlabels,'location','northwest')
set(gca,'fontweight','bold');
grid on 
xlabel('a');  
ylabel('v');
title('Stability along three branches (a,v)-plane')
%% 
%%
%% Define the flow map
%M=@(u)MyIVP(@(t,u)rhs(u,a,t),u,[0,T],N); 
%% Q1)-(c) Initialise 
%apd_loc=find(diff(sink)); % identify the location of the value of 'a' where the stability changes
%uf=Mlist(1:2,apd_loc);  % stable point when a=apd1 (as named in the cw2).
uh=[0.01;0.01]; 
u0=ud+uh; % use this initial point for iteration the map M, s.t u_j+1=M(u_j)
apd=Mlist(3,apd_loc); 
ad=2.55%apd+0.5; % chose a value of a such that apd < a < 3.
Nn=5000; % number of iteration the map M.
%Mu=@(u)MyIVRk(@(t,u)rhs(u,ad,t),u,[0,T],300); 
% Lyapunov Exponents 
[lambda,nRdiag,xy]=LyapunovQR(@(u)M(u,ad),u0,Nn);
Rdiag=nRdiag(:,101:Nn); % to remove the first 100 points
x=xy(:,101:Nn);
N=Nn-100  % to remove the first 100 points
figure(89)
clf;
hold on
%plot trajectories
plot(x(1,:),x(2,:),'.',ldeco{:})
grid on
set(gca,'fontweight','bold');
xlabel('r');
ylabel('v');
title(' trajectory in (r,v)-plane')
%%
figure(890)
clf;
hold on
%for i=1:N
%plot trajectories
plot(1:N,x)
grid on
set(gca,'fontweight','bold');
xlabel('iter');
ylabel('r-v');
title(' trajectory againest iterate number ')
%% the logarithms of the  entries of the Rdiag
figure(40)
clf;
hold on
plot(1:N,log(Rdiag(1,1:N)),'.-',1:N,log(Rdiag(2,1:N)),'.-')
grid on
xlabel('iterate');
ylabel('log(Rdiag)');
%% The cumulative means
figure(50)
clf;
hold on
plot(1:N,cumsum(log(Rdiag),2)./repmat((1:N),1))
grid on
set(gca,'fontweight','bold');
xlabel('iterate');
ylabel('comulative mean of log(Rdiag)');
title('Lyapunov exponent')
%%
save('cw21.mat')
%%
figure(500)
clf;
hold on
plot(1:N,cumsum(log(Rdiag),2)./repmat((1:N),1))%cumsum(log(Rdiag),2)./repmat((1:N),1))
grid on
set(gca,'fontweight','bold');
xlabel('iterate');
ylabel('comulative mean of log(Rdiag)');
title('Lyapunov exponent')
