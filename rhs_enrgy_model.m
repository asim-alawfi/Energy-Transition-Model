%%
% System of ODEs for the energy model equations 
parnames = {'a','b','sigma','mu_n','mu_r','tau_delay'};
cind = [parnames; num2cell(1:length(parnames))];
in = struct(cind{:});
par([in.a,      in.b,     in.sigma,         in.mu_n,     in.mu_r,  in.tau_delay])=... 
    [0.8,         1,        0.0,             0.2,           0.4,        0  ];
ndim=2;
%% 
format short

%%y=[N;R]
%Define the Energy model, modified from Montbrio et al 2015 equations E1
% Gamma=sqrt(Delta);
% dNdt= N.*( b- b.*N -b.*R - mu_n) - sigma.*N.*R;%Lambda.*((1-w)+(1-gamma).*(w-R))-N; %% need to check the term gamma*w or just gamma
% dRdt= R.*(a- a.*N -a.*R - mu_r)  + sigma.*N.*R;%Lambda.*(1-R)-e*R;
F=@(t,y,par,sigma) [y(1,:).*( par(in.b)- par(in.b).*y(1,:) -par(in.b).*y(2,:) - par(in.mu_n)) - sigma.*y(1,:).*y(1,:);...
    y(2,:).*(par(in.a)- par(in.a).*y(1,:) -par(in.a).*y(2,:) - par(in.mu_r))  + sigma.*y(1,:).*y(2,:)];

% to use with your functions:
rhs=@(y,sigma,t) F(t,y,par,sigma); 

