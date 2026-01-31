function dydt = rhs_dde23_enrgy(t, y, Z,par)
% y(1)=N(t), y(2)=R(t)
% Z(:,1)=y(t-tau), so Z(2,1)=R(t-tau)
parnames = {'a','lambda','sigma','mu_n','mu_r','k','tau_delay'};
cind = [parnames; num2cell(1:length(parnames))];
in = struct(cind{:});
par([in.a,      in.lambda,     in.sigma,          in.mu_n,           in.mu_r,   in.k   in.tau_delay])=... 
    [0.35,         0.15,            2,              0.3,               0.6,            0.9       3 ];

N    = y(1);
R    = y(2);
Rtau = Z(2,1);
dNdt = N .* ( par(in.lambda) .* (1 - (N + R)./par(in.k)) - par(in.mu_n) ) ...
       - par(in.sigma) .* N .* Rtau;

dRdt = R .* ( par(in.a) .* (1 - (N + R)./par(in.k)) - par(in.mu_r) ) ...
       + par(in.sigma) .* N .* Rtau;
dydt = [dNdt; dRdt];
end
