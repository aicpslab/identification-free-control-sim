function Ctrl = ControlAlgorithm1(Xf,Pu,Px,C,U0,X0,Z0,X1,c3,gamma0,beta0,beta1,flag)
yalmip('clear');
[nu,ell] = deal(size(U0,1), size(U0,2));

% Parameters (runtime inputs)
x = sdpvar(size(X0,1),1);
z = sdpvar(size(Z0,1),1);
r = sdpvar(size(C,1),1);
% Decision variables
delta0 = sdpvar(1,1);
delta1 = sdpvar(1,1);

Wg     = sdpvar(ell,ell,'symmetric');      % lifting for CLF
pk     = sdpvar(ell,1);
gk     = sdpvar(ell,1);
u      = sdpvar(nu,1);

Constr = [];
%% ---- (a) Cost: squared norm (fast; avoids SOCP norm) ----

J =(C*X1*pk-r)'*(C*X1*pk-r) + beta0*delta0^2+ beta1*delta1^2;

%% ---- (b) CLF constraint: squared form ----
V0 = trace((C*X1)'*(C*X1)*Wg)- 2*r'*C*X1*gk + r'*r;
Constr = [Constr; [Wg, gk; gk', 1] >= 0];
Constr = [Constr; V0 - (1-c3)*(C*x - r)'*(C*x - r) <= delta0 ];

%% ---- (c) CBF constraint ----
Constr = [Constr; -Px.A*X1*gk + Px.b >= (1-gamma0+delta1)*(-Px.A*x + Px.b)];
Constr = [Constr; gamma0-1 <= delta1; delta1 <= gamma0];

%% ---- (d) CIS constraint (optional) ----
if flag
    Constr = [Constr; Xf.A*X1*gk <= Xf.b];
end

%% ---- (e) Input constraint ----
Constr = [Constr; Pu.A*u <= Pu.b];

%% ---- (f) Data-driven equalities ----
Constr = [Constr; [U0;X0;Z0;ones(1,ell)]*gk == [u;x;z;1]];
Constr = [Constr; [U0;C*X0;Z0;ones(1,ell)]*pk == [u;r;z;1]];
%%
ops = sdpsettings('solver','mosek','verbose',0);

Ctrl=optimizer(Constr, J, ops, {x,z,r}, u);
end
