function Ctrl = ControlAlgorithm2( ...
            Xf,Pu,Px,C, ...              % polyhedra & CLF matrices
            U0,X0,Z0,X1, ...               % data matrices
            c3,gamma0,epsilon0, ...        % CLF/CBF & margin
            beta0,beta1,beta2, ...         % weights for δ0,δ1,δ2
            Lx_cell,Lc,W_scale, ...        % {L_r} for CIS, Lc for CBF (can be [])
            flag)

    % Dimensions
    [nu,t] = deal(size(U0,1), size(U0,2));  % nu: input dim, t: #columns
    qX     = size(Xf.A,1);                  % #CIS inequalities
    ny     = size(C,1);                     % output dim for cost/CLF

    % Parameters (runtime inputs)
    x = sdpvar(size(X0,1),1);               % xk
    z = sdpvar(size(Z0,1),1);               % zk
    r = sdpvar(ny,1);                       % rk

    % Decision variables
    delta0 = sdpvar(1,1);      % CLF slack
    delta1 = sdpvar(1,1);      % CBF slope slack
    delta2 = sdpvar(1,1);      % CBF margin slack

    Wg     = sdpvar(t,t,'symmetric');      % lifting for CLF
    pk     = sdpvar(t,1);      % cost regressor
    gk     = sdpvar(t,1);      % prediction regressor

    u_var  = sdpvar(nu,1);     % input

    Constr = [];

%% (a) Cost function (nominal lifted)
    J0 = (C*X1*pk-r)'*(C*X1*pk-r);
    % Total cost
    J = J0 + beta0*delta0^2 + beta1*delta1^2 + beta2*delta2^2;

    %% (b) CLF constraint (nominal lifted, with δ0)
    V0 = trace(X1'*(C'*C)*X1*Wg) - 2*r'*C*X1*gk + r'*r;
    Constr = [Constr; [Wg, gk; gk', 1] >= 0];
    Constr = [Constr; V0 - (1-c3)*(C*x - r)'*(C*x - r) <= delta0 ];
    
    %% (c) Robust CBF constraint based on X (Px plays role of C)
    wj=W_scale*[U0*gk;X0*gk;Z0*gk;1];

    if ~flag(3),epsilon0=0;end
    hxk    = -Px.A*x   + Px.b;   % h(x_k)
    
    if ~isempty(Lc)&&flag(1)
        hxnext = -Px.A*X1*gk + Px.b - norm(Lc*wj,2);
    else
        hxnext = -Px.A*X1*gk + Px.b;  % nominal h(x_{k+1})
    end

    Constr = [Constr;hxnext >= (1-gamma0+delta1)*hxk + (epsilon0+delta2)];

    % Bounds for δ1, δ2  (Algorithm 2: δ1∈[γ0-1,γ0], δ2∈[-ε0,0])
    Constr = [Constr;gamma0-1 <= delta1; delta1 <= gamma0];
    Constr = [Constr;-epsilon0 <= delta2; delta2 <= 0];
    %% (d) Robust CIS constraints (row-wise), controlled by flag
    if flag(2)
        if ~isempty(Lx_cell)&&flag(1)
            % Robust CIS:  A_X x_{k+1} + sup_Δ A_X Δ g <= b_X
            % with sup_Δ term bounded by ||L_r g||_2 ≤ s_x(r)
            for rr = 1:qX
                Constr = [Constr; Xf.A(rr,:)*X1*gk + norm(Lx_cell{rr}*wj,2) <= Xf.b(rr)];
            end
        else
            % If no Lx_cell is provided, use nominal CIS constraints
            Constr = [Constr; Xf.A*X1*gk <= Xf.b];
        end
    end

    %% (e) Input constraints
    Constr = [Constr; Pu.A*u_var <= Pu.b];

    %% (f) Data-driven equalities (Theorem 1 relations)
    Constr = [Constr; [U0; X0; Z0; ones(1,t)]*gk == [u_var; x; z; 1]];
    Constr = [Constr; [U0; C*X0; Z0; ones(1,t)]*pk == [u_var; r; z; 1]];
    %% Compile to reusable parametric solver
    ops = sdpsettings('solver','mosek','verbose',0,'cachesolvers',1);
    
    Ctrl = optimizer(Constr, J, ops, {x,z,r}, u_var);
end