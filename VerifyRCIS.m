function report = VerifyRCIS(Xf, U0, X0, Z0, X1, Pu, Pz, Lx_cell, W_scale, eta)
%% VerifyRCIS (simplified)
% Check (R)CIS feasibility over all vertex pairs of (Pgx, Pgz) in g-space.

%% Dimensions
[ell, nu, nx, nz] = deal(size(U0,2), size(U0,1), size(X0,1), size(Z0,1));

%% Pseudoinverse map
W = [U0; X0; Z0; ones(1,ell)];
%% Ray-based offset deltaRay
deltaRay = zeros(size(Xf.A,1),1);

if eta > 0 && ~isempty(Xf.R)
    Rx = Xf.R;
    denom = vecnorm(Rx, 2, 2);
    Rx = Rx ./ denom;

    for i = 1:size(Rx,1)
        xj = eta * Rx(i,:).';

        gj =sdpvar(ell,1);
        Constr=W*gj==[zeros(nu,1); xj; zeros(nz,1);0];

        ops = sdpsettings('solver','mosek','verbose',0);

        if ~isempty(Lx_cell)
            w_raw = [zeros(nu,1); xj; zeros(nz,1); 0];
            wj    = W_scale * w_raw;
            for r = 1:size(Xf.A,1)
                J=Xf.A(r,:)*X1*gj + norm(Lx_cell{r}*wj,2);
                optimize(Constr,-J,ops);
                deltaRay(r) = max(deltaRay(r), value(J));
            end
        else
            for r = 1:size(Xf.A,1)
                J=Xf.A(r,:)*X1*gj;
                optimize(Constr,-J,ops);
                deltaRay(r) = max(deltaRay(r), value(J));
            end
        end
    end
end
deltaRay = max(deltaRay, 0);

%% Build reusable feasibility optimizer: given (gx,ge), find feasible u
x = sdpvar(nx,1);
z = sdpvar(nz,1);

u = sdpvar(nu,1);
g = sdpvar(ell,1);

Constr = [];
Constr = [Constr, W*g == [u;x;z;1]];
Constr = [Constr, Pu.A*u <= Pu.b];

if ~isempty(Lx_cell)
    for r = 1:size(Xf.A,1)
        w_raw = [u; x; z; 1];
        wj    = W_scale * w_raw;
        Constr = [Constr, Xf.A(r,:)*X1*g + norm(Lx_cell{r}*wj,2) + deltaRay(r) <= Xf.b(r)];
    end
else
    Constr = [Constr, Xf.A*X1*g + deltaRay <= Xf.b];
end

ops = sdpsettings('solver','mosek','verbose',0);
verifyOpt = optimizer(Constr, [], ops, {x, z}, {u});

%% Loop over vertex pairs
nVx = size(Xf.V,1);
nVz = size(Pz.V,1);
totalPairs = nVx*nVz;

feasibleCount = 0;
infeasiblePairs = [];  % [j,k]

for j = 1:nVx
    xv = Xf.V(j,:).';
    for k = 1:nVz
        zv = Pz.V(k,:).';

        [u, diag] = verifyOpt{ xv, zv };
        isFeas = (diag == 0)&&(~isnan(u));

        if isFeas
            feasibleCount = feasibleCount + 1;
        else
            infeasiblePairs(end+1,:) = [j,k]; 
        end
    end
end

%% Report
report = struct();
report.isRCIS          = (feasibleCount == totalPairs);
report.nVertexX        = nVx;
report.nVertexZ        = nVz;
report.totalPairs      = totalPairs;
report.feasiblePairs   = feasibleCount;
report.infeasiblePairs = totalPairs - feasibleCount;
report.infeasibleIndex = infeasiblePairs;

fprintf('\n=== RCIS Verification Report ===\n');
fprintf('Pgx vertices: %d, Pgz vertices: %d\n', nVx, nVz);
fprintf('Total pairs: %d | Feasible: %d | Infeasible: %d\n', ...
        totalPairs, feasibleCount, totalPairs-feasibleCount);
fprintf('RCIS satisfied: %s\n', string(report.isRCIS));

if ~report.isRCIS && ~isempty(infeasiblePairs)
    fprintf('First infeasible pair: (j,k)=(%d,%d)\n', infeasiblePairs(1,1), infeasiblePairs(1,2));
end
fprintf('================================\n\n');

end
