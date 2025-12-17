function L_cell = EstimateUncertaintyBoundsFromData(A_all, U0, X0, Z0, X1, U0_full, X0_full, Z0_full,  X1_full, W_scale, alpha)
% BUILDUNCERTAINTYBOUNDSFROMDATA
% For each row of A_all, build an SDP-based factor L{i} such that
% |A_all(i,:) * (X1_full - X1_hat)(:,j)| <= ||L{i} * (W_scale * w_j)||_2.

tic
[nu, ell] = size(U0); nx = size(X0, 1); nz = size(Z0, 1); nW = nu + nx + nz + 1; N = size(U0_full, 2); nRow = size(A_all, 1);

% Basis and pseudoinverse
W_basis = [U0; X0; Z0; ones(1, ell)];    % (nW × ell)
W_pinv  = pinv(W_basis);                 % (ell × nW)

% Full data and scaled version
W_full   = [U0_full; X0_full; Z0_full; ones(1, N)];   % (nW × N)
W_scaled = W_scale * W_full;                          % (nW × N)

% All g_j, predicted X1_hat, and residuals
G_full = W_pinv * W_full;              % (ell × N)
X1_hat = X1 * G_full;                  % (nx × N)
E_all  = A_all * (X1_full - X1_hat);   % (nRow × N)

% SDP setup
opts  = sdpsettings('solver','mosek','verbose',0,'cachesolvers',1);
R_row = cell(nRow, 1);
I_reg = 1e-16 * eye(nW);

for i = 1:nRow
    Qi     = sdpvar(nW, nW, 'symmetric');
    Constr = (Qi >= I_reg);
    e_i    = E_all(i, :);          % 1 × N
    idx    = 1:N;                  % use all samples (most conservative)
    for j = idx
        wj     = W_scaled(:, j);   % scaled regressor w_j
        Constr = [Constr, e_i(j)^2 <= wj.' * Qi * wj];
    end
    sol = optimize(Constr, trace(Qi), opts);
    if sol.problem ~= 0
        warning('Row %d: solver issue: %s', i, sol.info);
    end
    Qi_val = alpha * value(Qi);
    Qi_sym = 0.5 * (Qi_val + Qi_val.');
    [Rtmp, p] = chol(Qi_sym + 1e-12 * eye(nW), 'lower');
    if p == 0
        R_row{i} = Rtmp;                 % R'*R ≈ Qi
    else
        [Uq, Sq, ~] = svd(Qi_sym);
        R_row{i}    = sqrt(Sq) * Uq.';   % ||R w||₂² = w'Qi w
    end
end

L_cell = R_row;
toc
end
