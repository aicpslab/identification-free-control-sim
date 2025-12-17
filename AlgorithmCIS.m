function Xf = AlgorithmCIS(U0,X0,Z0,X1,Pu,Px,Pz)
% Data-driven construction of a nominal control-invariant set Xf
% Stable/safe version (H-rep-centric; minimal V usage)

    [ell,nu,nx,nz] = deal(size(U0,2),size(U0,1),size(X0,1),size(Z0,1));
    pUXZ1 = pinv([U0; X0; Z0; ones(1,ell)]);

    canonH = @(P) minHRep(normalize(P));
    canonV = @(P) minVRep(P);

    function P = ensureH(P)
        P = canonH(P);
        if isempty(P.A)
            P = canonH(computeHRep(canonV(P)));
        end
    end

    % Initial set in x-space (H as truth)
    Xi = ensureH(Px);

    % g-space polyhedra (build in H, then compute V once where needed)
    Pgu = Polyhedron('V',(pUXZ1*[Pu.V';zeros(nx+nz+1,size(Pu.V,1))])');
    Pgxi = Polyhedron('V',(pUXZ1*[zeros(nu,size(Xi.V,1));Xi.V';zeros(nz+1,size(Xi.V,1))])','R',(pUXZ1*[zeros(nu,size(Xi.R,1));Xi.R';zeros(nz+1,size(Xi.R,1))])');
    Pgz = Polyhedron('V',(pUXZ1*[zeros(nu+nx,size(Pz.V,1));Pz.V';ones(1,size(Pz.V,1))])');


    Pgu  = ensureH(Pgu);  Pgu  = canonV(Pgu);
    Pgxi = ensureH(Pgxi); Pgxi = canonV(Pgxi);
    Pgz  = ensureH(Pgz);  Pgz  = canonV(Pgz);

    % One-step images under X1 (construct via vertices/rays, then fix to H once)
    iTu = Polyhedron('V', -Pgu.V*X1');
    Te  = Polyhedron('V',  Pgz.V*X1');
    if ~isempty(Pgxi.R)
        Txi = Polyhedron('V', Pgxi.V*X1', 'R', Pgxi.R*X1');
    else
        Txi = Polyhedron('V', Pgxi.V*X1');
    end
    iTu = ensureH(iTu);
    Te  = ensureH(Te);
    Txi = ensureH(Txi);

    % Candidate post-set
    Pi = ensureH(Txi + Te);
    Si = ensureH(Xi + iTu);
    i  = 0;

    options = sdpsettings('solver','mosek','verbose',0);

    while ~(Si >= Pi) && i <= 50

        PiTilde = ensureH(Si & Pi);

        % rho via support in Te along normals of PiTilde
        rho = zeros(size(PiTilde.A,1),1);
        x   = sdpvar(nx,1);

        Te = ensureH(Te);

        for j = 1:size(PiTilde.A,1)
            Constr = (Te.A*x <= Te.b);
            if ~isempty(Te.Ae)
                Constr = [Constr, Te.Ae*x == Te.be];
            end
            sol = optimize(Constr, -PiTilde.A(j,:)*x, options);
            if sol.problem ~= 0
                rho(j) = inf;
            else
                rho(j) = value(PiTilde.A(j,:)*x);
            end
        end

        % Update Txi and Pi (H only)
        Txi = Polyhedron('A', PiTilde.A, 'b', PiTilde.b - rho);
        Txi = ensureH(Txi);
        Pi  = ensureH(Txi + Te);

        % Update Pgxi in g-space (H only; avoid V in the loop)
        Pgxi_new = ensureH(Polyhedron('A', Txi.A*X1, 'b', Txi.b));
        Pgxi     = ensureH(Pgxi);    % ensure Pgxi has H after earlier V usage
        Pgxi     = ensureH(Pgxi_new & Pgxi);

        % Map back to x-space
        Xi = ensureH(Pgxi.affineMap(X0));

        i  = i + 1;
        Si = ensureH(Xi + iTu);
    end

    fprintf('Iteration: %d\n', i);
    Xf = ensureH(Xi).computeVRep
end
