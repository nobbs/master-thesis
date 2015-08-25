% This script calculates the inf-sup-constant `\beta_{\mathcal N}` and the bound
% for this constant given in Section 4.2.

% matrices to save the calculated data
bounds = zeros(nt, ns);
exact  = zeros(nt, ns);
cfls   = zeros(nt, ns);

% progress indicator
textprogressbar('calculating: ');

% iterate!
for tdx = 1:nt
  for sdx = 1:ns
    textprogressbar(((tdx - 1) * ns + sdx) / (nt * ns) * 100);

    % set the remaining settings for this iteration
    pd.tgrid = linspace(pd.tspan(1), pd.tspan(2), tvalues(tdx));
    % create the assembly objects
    temporal = TemporalAssemblyNodal(pd, 0);
    spatial  = SpatialAssembly(pd, svalues(sdx));
    solver   = SolverNodal(pd, spatial, temporal);

    % if we have more than one field, then we have to decompose the temporal
    % matrices into parts corresponding to the temporal subintervals. since this
    % is not easily implemented into the temporal assembly objects, we will use
    % some hacky stuff to find the correct decomposition.

    % first, assemble the matrices
    AtE = temporal.stiffnessMatrix();
    MtF = temporal.massMatrix('test');
    % and generate the cells to hold the decompositions
    AtEc = cell(pd.nF, 1);
    MtFc = cell(pd.nF, 1);

    % second, search for the correct time grid indexes
    tpoints = [pd.tspan(1), pd.f, pd.tspan(2)];
    [check, tindex] = ismember(tpoints, pd.tgrid);
    assert(sum(check) == pd.nF + 1);

    % third, iterate over the subintervals and decompose the matrices
    for fdx = 1:pd.nF
      % get the right indexes
      idxf = tindex(fdx);
      idxt = tindex(fdx + 1);
      rekt = idxf:idxt;

      % start with the stiffness matrix. here we have to split the first / last
      % entry on the main diagonal (only if it's not start or end point of the
      % complete time interval!)
      AtEc{fdx} = sparse(size(AtE, 1), size(AtE, 2));
      AtEc{fdx}(rekt, rekt) = AtE(rekt, rekt);
      if fdx ~= 1 && fdx ~= pd.nF
        AtEc{fdx}(rekt(1), rekt(1)) = AtEc{fdx}(rekt(1), rekt(1)) / 2;
        AtEc{fdx}(rekt(end), rekt(end)) = AtEc{fdax}(rekt(end), rekt(end)) / 2;
      elseif fdx ~= 1
        AtEc{fdx}(rekt(1), rekt(1)) = AtEc{fdx}(rekt(1), rekt(1)) / 2;
      elseif fdx ~= pd.nF
        AtEc{fdx}(rekt(end), rekt(end)) = AtEc{fdx}(rekt(end), rekt(end)) / 2;
      end

      % and now the mass matrix of the test space. we simply split it and are
      % done.
      rekt = rekt(1):(rekt(end) - 1);
      MtFc{fdx} = sparse(size(MtF, 1), size(MtF, 2));
      MtFc{fdx}(rekt, rekt) = MtF(rekt, rekt);
    end

    % fourth, assemble space matrices
    Ax = spatial.stiffnessMatrix();
    Hx = spatial.massMatrix();
    Vx = Hx + Ax;
    if useSineExpansion
      Wx = spatial.fieldDependentSine();
    else
      Wx = spatial.fieldDependentFourier();
    end
    % and the spatial field dependent part
    Ac = cell(pd.nF, 1);
    for fdx = 1:pd.nF
      Ac{fdx} = pd.laplacian * Ax + pd.offset * Hx;
      for cdx = 1:pd.nC
        pos     = (fdx - 1) * pd.nC + cdx;
        Ac{fdx} = Ac{fdx} + fieldcoeffs(pos) * Wx{cdx};
      end
    end

    % fifth, assemble the norms.
    % the plus-Norm is straight forward
    Mplus = kron(MtFc{1}, Ac{1});
    for fdx = 2:pd.nF
      Mplus = Mplus + kron(MtFc{fdx}, Ac{fdx});
    end
    % for the minus-Norm we have to assure, that we check only elements where
    % the norm is not zero (since it's part of the denominator). since we
    % evaluate the minus-Norm of the time derivative of an ansatz function, we
    % have to check, for which functions this time derivative is zero. This is
    % the case iff the ansatz function is constant. In the end, we can obtain
    % the wanted result by calculating all eigenvalues and eigenvectors of AtE,
    % ignoring the eigenvector of the zero eigenvalue and multiplying AtE with
    % the remaining eigenvectors. This removes the zero eigenvalue and allows us
    % to calculate the inf-sup-constant as a generalized eigenvalue problem.
    [V, d] = eigs(AtE, size(AtE, 1) - 1);
    seye   = speye(size(Ac{1}));
    Mminus = kron(V.' * AtEc{1} * V, Ac{1} \ seye);
    for fdx = 2:pd.nF
      Mminus = Mminus + kron(V.' * AtEc{fdx} * V, Ac{fdx} \ seye);
    end
    Mminus = (Mminus + Mminus.') / 2;

    % sixth, assemble the needed matrix to compute the supremizer
    CtFE = temporal.halfStiffnessMatrix();
    Cl   = kron(CtFE * V, Hx);

    % seventh, finally, construct the generalized eigenvalue problem and solve it
    Yl         = Cl.' * (Mplus \ Cl);
    Yl         = (Yl + Yl.') / 2;
    [mi, ~, ~] = computeMinMaxEv(Yl, Mminus);
    betapm     = sqrt(mi);

    % calculate the CFL number
    maxTk      = max(diff(pd.tgrid));
    ma         = eigs(Vx, Hx.' * (Vx \ Hx), 1, 'lm');
    cfl        = maxTk * sqrt(ma);

    % calculate the bound for the inf-sup-constant
    infsup_bound = min(1, betapm) * min(1, 1 / cfl);

    % and now calculate the real inf-sup-constant
    Lhs          = solver.spacetimeSystemMatrix(fieldcoeffs);
    Ynorm        = solver.TeNorm;
    Xnorm        = solver.TrNorm;
    supNorm      = Lhs.' * (Ynorm \ Lhs);
    [mi, ~, ~]   = computeMinMaxEv(supNorm, Xnorm);
    infsup_exact = sqrt(mi);

    % save the data
    bounds(tdx, sdx) = infsup_bound;
    exact(tdx, sdx)  = infsup_exact;
    cfls(tdx, sdx)   = cfl;
  end
end

textprogressbar(' done!');
