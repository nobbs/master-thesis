classdef SCM < handle
  % Successive constraint method.
  %
  % Allows the computation of a lower bound for the inf-sup-condition of a
  % variational problem.
  %
  % @todo refactor
  % @todo describe
  % @todo holy macaroni!

  properties
    % Total number of terms in the affine representation of the supremizers @type integer
    nQt;

    % Max error in each iteration @type vector
    errors;

    offC;% #
    offD; % #
    offThetaC;
    offThetaD; % #
    offUb;% #
    offLb;% #
    offalphaC;% #
    offystarC; % #

    % Struct that holds all the data needed from the offline stage needed for
    % the online stage. @type stuct
    offData;
  end

  properties(Access = 'protected')
    % Reference to the reduced basis solver object from which this object was
    % created @type RBMSolverAbstract;
    rbm;

    % Holds the affine representation of the Y-norm of the supremizers @type
    % cellarray
    Aq;
  end


  methods

    function obj = SCM(rbm)
      obj.rbm = rbm;
    end

    function prepare(obj)
      % number of terms in the affine representation
      obj.nQt = obj.rbm.nQb * (obj.rbm.nQb + 1) / 2;

      Aq = cell(obj.nQt, 1);
      Yn = obj.rbm.solver.TeNorm;
      Bq = obj.rbm.solver.Lhs;

      idx = 1;
      for pdx = 1:obj.rbm.nQb
        Aq{idx} = Bq{pdx}.' * (Yn \ Bq{pdx});
        idx = idx + 1;
      end

      for pdx = 1:obj.rbm.nQb
        for qdx = (pdx + 1):obj.rbm.nQb
          Aq{idx} = (Bq{pdx} + Bq{qdx}).' * (Yn \ (Bq{pdx} + Bq{qdx}));
          idx = idx + 1;
        end
      end

      obj.Aq = Aq;
    end

    function newparam = mapParam(obj, param)
      newparam = zeros(obj.nQt, 1);

      % pad the param with 1 (field independent part)
      param = [1; shiftdim(param)];

      idx = 1;
      for pdx = 1:obj.rbm.nQb
        tmp = param(pdx)^2;
        for jdx = 1:obj.rbm.nQb
          if jdx ~= pdx
            tmp = tmp - param(pdx) * param(jdx);
          end
        end
        newparam(idx) = tmp;
        idx = idx + 1;
      end

      for pdx = 1:obj.rbm.nQb
        for qdx = (pdx + 1):obj.rbm.nQb
          newparam(idx) = param(pdx) * param(qdx);
          idx = idx + 1;
        end
      end
    end

    function M = assemble(obj, param)
      param = obj.mapParam(param);
      M = sparse(obj.rbm.nTrialTruth, obj.rbm.nTrialTruth);
      for idx = 1:obj.nQt
        M = M + param(idx) * obj.Aq{idx};
      end
    end

    function [mi, ma] = computeMinMaxEV(obj, A, B)
      lm = obj.computeEV(A, B, 0);
      if lm >= 0
        ma = lm;
        mi = obj.computeEV(A, B, lm);
      else
        mi = lm;
        ma = obj.computeEV(A, B, lm);
      end
    end

    function [lb, ub] = bounds(obj)
      % @todo Not yet implemented!

      Xn = obj.rbm.solver.TrNorm;

      % Berechnung der Bounds für die Variablen
      lb = zeros(obj.nQt, 1);
      ub = zeros(obj.nQt, 1);
      for idx = 1:obj.nQt
        [mi, ma] = obj.computeMinMaxEV(obj.Aq{idx}, Xn);
        ub(idx) = ma;
        lb(idx) = mi;
      end
    end

    function offlinePhase(obj)
      % @todo Not yet implemented!
      % error('Not yet implemented!');

      Xn = obj.rbm.solver.TrNorm;


      % Calculate bounds for the variables
      [lb, ub] = obj.bounds()
      obj.offLb = lb;
      obj.offUb = ub;

      % Set Xi definieren (sollte groß sein)
      obj.offD = [-3:0.1:3];
      obj.offThetaD = zeros(obj.nQt, size(obj.offD, 2));
      for idx = 1:size(obj.offD, 2)
        obj.offThetaD(:, idx) = obj.mapParam(obj.offD(:, idx));
      end

      % Training set for C_K
      CkTrain = [-3:0.4:3];

      % keyboard

      % preperation for the greedy algorithm
      isDone = false;
      alphaC = [];
      ystarC = [];
      C = [];
      thetaC = [];
      max_err = [];
      exflag = 0;

      % break conditions
      tolerance = 1e-10;
      maxiter = 30;

      curparam = CkTrain(:, 1);
      while ~isDone


        % assemble the system for the chosen parameter and calculate the
        % smallest eigenvalue and it's eigenvector
        M = obj.assemble(curparam);
        [mi, mivec] = obj.computeEV(M, Xn, 0);
        % if >= then we have the largest ev and do it again!!!!!
        if mi >= 0
          [mi, mivec] = obj.computeEV(M, Xn, mi);
        end

        % calculate ystar
        ystar = zeros(obj.nQt, 1);
        for idx = 1:obj.nQt
          ystar(idx) = mivec.' * obj.Aq{idx} * mivec;
        end

        C = [C, curparam];
        alphaC = [alphaC; mi];
        ystarC = [ystarC, ystar];
        thetaC = [thetaC, obj.mapParam(curparam)];

        % save for online use
        obj.offC = C;
        obj.offalphaC = alphaC;
        obj.offystarC = ystarC;
        obj.offThetaC = thetaC;


        % @todo feasability check for the linprog

        % call online solve
        quots = zeros(size(CkTrain, 2), 1);
        for idx = 1:size(CkTrain, 2)
          chkparam = CkTrain(:, idx);
          [l, u] = obj.onlineSolve(chkparam, 50, 200)
          quots(idx) = (u - l) / u;
        end

        % look for the max quot
        [maxquot, maxdx] = max(quots);
        max_err = [max_err, maxquot];
        obj.errors = max_err;

        % abbruchbedingungen
        if maxquot < tolerance
          exflag = 1;
          isDone = true;
        elseif size(C, 2) > maxiter
          exflag = 2;
          isDone = true;
        elseif isempty(CkTrain) == 1
          exflag = 3;
          isDone = true;
        else
          isDone = false
          curparam = CkTrain(:, maxdx);
          CkTrain(:, maxdx) = [];
        end
      end

      keyboard
    end

    function [lm, lmvec] = computeEV(obj, A, B, shift)
      % @todo describe whats going on...
      % @todo better error handling
      defshift = eps;
      warnerr = warning('error', 'MATLAB:eigs:NoEigsConverged');

      opts = struct();
      opts.tol = eps;

      try
        [lmvec, lm] = eigs(A - (shift + defshift) * B, B, 1, 'lm');
        lm = lm  + (shift + defshift);
      catch err
        if strcmp(err, 'MATLAB:eigs:NoEigsConverged') || strcmp(err.identifier, 'MATLAB:eigs:ARPACKroutineErrorMinus14')
          hasEV = false;
          while ~hasEV
            try
              opts.tol = opts.tol * 100;
              opts.p = 50;
              [lmvec, lm] = eigs(A - (shift + defshift) * B, B, 1, 'lm', opts);
              lm = lm  + (shift + defshift);
              hasEV = true;
            catch err2
              hasEV = false;
              if strcmp(err, 'MATLAB:eigs:NoEigsConverged') || strcmp(err2.identifier, 'MATLAB:eigs:ARPACKroutineErrorMinus14')
                opts.tol = opts.tol * 100;
              else
                rethrow(err2);
              end
            end
          end
        else
          rethrow(err);
        end
      end

      % Null is null!
      if abs(lm) < sqrt(eps)
        lm = 0;
      end
    end

    function [lb, ub] = onlineSolve(obj, param, Malpha, Mplus)
      % @todo Not yet implemented!

      [~, Idxalpha] = obj.getNeighbors(Malpha, param, obj.offC);
      [~, Idxplus] = obj.getNeighbors(Mplus, param, obj.offD);

      % @todo should perform a feasability check here

      f = obj.mapParam(param);
      % keyboard
      b = [obj.offUb; -obj.offLb; -obj.offalphaC(Idxalpha); zeros(size(Idxplus, 2), 1)];
      A = [eye(size(obj.offUb, 1)); - eye(size(obj.offLb, 1)); -obj.offThetaC(:, Idxalpha).'; -obj.offThetaD(:, Idxplus).'];

      [~, fval, flag, ~] = linprog(f, A, b, [],  [],  [],  [], obj.offUb, struct('LargeScale','off', 'Display', 'off'));

      if flag < 0
        warning('somethings fishy!');
        keyboard;
      end

      lb = sqrt(fval);
      ub = sqrt(min(obj.offystarC.' * f));
    end

    function [Pm, Idx] = getNeighbors(obj, M, param, C)
      % @todo document

      Pm = [];
      Idx = [];

      if M ~= 0 && isempty(C) ~= 1
        % number of params in c
        nc = size(C, 2);
        % calculate euclidian distance from param
        Param = repmat(param, 1, nc);
        [~, order] = sort(sqrt(sum((Param - C).^2, 1)));
        Idx = order(1:min(nc, M));
        Pm = C(:, Idx);
      end
    end

  end

end
