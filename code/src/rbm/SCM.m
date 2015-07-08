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
    % Total number of addends in the affine representation of the supremizers
    % @type integer
    nQt;

    % Total number of addends in the affine representation of the underlying
    % variational problem @type integer
    nQb;

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
    affineT;
  end


  methods

    function obj = SCM(rbm)
      % Constructor for this class.
      %
      % Parameters:
      %   rbm: Reference to the reduced basis solver that created this scm
      %     object. @type RBMSolverAbstract

      % first set the reference
      obj.rbm = rbm;

      % and now prepare the scm offline stage
      obj.prepare();
    end

  end

  methods

    % Main methods of the successive constraint method algorithm

    function offlinePhase(obj)
      % @todo Not yet implemented!
      % error('Not yet implemented!');

      Xn = obj.rbm.solver.TrNorm;


      % Calculate bounds for the variables
      [lb, ub] = obj.bounds()
      obj.offLb = lb;
      obj.offUb = ub;

      % Set Xi definieren (sollte groß sein)
      % obj.offD = [-3:0.1:3];
      obj.offD = rand(1, 10000) * 10 - 5;
      obj.offThetaD = zeros(obj.nQt, size(obj.offD, 2));
      for idx = 1:size(obj.offD, 2)
        obj.offThetaD(:, idx) = obj.mapParam(obj.offD(:, idx));
      end

      % Training set for C_K
      % CkTrain = [-3:0.4:3];
      CkTrain = rand(1, 100) * 10 - 5;

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
      mudx = 1;
      while ~isDone
        mudx

        % assemble the system for the chosen parameter and calculate the
        % smallest eigenvalue and it's eigenvector
        M = obj.assembleAffineT(curparam);
        [mi, mivec] = obj.computeEV(M, Xn, 0);
        % if >= then we have the largest ev and do it again!!!!!
        if mi >= 0
          [mi, mivec] = obj.computeEV(M, Xn, mi);
        end

        % calculate ystar
        ystar = zeros(obj.nQt, 1);
        for idx = 1:obj.nQt
          ystar(idx) = mivec.' * obj.affineT{idx} * mivec;
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

        mudx = mudx + 1;
      end

      keyboard
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

  methods

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
        [mi, ma] = obj.computeMinMaxEV(obj.affineT{idx}, Xn);
        ub(idx) = ma;
        lb(idx) = mi;
      end
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

  end

  methods(Access = 'private')

    % Helper methods

    function newparam = mapParam(obj, param)
      % Remaps the parameter for the affine representation.
      %
      % This converts a parameter from the variational problem to a modified
      % parameter needed for the final assembly of the Y-norm of the supremizer.
      % For a motivation of this, see my master thesis.
      %
      % @todo Reference!
      %
      % Parameters:
      %   param: parameter from the variational problem @type vector
      %
      % Return values:
      %   newparam: parameter for the supremizer @type vector

      % create the needed vector
      newparam = zeros(obj.nQt, 1);

      % pad the parameter from the variational problem with 1 for the field
      % independent part
      param = [1; shiftdim(param)];

      % and now we handle the assembly of the new parameter the same way we
      % handled the assembly of the addends in SCM@prepare.

      idx = 1;
      % first kind is handled first
      for qdx = 1:obj.rbm.nQb
        val = param(qdx)^2;
        for jdx = 1:obj.rbm.nQb
          if jdx ~= qdx
            val = val - param(qdx) * param(jdx);
          end
        end
        newparam(idx) = val;
        idx = idx + 1;
      end

      % and now the simpler second kind
      for qdx = 1:obj.rbm.nQb
        for pdx = (qdx + 1):obj.rbm.nQb
          newparam(idx) = param(qdx) * param(pdx);
          idx = idx + 1;
        end
      end
    end

    function T = assembleAffineT(obj, param)
      % Assemble the supremizer for the given parameter.
      %
      % Parameters:
      %   param: parameter from the variational problem @type vector
      %
      % Return values:
      %   T: supremizer for the given parameter @type matrix

      % remap the parameter
      param = obj.mapParam(param);

      % and sum up the addends of the affine representation
      T = sparse(obj.rbm.nTrialTruth, obj.rbm.nTrialTruth);
      for idx = 1:obj.nQt
        T = T + param(idx) * obj.affineT{idx};
      end
    end

    function prepare(obj)
      % Prepare the SCM object for the offline phase.
      %
      % This mainly consists of the creation of the addends for the affine
      % representation of the Y-norm of the supremizers.

      % get the needed matrices from the truth solver (through the rbm solver)
      Ry = obj.rbm.solver.TeNorm;
      Bq = obj.rbm.solver.Lhs;

      % set the total number of needed addends for the affine representation
      obj.nQt = obj.rbm.nQb * (obj.rbm.nQb + 1) / 2;

      % create the needed structure
      obj.affineT = cell(obj.nQt, 1);

      % iterate over the different addends of the affine representation. as we
      % are using an alternative expansion of the field independent operators to
      % ensure that they are all symmetric and positiv semi-definite, we have to
      % deal with the two following types of addends:
      %  1. <T_q u, T_q u>_Y for q = 1..Qb,
      %  2. <T_q u + T_p u, T_q u + T_p u>_Y for q = 1..Qb and p = (q+1)..Qb.

      % index for the affineT structure (a mapping of the index pair (q, p) to a
      % useful 1d index would be nice, but... meh).
      idx = 1;

      % assemble the affine addends of the first kind
      for qdx = 1:obj.rbm.nQb
        obj.affineT{idx} = Bq{qdx}.' * (Ry \ Bq{qdx});
        idx = idx + 1;
      end

      % and now assemble the affine addends of the second kind
      for qdx = 1:obj.rbm.nQb
        for pdx = (qdx + 1):obj.rbm.nQb
          obj.affineT{idx} = (Bq{qdx} + Bq{pdx}).' * (Ry \ (Bq{qdx} + Bq{pdx}));
          idx = idx + 1;
        end
      end
    end

  end

end
