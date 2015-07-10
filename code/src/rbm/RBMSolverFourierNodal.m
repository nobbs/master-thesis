classdef RBMSolverFourierNodal < RBMSolverAbstract
  % Solver based on the reduced basis method.

  properties
    verbose = true;

    % @todo Holds the affine decompositon of somethign @type cell
    affineTestshotBase;
  end

  properties %(Access = 'protected')
    G;
    scm;
    isFreshStart;
    errors;
    params;
  end

  properties(Dependent)
    nQf;
    nQb;
    nQr;
    nTrialTruth;
    nTestTruth;
    nTrialRB;
    nTestRB;
  end

  methods
    function val = get.nQf(~)
      val = 1;
    end

    function val = get.nQb(obj)
      val = 1 + obj.pd.nP;
    end

    function val = get.nQr(obj)
      val = obj.nQf + obj.nTrialRB * obj.nQb;
    end

    function val = get.nTrialTruth(obj)
      val = obj.solver.nTrialDim;
    end

    function val = get.nTestTruth(obj)
      val = obj.solver.nTestDim;
    end

    function val = get.nTrialRB(obj)
      val = size(obj.trialSnapshots, 2);
    end

    function val = get.nTestRB(obj)
      val = obj.nTrialRB;
    end
  end

  methods

    function obj = RBMSolverFourierNodal(problem)
      % Constructor for this class.
      %
      % Parameters:
      %   problem: reference to a problem data object. @type ProblemData

      % first call the constructor of the superclass
      obj@RBMSolverAbstract(problem);

      % save the problem data object reference
      % obj.pd = problem;
    end

    function [rbsolvec, errest] = onlineSolve(obj, param)
      % Solve online for a given Parameter with the reduced basis solver.
      %
      % Parameters:
      %   param: parameter to use. @type colvec
      %
      % Return values:
      %   rbsolvec: solution vector for the rb truth trial space @type colvec
      %   errest: reduced basis error estimate @type double

      % construct the reduced basis system
      testspace = obj.constructTestSpace(param);
      Lhs       = testspace.' * obj.solver.spacetimeSystemMatrix(param) * ...
                    obj.trialSnapshots;
      Rhs       = testspace.' * obj.solver.spacetimeLoadVector;
      % now solve it
      rbsolvec   = Lhs \ Rhs;
      % compute the error estimate
      errest     = obj.estimateError(param, rbsolvec);
    end

    function errest = estimateError(obj, param, rbsolvec)
      % Estimate the error of a reduced basis solution.
      %
      % Parameters:
      %   param: parameter to use. @type colvec
      %   rbsolvec: solution vector of the reduced system. @type colvec
      %
      % Return values:
      %   errest: estimated error @type double

      % @todo move this somewhere else
      Malpha = 10;
      Mplus = 100;

      % calculate the residual norm
      res = obj.residual(rbsolvec, param);
      % calculate the bounds of the inf-sup-constant
      [lb, ~] = obj.scm.onlineQuery(param, Malpha, Mplus);

      % finally compute the error estimate
      errest = res / lb;
    end

    function resetTraining(~)
    end

    function prepare(obj)
      % Prepare the reduced basis method solver.
      %
      % This mainly consists of forwarding needed values to underlying classes
      % and setting up the galerkin solver.
      %
      % @todo refactor this pile of junk...


      obj.solver = SolverFourierNodal(obj.pd, 10, 0);

      obj.trialSnapshots = sparse(obj.solver.nTrialDim, 1);
      obj.isFreshStart = true;
      obj.affineTestshotBase = cell(0);
      obj.errors = [];
      obj.params = [];
    end

    function testshots = constructTestSpace(obj, param)
      % @todo describe it

      testshots = sparse(obj.nTestTruth, obj.nTrialRB);

      % pad with 1 for the field independent part
      param = [1; shiftdim(param)];

      for ndx = 1:obj.nTrialRB
        testshots(:, ndx) = obj.affineTestshotBase{ndx} * param;
      end
    end

    function offlineStage(obj, paramTrain)
      % Offline stage of the reduced basis method.
      %
      % This performs a greedy style algorithm to "train" the reduced basis
      % system. The greedy loop works as follows:
      % - in each cycle we compute error estimates for all the remaining
      %   parameters in our training set and choose the one with the largest
      %   error.
      % - if the largest error is above a given tolerance, then we use the truth
      %   solver to get a solution for this parameter and add this solution to
      %   our reduced basis trial space basis and go to the next cycle.
      % - else we abort and are done with the offline stage.
      %
      % Parameters:
      %   paramTrain: parameter training set. @type matrix

      if obj.verbose
        fprintf('# RBM: starting offline stage.\n');
      end

      % @todo move it
      ptest = rand(1, 10000) * 6 - 3;

      % first we have to set up the successive constraint method and start it's
      % offline stage, as we are going to need a fast way to calculate the inf-
      % sup-constant for lots of parameters
      obj.scm = SCM(obj);
      % let's start the offline stage of the scm. this could take a while...
      obj.scm.offlineStage(paramTrain, ptest);

      % now we are ready to start the greedy loop! first some more preparation
      isDone = false;
      exflag = 0;
      tol = 1e-6;
      
      maxerr = 0;

      if obj.verbose
        fprintf('# RBM: starting greedy loop ');
      end

      % select the first parameter by random
      curIdx = randi(size(paramTrain, 2), 1, 1);
      curTrain = paramTrain(:, curIdx);
      paramTrain(:, curIdx) = [];
      obj.params(:, 1) = curTrain;

      % and now lets move!
      while ~isDone
        if obj.verbose
          fprintf('.');
          maxerr
        end

        % debugging
%         keyboard;

        % first we compute the truth solution for the current parameter
        solvec = obj.solver.solve(curTrain);

        % now we use a gram-schmidt-orthonormalization
        z = solvec;
        if ~obj.isFreshStart
          z = solvec - obj.trialSnapshots * (solvec.' * obj.solver.TrNorm * obj.trialSnapshots).';
          z = z / (sqrt(z' * obj.solver.TrNorm * z));
          obj.trialSnapshots(:, end + 1) = z;
        else
          z = z / (sqrt(z' * obj.solver.TrNorm * z));
          obj.trialSnapshots(:, 1) = z;
          obj.isFreshStart = false;
        end
        
        % and compute the needed affine decomposition for the on-demand assembly
        % of the reduced basis test space

        % iterate over the space time matrix parts
        Zqn = zeros(obj.nTestTruth, obj.nQb);
        for pdx = 1:obj.nQb
          Zqn(:, pdx) = obj.solver.TeNorm \ (obj.solver.Lhs{pdx} * solvec);
        end
        obj.affineTestshotBase{end + 1} = Zqn;

        % next step: compute all the error estimates, so let's prepare the
        % computation of the residual
        obj.prepareResidual;
        % and now lets iterate
        errests = zeros(size(paramTrain, 2), 1);
        for idx = 1:size(paramTrain, 2)
          [rbsolvec, errest] = obj.onlineSolve(paramTrain(:, idx));
          errests(idx) = errest;
        end

        % get the parameter with the largest error estimate
        [maxerr, maxdx] = max(errests);
        obj.errors(end + 1) = maxerr;

        % check the breaking conditions
        % @todo add more of 'em
        if maxerr < tol
          % if the maximum error is smaller than the given tolerance, then we
          % are done here!
          exflag = 1;
          isDone = true;
        else
          % looks like we have to do more cycles, so select the right parameter
          curTrain = paramTrain(:, maxdx);
          paramTrain(:, maxdx) = [];
          obj.params(:, end + 1) = curTrain;
        end
      end

      if obj.verbose
        fprintf(' done!\n');
      end
    end

    function prepareResidual(obj)
      % Set up the needed matrix for the computation of the residual.
      %
      % This assembles the solution and field independent part of the
      % computation of the Y'-norm of the residual through it's riesz
      % representation.
      %
      % @todo describe it more!

      % create the matrix
      H = zeros(obj.nTestTruth, obj.nQr);

      % the first column is the rhs of the truth system
      H(:, 1) = obj.solver.spacetimeLoadVector;;

      % the remaining columns are the products of the affine composition of the
      % truth system matrix and the chosen trial truth space snapshots
      for pdx = 1:obj.nQb
        H(:, (pdx + 1):obj.nQb:end) = - obj.solver.Lhs{pdx} * obj.trialSnapshots;
      end

      % and now we multiply that all together and the preparation is done
      obj.G = H.' * (obj.solver.TeNorm.' \ H);
    end

    function val = residual(obj, uN, param)
      % Calculate the residual of the rb solution for the given parameter.
      %
      % This is done by evaluation of the y-norm of the riesz representation of
      % the residual.
      %
      % Parameters:
      %   uN: rb trial space coefficient vector of the rb solution @type vector
      %   param: parameters @type matrix
      %
      % Return values:
      %   val: residual of the rb solution

      % create the vector for the construction of the field dependent needed for
      % the computation
      Eps    = zeros(obj.nQr, 1);
      % the first component is the parameter of the rhs. (we don't have one!)
      Eps(1) = 1;

      % the next components (spaced over the vector) is the parameter of the
      % field independent part of the system matrix
      Eps(2:obj.nQb:end) = 1 * uN;

      % pad the parameter with 1 for field independent part
      param = [1; shiftdim(param)];

      % and the remaining parts are the field dependent parts, also spaced over
      % the whole vector
      for pdx = 1:obj.nQb
        Eps((pdx + 1):obj.nQb:end) = param(pdx) * uN;
      end

      % and now calculate the y-norm of the riesz representation
      val  = sqrt(abs(Eps.' * obj.G * Eps));
    end

    function solval = evaluateSolutionRb(obj, rbsolvec)
      % Evaluate a solution of the reduced basis system.
      %
      % Parameters:
      %   rbsolvec: solution vector of the rb system @type vector
      %
      % Return values:
      %   solval: values of the solution on the grid points @type matrix

      % forward the resulting vector of the truth trial space to the
      % corresponding evaluation method
      solval = obj.evaluateSolutionTruth(obj.trialSnapshots * rbsolvec);
    end

    function solval = evaluateSolutionTruth(obj, gsolvec)
      % Evaluate a solution of the truth solver.
      %
      % Parameters:
      %   gsolvec: solution vector of the truth trial space @type vector
      %
      % Return values:
      %   solval: values of the solution on the grid points @type matrix

      % Just forward it to the truth solver method
      solval = obj.solver.evaluateSolution(gsolvec);
    end

  end

  % collection of maybe necessary stuff

  methods%(Access = 'protected')
    function [beta, gamma] = calcDiscreteInfSupAndContinuityRB(obj, param)
      % Calculate the discrete inf-sup-constant of the rb system for the given
      % parameter.
      %
      % Parameters:
      %   param: field parameters @type matrix
      %
      % Return values:
      %   beta: discrete inf-sup-constant @type double
      %   gamma: discrete continuity constant @type double
      %
      % @todo eventuell affine zerlegbarkeit ausnutzen!

      testshots = obj.constructTestSpace(param);

      % first we have to construct the rb system matrix for this parameter and
      % the y- and x-norm matrices.
      Lhs   = full(testshots.' * obj.solver.spacetimeSystemMatrix(param) * obj.trialSnapshots);
      Ynorm = full(testshots.' * obj.solver.TeNorm * testshots);
      Xnorm = full(obj.trialSnapshots.' * obj.solver.TrNorm * obj.trialSnapshots);

      % now we solve the generalized eigenvalue problem, the inf-sup and continuity
      % constants are then the square roots of the smallest respectively largest
      % generalized eigenvalue
      supNorm = Lhs.' * (Ynorm \ Lhs);
      ev      = eig(supNorm, Xnorm);
      beta    = sqrt(min(ev));
      gamma   = sqrt(max(ev));
      cond_ = gamma / beta
    end

    function [beta, gamma] = calcDiscreteInfSupAndContinuityTruth(obj, param)
      % Calculate the discrete inf-sup-constant of the truth system for the
      % given parameter.
      %
      % Parameters:
      %   param: field parameters @type matrix
      %
      % Return values:
      %   beta: discrete inf-sup-constant @type double
      %   gamma: discrete continuity constant @type double
      %
      % @todo eventuell affine zerlegbarkeit ausnutzen!

      % get all the needed matrices from the truth solver
      Lhs   = obj.solver.spacetimeSystemMatrix(param);
      Ynorm = obj.solver.TeNorm;
      Xnorm = obj.solver.TrNorm;

      % solve the generalized eigenvalue problem, the inf-sup and continuity
      % constants are the square root of the smallest respectively largest
      % generalized eigenvalue.
      supNorm = Lhs.' * (Ynorm \ Lhs);
      beta    = sqrt(eigs(supNorm, Xnorm, 1, 'sm'));
      gamma   = sqrt(eigs(supNorm, Xnorm, 1, 'lm'));
      cond_ = gamma / beta
    end
  end

end
