classdef RBM < handle
  % Reduced basis method.
  %
  % This class implements the reduced basis method solver. It theoretically
  % supports different galerkin-based truth solvers as long as these implement
  % all the needed properties and methods defined in SolverAbstract.

  properties
    % Toggle verbosity of the solver. @type logical
    verbose = true;

    % Reference to the given problem data object. @type ProblemData
    pd;
    % Reference to the galerkin truth solver @type SolverAbstract
    solver;
    % Reference to the scm object @type SCM
    scm;

    % Holds the largest error estimate obtained by the greedy procedure @type vector
    errors = [];
    % Holds the selected parameters @type matrix
    params= [];

    % Number of reduced basis trial basis vectors @type integer
    nTrialRB = 0;
    % Number of reduced basis test basis vectors @type integer
    nTestRB = 0;
  end

  properties%(Access = 'protected')
    % Holds the truth solution vectors which form the basis of the reduced basis
    % trial space. @type matrix
    trialSnapshots;
    % @todo Holds the affine decompositon of somethign @type cell
    affineTestshotBase;
    % Needed for the computation of the residual @type matrix
    G;
    % Toggle, whether the offline data was reset recently @type logical
    % isFreshStart = true;
  end

  properties(Dependent)
    nQf;
    nQb;
    nQr;
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
  end

  methods
    function obj = RBM(problem, truthsolver, scm)
      % Default constructor.
      %
      % Parameters:
      %   problem: reference to a problem data object. @type ProblemData
      %   truthsolver: reference to the truth solver @type SolverAbstract
      %   scm: reference to the precomputed scm object @type SCM

      % save the references, nothing more
      obj.pd     = problem;
      obj.solver = truthsolver;
      obj.scm    = scm;

      % reset all variables to default values
      obj.resetTraining;
    end

    function [rbsolvec, errest] = onlineQuery(obj, param)
      % Use the reduced basis offline generated data to solve the problem for
      % the given parameter.
      %
      % Parameters:
      %   param: parameter to use. @type colvec
      %
      % Return values:
      %   rbsolvec: solution vector for the rb truth trial space @type colvec
      %   errest: reduced basis error estimate @type double

      % construct the reduced basis system
      testSnapshots = obj.assembleTestSnapshots(param);
      Lhs = testSnapshots.' * obj.solver.spacetimeSystemMatrix(param) * ...
              obj.trialSnapshots;
      Rhs = testSnapshots.' * obj.solver.spacetimeLoadVector;
      % now solve it
      rbsolvec = Lhs \ Rhs;
      % compute the error estimate
      errest = obj.estimateError(param, rbsolvec);
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

      % we are ready to start the greedy loop! first some more preparation
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
        if obj.nTrialRB > 0
          z = solvec - obj.trialSnapshots * (solvec.' * obj.solver.TrNorm * obj.trialSnapshots).';
          z = z / (sqrt(z' * obj.solver.TrNorm * z));
          obj.trialSnapshots(:, end + 1) = z;
        else
          z = z / (sqrt(z' * obj.solver.TrNorm * z));
          obj.trialSnapshots(:, 1) = z;
        end

        obj.nTrialRB = obj.nTrialRB + 1;

        % and compute the needed affine decomposition for the on-demand assembly
        % of the reduced basis test space

        % iterate over the space time matrix parts
        Zqn = zeros(obj.solver.nTestDim, obj.nQb);
        for pdx = 1:obj.nQb
          Zqn(:, pdx) = obj.solver.TeNorm \ (obj.solver.Lhs{pdx} * solvec);
        end
        obj.affineTestshotBase{end + 1} = Zqn;

        obj.nTestRB = obj.nTestRB + 1;

        % next step: compute all the error estimates, so let's prepare the
        % computation of the residual
        obj.prepareResidual;
        % and now lets iterate
        errests = zeros(size(paramTrain, 2), 1);
        for idx = 1:size(paramTrain, 2)
          [rbsolvec, errest] = obj.onlineQuery(paramTrain(:, idx));
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

    function snapshots = assembleTestSnapshots(obj, param)
      % Assemble the reduced basis test space for the given parameter.
      %
      % Parameters:
      %   param: parameter for the problem @type colvec
      %
      % Return values:
      %   snapshots: reduced basis test space @type matrix
      %
      % @todo describe it

      % prepad with 1 for the field independent part
      param = [1; shiftdim(param)];

      % evaluate the affine decomposition
      snapshots = sparse(obj.solver.nTestDim, obj.nTrialRB);
      for ndx = 1:obj.nTrialRB
        snapshots(:, ndx) = obj.affineTestshotBase{ndx} * param;
      end
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

    function resetTraining(obj)
      % Reset (or prepare) the object for the execution of the offline stage.
      %
      % Warning:
      %   This deletes all offline data!

      obj.trialSnapshots = sparse(obj.solver.nTrialDim, 1);
      obj.affineTestshotBase = cell(0);
      obj.errors = [];
      obj.params = [];
      obj.nTrialRB = 0;
      obj.nTestRB = 0;
    end
  end

  % collection of maybe necessary stuff

  methods(Access = 'protected')
    function prepareResidual(obj)
      % Set up the needed matrix for the computation of the residual.
      %
      % This assembles the solution and field independent part of the
      % computation of the Y'-norm of the residual through it's riesz
      % representation.
      %
      % @todo describe it more!

      % create the matrix
      H = zeros(obj.solver.nTestDim, obj.nQr);

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
  end

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

      testshots = obj.assembleTestSnapshots(param);

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

end % classdef
