classdef RBMSolverFourierNodal < RBMSolverAbstract
  % Solver based on the reduced basis method.

  properties
    % nC;

    % @type cellarray
    affineTestshotBase;
  end

  properties %(Access = 'protected')
    G;
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
      val = 1 + obj.nFields * obj.nC;
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
      val = size(obj.trialshots, 2);
    end

    function val = get.nTestRB(obj)
      val = obj.nTrialRB;
    end
  end

  methods

    function obj = RBMSolverFourierNodal()
      % Constructor for this class.

      % first call the constructor of the superclass
      obj@RBMSolverAbstract();

      % and now the rest
      obj.solver = SolverFourierNodal();
    end

    function prepare(obj)
      % Prepare the reduced basis method solver.
      %
      % This mainly consists of forwarding needed values to underlying classes
      % and setting up the galerkin solver.
      %
      % @todo refactor this pile of junk...


      obj.solver.nTrialS    = 50;
      obj.solver.nTrialT    = 100;
      obj.solver.nTestS     = obj.solver.nTrialS;
      obj.solver.nTestT     = obj.solver.nTrialT - 1;
      obj.solver.nTestSic   = obj.solver.nTrialS;
      obj.solver.tgrid      = linspace(0, 1, obj.solver.nTrialT);
      obj.solver.tspan      = [0 1];
      obj.solver.xspan      = [0 1];
      obj.solver.cLaplacian = 3.3333;
      obj.solver.cOffset    = 0;
      obj.solver.useRefinement = true;
      obj.solver.breakpoints = obj.breakpoints;
      % obj.nC = 1;
      obj.solver.nC = obj.nC;

      obj.solver.prepare;

      obj.trialshots = sparse(obj.solver.nTrialDim, 1);
      obj.affineTestshotBase = cell(1, 1);

    end

    function [B, F] = rbSystemMatrixAndLoadVector(obj, param)
      % Assemble the rb system matrix and load vector.
      %
      % This relies on the affine representation of the space time system
      % matrix.
      %
      % Parameters:
      %   param: field coefficients. @type struct
      %
      % Return values:
      %   B: rb system matrix for the given parameter @type matrix
      %   F: rb load vector @type vector

      testshots = obj.constructTestSpace(param);
      B = testshots.' * obj.solver.spacetimeSystemMatrix(param) * obj.trialshots;
      F = testshots.' * obj.solver.spacetimeLoadVector;
    end

    function rbsolvec = rbSolve(obj, param)
      [Lhs, Rhs] = obj.rbSystemMatrixAndLoadVector(param);
      rbsolvec   = Lhs \ Rhs;
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

    function solveAndAdd(obj, param, isFirst)
      if nargin == 2
        isFirst = false;
      end

      if isFirst
        pos = 1;
      else
        pos = size(obj.affineTestshotBase, 2) + 1;
      end

      [gsolvec, ~, ~] = obj.solver.solve(param);

      %| @todo aufräumen
      % gram-schmidt-verfahren laufen lassen
      z = gsolvec;
      if ~isFirst
        z = gsolvec - obj.trialshots * (gsolvec.' * obj.solver.TrNorm * obj.trialshots).';
      end
      z = z / (sqrt(z' * obj.solver.TrNorm * z));
      obj.trialshots(:, pos) = z;

      % zqn update
      % iterate over the space time matrix parts
      Zqn = zeros(obj.nTestTruth, obj.nQb);
      for pdx = 1:obj.nQb
        Zqn(:, pdx) = obj.solver.TeNorm \ (obj.solver.Lhs{pdx} * gsolvec);
      end
      obj.affineTestshotBase{pos} = Zqn;
    end

    function gramSchmidtOrtho(obj)
      % Execute the gram schmidt orthonormalization.
      %
      % @deprecated Not working correctly for more that a handful trial
      %     snapshots

      % compute a cholesky decomposition of the X-norm on the ansatz space. the
      % coefficients for the gram schmidt orthonormalization are found in the
      % inverse-transpose of the lower left triangle.
      K = obj.trialshots.' * obj.solver.TrNorm * obj.trialshots;
      L = chol(K, 'lower');
      C = inv(L.');
      % C = L.' \ eye(size(L));

      % full(C)

      % orthonormalize the trial rb basis
      trialshots = obj.trialshots * C;

      obj.trialshots = trialshots;
    end


    function offlinePhase(obj)
      % solve the system for some parameters
      % [gsolvec, LhsPre, RhsPre] = obj.solver.solve({[0 0 0]});

      % obj.solveAndAdd([0], true)
      % obj.solveAndAdd([2])
      % obj.solveAndAdd([4])
      [b, g] = obj.calcDiscreteInfSupAndContinuityTruth([-3]);
      b = b^2

      scm = SCM(obj);

      % [b, g] = obj.calcDiscreteInfSupAndContinuityTruth([3/2]);
      % b = b^2
      % g = g^2
      scm.offlinePhase


      return;
      obj.solveAndAdd([0], true)
      obj.solveAndAdd([0], true)
      % obj.gramSchmidtOrtho
      for idx = 1:2:10
        obj.solveAndAdd([idx]);
      end
      obj.solveAndAdd([8]);
      obj.solveAndAdd([8.5]);
      obj.solveAndAdd([9.5]);
      obj.solveAndAdd([10]);

      % return

      obj.prepareResidual;

      realerr = [];
      bounderr = [];

      for cof = 0:0.1:10
        rbsolvec = obj.rbSolve([cof]);
        res = obj.residual(rbsolvec, [cof]);
        [gsolvec2, ~, ~] = obj.solver.solve([cof]);

        [b, g] = obj.calcDiscreteInfSupAndContinuityTruth([cof])
        bounderr(end + 1) = res / b;
        realerr(end + 1) = sqrt((obj.trialshots * rbsolvec - gsolvec2).' * obj.solver.TrNorm * (obj.trialshots * rbsolvec - gsolvec2));
        [b, g] = obj.calcDiscreteInfSupAndContinuityRB([cof])
      end

      figure(10)
      semilogy(0:0.1:10, realerr, 0:0.1:10, bounderr)
      legend('real', 'bound')

      return;

      testparam = [2; 0; 0];
      obj.testshots = obj.constructTestSpace(testparam);
      [bta, gma] = obj.calcDiscreteInfSupAndContinuityRB(testparam)

      testparam = [2; 0; 0];
      rbsolvec = obj.rbSolve(testparam);
      res = obj.residual(rbsolvec, testparam)

      testparam = [2; 0; 0];
      rbsolvec = obj.rbSolve(testparam);
      res = obj.residual(rbsolvec, testparam)


      % testparam = [2; 0; 0];
      % rbsolvec = obj.rbSolve(testparam);
      % res = obj.residual(rbsolvec, testparam)

      xg = linspace(obj.solver.xspan(1), obj.solver.xspan(2), 50);
      [gsolvec2, ~, ~] = obj.solver.solve(testparam);

      % estimate error
      errbound = res / obj.calcDiscreteInfSupAndContinuityTruth(testparam)
      errreal = sqrt((obj.trialshots * rbsolvec - gsolvec2).' * obj.solver.TrNorm * (obj.trialshots * rbsolvec - gsolvec2))

      figure(1);
      mesh(obj.evaluateSolutionTruth(gsolvec2, xg));
      figure(2)
      mesh(obj.evaluateSolutionRb(rbsolvec, xg));
      figure(3)
      mesh(abs(obj.evaluateSolutionRb(rbsolvec, xg) - obj.evaluateSolutionTruth(gsolvec2, xg)));

    end

    function gsolvec = offlineSolve(obj, param)
      % Solve the pde for the given parameters with the truth galerkin solver.
      %
      % Parameters:
      %   param: cellarray containing the coefficients of the series expansions
      %     of the fields. @type cellarray
      %
      % Return values:
      %   gsolvec: truth solution vector @type vector
      %
      % @deprecated replaced with solveAndAdd

      [gsolvec, LhsPre, RhsPre] = obj.solver.solve(param);

      % calculate test vector through supremizer
      gtestvec = obj.solver.TeNorm \ (LhsPre * gsolvec);

      % save both vectors
      trialshots(:, end + 1) = gsolvec;
      testshots(:, end + 1)  = gtestvec;

      % Arb = gtestvec.' * LhsPre * gsolvec;
      % Frb = gtestvec.' * RhsPre;

      % Arb \ Frb
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
      H = zeros(obj.nTrialTruth, obj.nQr);

      % the first column is the rhs of the truth system
      H(:, 1) = obj.solver.spacetimeLoadVector;;

      % the remaining columns are the products of the affine composition of the
      % truth system matrix and the chosen trial truth space snapshots
      for pdx = 1:obj.nQb
        H(:, (pdx + 1):obj.nQb:end) = - obj.solver.Lhs{pdx} * obj.trialshots;
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

    function solval = evaluateSolutionRb(obj, rbsolvec, xgrid)
      % Evaluate a solution of the reduced basis system.
      %
      % Parameters:
      %   rbsolvec: solution vector of the rb system @type vector
      %   xgrid: spatial grid @type vector
      %
      % Return values:
      %   solval: values of the solution on the grid points @type matrix

      % forward the resulting vector of the truth trial space to the
      % corresponding evaluation method
      solval = obj.evaluateSolutionTruth(obj.trialshots * rbsolvec, xgrid);
    end

    function solval = evaluateSolutionTruth(obj, gsolvec, xgrid)
      % Evaluate a solution of the truth solver.
      %
      % Parameters:
      %   gsolvec: solution vector of the truth trial space @type vector
      %   xgrid: spatial grid @type vector
      %
      % Return values:
      %   solval: values of the solution on the grid points @type matrix

      % Just forward it to the truth solver method
      solval = obj.solver.evaluateSolution(gsolvec, xgrid);
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
      Lhs   = full(testshots.' * obj.solver.spacetimeSystemMatrix(param) * obj.trialshots);
      Ynorm = full(testshots.' * obj.solver.TeNorm * testshots);
      Xnorm = full(obj.trialshots.' * obj.solver.TrNorm * obj.trialshots);

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
