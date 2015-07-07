classdef RBMSolverFourierNodal < RBMSolverAbstract
  % Solver based on the reduced basis method.

  properties
    % nC;

    % @type cellarray
    affineTestshotBase;
  end

  properties %(Access = 'protected')
  end

  properties(Dependent)
    nQf;
    nQb;
    nQr;
    nNrb;
    nNfeX;
    nNfeY;
  end

  methods
    function val = get.nQf(obj)
      val = 1;
    end

    function val = get.nQb(obj)
      val = 1 + obj.nFields * obj.nC;
    end

    function val = get.nQr(obj)
      val = obj.nQf + obj.nNrb * obj.nQb;
    end

    function val = get.nNrb(obj)
      val = size(obj.trialshots, 2);
    end

    function val = get.nNfeX(obj)
      val = obj.solver.nTrialDim;
    end

    function val = get.nNfeY(obj)
      val = obj.solver.nTestDim;
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


      obj.solver.nTrialS    = 20;
      obj.solver.nTrialT    = 100;
      obj.solver.nTestS     = obj.solver.nTrialS;
      obj.solver.nTestT     = obj.solver.nTrialT - 1;
      obj.solver.nTestSic   = obj.solver.nTrialS;
      obj.solver.tgrid      = linspace(0, 1, obj.solver.nTrialT);
      obj.solver.tspan      = [0 1];
      obj.solver.xspan      = [0 10];
      obj.solver.cLaplacian = 3.3333;
      obj.solver.cOffset    = 0;
      obj.solver.breakpoints = obj.breakpoints;
      obj.solver.nFieldCoeffs = obj.nC;

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

      testshots = obj.constructRBTestSpace(param);
      B = testshots.' * obj.solver.spacetimeSystemMatrix(param) * obj.trialshots;
      F = obj.testshots.' * obj.solver.spacetimeLoadVector;
    end

    function rbsolvec = solveRB(obj, param, RhsPre)
      [Lhs, Rhs] = obj.rbSystemMatrixAndLoadVector(param);
      rbsolvec   = Lhs \ Rhs;
    end

    function testshots = constructRBTestSpace(obj, param)
      testshots = sparse(obj.solver.nTestDim, size(obj.trialshots, 2));
      for ndx = 1:size(obj.trialshots, 2)
        tmp = obj.affineTestshotBase{ndx}(:, 1);
        for fdx = 1:obj.nFields
          for cdx = 1:obj.nC
            kdx = 1 + (fdx - 1) * obj.nC + cdx;
            tmp = tmp + param(cdx, fdx) * obj.affineTestshotBase{ndx}(:, kdx);
          end
        end
        testshots(:, ndx) = tmp;
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

      % gram-schmidt-verfahren laufen lassen
      z = gsolvec;
      if ~isFirst
        for j = 1:size(obj.trialshots, 2)
            z = z - (gsolvec.' * obj.solver.TrNorm * obj.trialshots(:, j)) * obj.trialshots(:, j);
        end
      end
      z = z / (sqrt(z' * obj.solver.TrNorm * z));
      obj.trialshots(:, pos) = z;

      % zqn update
      Zqn = sparse(obj.solver.nTestDim, obj.nQb);
      Zqn = obj.solver.TeNorm \ (obj.solver.LhsFI * gsolvec);
      for fdx = 1:obj.nFields
        for cdx = 1:obj.nC
          kdx = 1 + (fdx - 1) * obj.nC + cdx;
          Zqn(:, kdx) = obj.solver.TeNorm \ (obj.solver.LhsFD{cdx, fdx} * gsolvec);
        end
      end
      obj.affineTestshotBase{pos} = Zqn;
    end


    function offlinePhase(obj)
      % solve the system for some parameters
      % [gsolvec, LhsPre, RhsPre] = obj.solver.solve({[0 0 0]});

      obj.solveAndAdd([0; 0; 0], true)
      obj.solveAndAdd([1; 0; 0]);
      obj.solveAndAdd([0; 2; 0]);
      obj.solveAndAdd([3; 0; 0]);
      obj.solveAndAdd([0; 4; 0]);
      obj.solveAndAdd([5; 0; 0]);
      obj.solveAndAdd([1; 0; 3]);
      obj.solveAndAdd([7; 0; 0]);
      % obj.solveAndAdd({[8 0 0]});
      % obj.solveAndAdd({[9 0 0]});
      % obj.solveAndAdd({[10 0 0]});
      % obj.solveAndAdd({[0 0 0]});
      % obj.solveAndAdd({[0 0 10]});
      % obj.solveAndAdd({[0 0 0]});

      testparam = [1; 1; 0];
      obj.testshots = obj.constructRBTestSpace(testparam);
      obj.calcDiscreteInfSupAndContinuityRB(testparam)



      xg = linspace(obj.solver.xspan(1), obj.solver.xspan(2), 50);
      rbsolvec = obj.solveRB(testparam, obj.solver.spacetimeLoadVector);
      [gsolvec2, ~, ~] = obj.solver.solve(testparam);

      figure(1);
      mesh(obj.evaluateSolution(gsolvec2, xg));
      figure(2)
      mesh(obj.evaluateRBMSolution(rbsolvec, xg));
      figure(3)
      mesh(abs(obj.evaluateRBMSolution(rbsolvec, xg) - obj.evaluateSolution(gsolvec2, xg)));

      % check

      return;

      % [gsolvec, LhsPre, RhsPre] = obj.solver.solve({[1 0 0]});
      % gtestvec = obj.solver.TeNorm \ (LhsPre * gsolvec);
      % obj.trialshots(:, end + 1) = gsolvec;
      % obj.testshots(:, end + 1)  = gtestvec;

      [gsolvec, LhsPre, RhsPre] = obj.solver.solve({[0 1 0]});
      gtestvec = obj.solver.TeNorm \ (LhsPre * gsolvec);
      obj.trialshots(:, end + 1) = gsolvec;
      obj.testshots(:, end + 1)  = gtestvec;

      [gsolvec, LhsPre, RhsPre] = obj.solver.solve({[0 0 1]});
      gtestvec = obj.solver.TeNorm \ (LhsPre * gsolvec);
      obj.trialshots(:, end + 1) = gsolvec;
      obj.testshots(:, end + 1)  = gtestvec;

      return;
      % @todo fix that stuff!

      % save the solution as a snapshot
      % @todo normalize that vector!
      obj.trialshots(:, 1) = gsolvec;
      % set up the stuff for the residual
      H = zeros(obj.nNfeY, obj.nQr);
      H(:, 1) = obj.solver.rhs();

      ynorm = obj.solver.TeNorm;

      H(:, 2:obj.nQb:end) = - obj.solver.LhsFI.F * obj.trialshots;

      for fdx = 1:obj.nFields
        for cdx = 1:obj.nC
          kdx = 1 + (fdx - 1) * obj.nC + cdx;
          H(:, (1+kdx):obj.nQb:end) = - obj.solver.LhsFD{cdx, fdx} * obj.trialshots;
        end
      end

      G = H.' * (ynorm \ H);

      %%

      u_rb = 1;

      Eps    = zeros(obj.nQr, 1);
      Eps(1) = 1;
      Eps(:, 2:obj.nQb:end) = 1 * u_rb;

      for fdx = 1:obj.nFields
        for cdx = 1:obj.nC
          kdx = 1 + (fdx - 1) * obj.nC + cdx;
          Eps((1+kdx):obj.nQb:end) = param2(cdx, fdx) * u_rb;
        end
      end

      size(Eps)

      e_hat_X  = sqrt(abs(Eps.' * G * Eps))
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

    end

    function val = residual(obj, gsolvec, param)
      % Compute the Y-Norm of the residual.
      %
      % Parameters:
      %   param: Parameter for which the residual should be computed.
      %
      % Return values:
      %   val: Y-Norm of the Residual @type double

      rhs = obj.solver.rhs();
      ynorm = obj.solver.TeNorm;

      % first part of the residual
      fst = rhs.' * (ynorm \ rhs);

      % second part of the residual
      % field independent part first
      snd = rhs.' * (ynorm \ (obj.solver.LhsFI.F * gsolvec));
      for fdx = 1:obj.nFields
        for cdx = 1:obj.nC
          snd = snd + param(cdx, fdx) * rhs.' * (ynorm \ (obj.solver.LhsFD{cdx, fdx} * gsolvec));
        end
      end

      % third part of the residual
      trd = (obj.solver.LhsFI.F * gsolvec).' * (ynorm \ (obj.solver.LhsFI.F * gsolvec));
      for fdx = 1:obj.nFields
        for cdx = 1:obj.nC
          trd = trd + 2 * param(cdx, fdx) * (obj.solver.LhsFI.F * gsolvec).' * (ynorm \ (obj.solver.LhsFD{cdx, fdx} * gsolvec));
        end
      end
      for fdx1 = 1:obj.nFields
        for fdx2 = 1:obj.nFields
          for cdx1 = 1:obj.nC
            for cdx2 = 1:obj.nC
              trd = trd + param(cdx1, fdx1) * param(cdx2, fdx2) * (obj.solver.LhsFD{cdx1, fdx1} * gsolvec).' * (ynorm \ (obj.solver.LhsFD{cdx2, fdx2} * gsolvec));
            end
          end
        end
      end

      % all together now
      val = sqrt(fst - 2 * snd + trd);
    end

    function solval = evaluateRBMSolution(obj, rbsolvec, xgrid)
      % @todo Not yet fully implemented!

      solval = obj.evaluateSolution(obj.trialshots * rbsolvec, xgrid);
    end

    function solval = evaluateSolution(obj, gsolvec, xgrid)
      % Evaluate a solution of the truth galerkin solver.
      %
      % Parameters:
      %   gsolvec: solution vector of the truth galerkin trial space.
      %     @type vector
      %   xgrid: grid of x values in which the solution should be evaluated.
      %     @type vector
      %
      % Return values:
      %   solval: evaluated solution @type matrix

      % Just forward it to the solver.
      solval = obj.solver.evaluateSolution(gsolvec, xgrid);
    end

  end

  % collection of maybe necessary stuff

  methods
    function [evs, evl] = calcDiscreteInfSupAndContinuityRB(obj, param)
      % Calculate the discrete inf-sup-constant of the rb system for the given
      % parameter.
      %
      % Parameters:
      %   param: field parameters @type matrix
      %
      % Return values:
      %   beta: discrete inf-sup-constant @type double
      %   gamma: discrete continuity constant @type double

      % first we have to construct the rb system matrix for this parameter and
      % the y- and x-norm matrices.
      Lhs   = full(obj.testshots.' * obj.solver.spacetimeSystemMatrix(param) * obj.trialshots);
      Ynorm = full(obj.testshots.' * obj.solver.TeNorm * obj.testshots);
      Xnorm = full(obj.trialshots.' * obj.solver.TrNorm * obj.trialshots);

      % now we solve the generalized eigenvalue problem, the inf-sup and continuity
      % constants are then the square roots of the smallest respectively largest
      % generalized eigenvalue
      supNorm = Lhs.' * (Ynorm \ Lhs);
      ev = eig(supNorm, Xnorm);
      beta = sqrt(min(ev));
      gamma = sqrt(max(ev));
    end

    function [beta, gamma] = calcDiscreteInfSupAndContinuity(obj, param)
      % Calculate the discrete inf-sup-constant of the truth system for the
      % given parameter.
      %
      % Parameters:
      %   param: field parameters @type matrix
      %
      % Return values:
      %   beta: discrete inf-sup-constant @type double
      %   gamma: discrete continuity constant @type double

      % get all the needed matrices from the truth solver
      Lhs = obj.solver.spacetimeSystemMatrix(param);
      Ynorm = obj.solver.TeNorm;
      Xnorm = obj.solver.TrNorm;

      % solve the generalized eigenvalue problem, the inf-sup and continuity
      % constants are the square root of the smallest respectively largest
      % generalized eigenvalue.
      supNorm = Lhs.' * (Ynorm \ Lhs);
      beta    = sqrt(eigs(supNorm, Xnorm, 1, 'sm'));
      gamma   = sqrt(eigs(supNorm, Xnorm, 1, 'lm'));
    end
  end

end
