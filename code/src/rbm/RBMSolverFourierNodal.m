classdef RBMSolverFourierNodal < RBMSolverAbstract
  % Solver based on the reduced basis method.

  properties
    % nC;

    % @type cellarray
    affineTestshotBase;
  end

  properties %(Access = 'protected')
    G;
    tmp = 0;
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

    function rbsolvec = solveRB(obj, param)
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

      %| @todo aufräumen
      % gram-schmidt-verfahren laufen lassen
      z = gsolvec;
      if ~isFirst
        % for j = 1:size(obj.trialshots, 2)
        %     z = z - (gsolvec.' * obj.solver.TrNorm * obj.trialshots(:, j)) * obj.trialshots(:, j);
        % end
        z = gsolvec - obj.trialshots * (gsolvec.' * obj.solver.TrNorm * obj.trialshots).';
      end
      z = z / (sqrt(z' * obj.solver.TrNorm * z));
      obj.trialshots(:, pos) = z;
      % obj.tmp = obj.tmp + toc;

      % zqn update
      Zqn = obj.solver.TeNorm \ (obj.solver.LhsFI * gsolvec);
      for fdx = 1:obj.nFields
        for cdx = 1:obj.nC
          kdx = 1 + (fdx - 1) * obj.nC + cdx;
          Zqn(:, kdx) = obj.solver.TeNorm \ (obj.solver.LhsFD{cdx, fdx} * gsolvec);
        end
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

      obj.solveAndAdd([0; 0; 0], true)
      % obj.gramSchmidtOrtho
      for idx = 1:10:100
        obj.solveAndAdd([idx; 0; 0]);
      end
      obj.tmp

      return

      obj.prepareResidual;

      testparam = [0; 0; 0];
      obj.testshots = obj.constructRBTestSpace(testparam);
      [bta, gma] = obj.calcDiscreteInfSupAndContinuityRB(testparam)

      testparam = [0; 0; 0];
      rbsolvec = obj.solveRB(testparam);
      res = obj.residual(rbsolvec, testparam)

      testparam = [1; 0; 0];
      rbsolvec = obj.solveRB(testparam);
      res = obj.residual(rbsolvec, testparam)

      testparam = [2; 0; 0];
      rbsolvec = obj.solveRB(testparam);
      res = obj.residual(rbsolvec, testparam)

      xg = linspace(obj.solver.xspan(1), obj.solver.xspan(2), 50);
      [gsolvec2, ~, ~] = obj.solver.solve(testparam);

      figure(1);
      mesh(obj.evaluateSolution(gsolvec2, xg));
      figure(2)
      mesh(obj.evaluateRBMSolution(rbsolvec, xg));
      figure(3)
      mesh(abs(obj.evaluateRBMSolution(rbsolvec, xg) - obj.evaluateSolution(gsolvec2, xg)));

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
      H(:, 2:obj.nQb:end) = - obj.solver.LhsFI * obj.trialshots;
      for fdx = 1:obj.nFields
        for cdx = 1:obj.nC
          kdx = 1 + (fdx - 1) * obj.nC + cdx;
          H(:, (1+kdx):obj.nQb:end) = - obj.solver.LhsFD{cdx, fdx} * obj.trialshots;
        end
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

      % and the remaining parts are the field dependent parts, also spaced over
      % the whole vector
      for fdx = 1:obj.nFields
        for cdx = 1:obj.nC
          kdx = 1 + (fdx - 1) * obj.nC + cdx;
          Eps((1+kdx):obj.nQb:end) = param(cdx, fdx) * uN;
        end
      end

      % and now calculate the y-norm of the riesz representation
      val  = sqrt(abs(Eps.' * obj.G * Eps));
    end

    function val = residual2(obj, uN, param)
      % Compute the Y-Norm of the residual.
      %
      % @deprecated somethings really wrong here...
      %
      % Parameters:
      %   param: Parameter for which the residual should be computed.
      %
      % Return values:
      %   val: Y-Norm of the Residual @type double


      % got everything right here
      Y = obj.solver.TeNorm;
      iY = inv(Y);
      F = obj.solver.spacetimeLoadVector;
      B0 = obj.solver.LhsFI;
      Bq = obj.solver.LhsFD;

      % vorberechnen...
      BqF = cell(obj.nQb, obj.nTrialRB);
      for m = 1:obj.nTrialRB
        BqF{1, m} = F.' * iY * B0 * obj.trialshots(:, m);
        for q = 2:obj.nQb
          BqF{q, m} = F.' * iY * Bq{q - 1, 1} * obj.trialshots(:, m);
        end
      end

      BqBq = cell(obj.nQb, obj.nQb, obj.nTrialRB, obj.nTrialRB);
      for m1 = 1:obj.nTrialRB
        for m2 = 1:obj.nTrialRB

          BqBq{1, 1, m1, m2} = (B0 * obj.trialshots(:, m1)).' * iY * (B0 * obj.trialshots(:, m2));

          for q = 2:obj.nQb
            BqBq{q, 1, m1, m2} = (B0 * obj.trialshots(:, m1)).' * iY * (Bq{q - 1, 1} * obj.trialshots(:, m2));
            BqBq{1, q, m1, m2} = (Bq{q - 1, 1} * obj.trialshots(:, m1)).' * iY * (B0 * obj.trialshots(:, m2));
          end

          for q1 = 2:obj.nQb
            for q2 = 2:obj.nQb
              BqBq{q1, q2, m1, m2} = (Bq{q1 - 1, 1} * obj.trialshots(:, m1)).' * iY * (Bq{q2 - 1, 1} * obj.trialshots(:, m2));
            end
          end
        end
      end


      % first part of the residual
      fst = F.' * iY * F;
      % fst = full(F.' * (Y \ F));

      % second part of the residual
      % field independent part first
      sndvec = zeros(1, 1);
      for ndx = 1:obj.nTrialRB
        sndvec(end+1) = uN(ndx) * BqF{1, ndx};
        % sndvec(end+1) = uN(ndx) * F.' * iY * B0 * obj.trialshots(:, ndx);
        % sndvec(end+1) = uN(ndx) * F.' * (Y \ (B0 * obj.trialshots(:, ndx)));
        for fdx = 1:obj.nFields
          for cdx = 1:obj.nC
            pos = 1 + (fdx - 1) * obj.nC + cdx;
            sndvec(end+1) = param(cdx, fdx) * uN(ndx) * BqF{pos, ndx};
          end
        end
      end

      snd = sum(sndvec)

      % third part of the residual
      trdvec = zeros(1, 1);
      for ndx1 = 1:obj.nTrialRB
        for ndx2 = 1:obj.nTrialRB

          trdvec(end+1) = uN(ndx1) * uN(ndx2) * BqBq{1, 1, ndx1, ndx2};
          % trdvec(end+1) = uN(ndx1) * uN(ndx2) * (B0 * obj.trialshots(:, ndx1)).' * iY * B0 * obj.trialshots(:, ndx2);

          for fdx = 1:obj.nFields
            for cdx = 1:obj.nC
              pos = 1 + (fdx - 1) * obj.nC + cdx;
              trdvec(end+1) = param(cdx, fdx) * uN(ndx1) * uN(ndx2) * BqBq{pos, 1, ndx1, ndx2};
              trdvec(end+1) = param(cdx, fdx) * uN(ndx1) * uN(ndx2) * BqBq{1, pos, ndx1, ndx2};
              % trdvec(end+1) = param(cdx, fdx) * uN(ndx1) * uN(ndx2) * (B0 * obj.trialshots(:, ndx1)).' * iY * Bq{cdx, fdx} * obj.trialshots(:, ndx2);
            end
          end

          % for fdx = 1:obj.nFields
          %   for cdx = 1:obj.nC
          %     trdvec(end+1) = param(cdx, fdx) * uN(ndx1) * uN(ndx2) * (Bq{cdx, fdx} * obj.trialshots(:, ndx1)).' * iY * B0 * obj.trialshots(:, ndx2);
          %   end
          % end

          for fdx1 = obj.nFields
            for fdx2 = obj.nFields
              for cdx1 = obj.nC
                for cdx2 = obj.nC
                  pos1 = 1 + (fdx1 - 1) * obj.nC + cdx1;
                  pos2 = 1 + (fdx2 - 1) * obj.nC + cdx2;
                  trdvec(end+1) = param(cdx1, fdx1) * param(cdx2, fdx2) * uN(ndx1) * uN(ndx2) * BqBq{pos1, pos2, ndx1, ndx2};
                end
              end
            end
          end
        end
      end

      trd = sum(trdvec)

      % all together now
      fst
      snd
      trd
      val = fst - 2 * snd + trd
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

  methods(Access = 'protected')
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
