classdef GalerkinSolver < handle
  % Galerkin-based Solver.
  %
  % @todo Not yet implemented.

  properties
    % time interval @type vector
    tspan = [0, 1];
    % spatial interval @type vector
    xspan = [0, 10];
    % field-switching point in time @type double
    fieldBreakpoint = 0.5;
    % multiplicative factor for the Laplacian @type double
    coeffLaplacian = 1;
    % additive field-offset `\mu`. @type double
    coeffOffset = 0;

    N = 10;
    M = 10;
    C = 10;
  end

  properties (Dependent)
    % checks if the field-offset is nonzero @type logical
    withOffset;
  end

  properties(Access = 'private')
    % Reference to assembly class @type AssemblyFourierLegendre
    assembly;

    partOne = {};
    partTwo = {};
  end

  methods

    % Constructor

    function obj = GalerkinSolver(N, M, C)
      % Parameters:
      %   N: number of spatial ansatz basis functions @type integer
      %   M: number of temporal ansatz basis functions @type integer
      %   N: number of spatial field basis functions @type integer

      obj.N = N;
      obj.M = M;
      obj.C = C;
    end


    % getters for dependent properties

    function val = get.withOffset(obj)
      % Prüfe, ob ein Offset bei den Feldern aktiviert ist.
      %
      % Return values:
      %   val: true genau dann wenn, die Felder um eine additive Konstante
      %     verschoben werden @type logical

      if obj.coeffOffset == 0
        val = false;
      else
        val = true;
      end
    end


    %% Preparation

    function preassemble(obj)
      obj.assembly                = AssemblyFourierLegendre();

      obj.assembly.setNumberOfAnsatzFuncs(obj.N, obj.M);
      obj.assembly.setNumberOfTestFuncsFromAnsatzFuncs();

      obj.assembly.xspan          = obj.xspan;
      obj.assembly.coeffLaplacian = obj.coeffLaplacian;
      obj.assembly.coeffOffset    = obj.coeffOffset;

      obj.assembly.tspan  = [obj.tspan(1), obj.fieldBreakpoint];
      [~, M1, M2, M3, M4F, M4B] = obj.assembly.assembleFieldIndependentMatrix();

      obj.partOne.TimeDerivativeMatrix = M1;
      obj.partOne.LaplacianMatrix      = M2;
      obj.partOne.OffsetMatrix         = M3;
      obj.partOne.ICForwardMatrix      = M4F;
      obj.partOne.ICBackwardMatrix     = M4B;
      obj.partOne.OmegaMatrices        = obj.assembly.assembleFieldDependentMatrixForFourierSeries(obj.C);


      obj.assembly.tspan  = [obj.fieldBreakpoint, obj.tspan(2)];
      [~, N1, N2, N3, N4F, N4B] = obj.assembly.assembleFieldIndependentMatrix();

      obj.partTwo.TimeDerivativeMatrix = N1;
      obj.partTwo.LaplacianMatrix      = N2;
      obj.partTwo.OffsetMatrix         = N3;
      obj.partTwo.ICForwardMatrix      = N4F;
      obj.partTwo.ICBackwardMatrix     = N4B;
      obj.partTwo.OmegaMatrices        = obj.assembly.assembleFieldDependentMatrixForFourierSeries(obj.C);
    end

    % Solver

    function solfun = solveTwoFieldsForward(obj, fieldCoeffsOne, fieldCoeffsTwo)
      % vorwärts
      Lhs = obj.partOne.TimeDerivativeMatrix ...
        + obj.coeffLaplacian * obj.partOne.LaplacianMatrix ...
        + obj.partOne.ICForwardMatrix;
      LhsComp = Lhs;
      for cdx = 1:obj.C
        LhsComp = LhsComp + fieldCoeffsOne(cdx) * obj.partOne.OmegaMatrices{cdx};
      end
      Rhs = obj.assembly.assembleVectorOnes();
      solCoeffOne = LhsComp \ Rhs;

      Lhs = obj.partTwo.TimeDerivativeMatrix ...
        + obj.coeffLaplacian * obj.partTwo.LaplacianMatrix ...
        + obj.partTwo.ICForwardMatrix;
      LhsComp = Lhs;
      for cdx = 1:obj.C
        LhsComp = LhsComp + fieldCoeffsTwo(cdx) * obj.partTwo.OmegaMatrices{cdx};
      end
      Rhs = obj.assembly.assembleVectorFromSolutionCoeffs(solCoeffOne);
      solCoeffTwo = LhsComp \ Rhs;

      solfun = @(t, x) ...
        (t < obj.fieldBreakpoint) .* obj.assembly.solutionFuncFromCoeffs(...
          solCoeffOne, t, x, [obj.tspan(1), obj.fieldBreakpoint]) ...
        + ~(t < obj.fieldBreakpoint) .* obj.assembly.solutionFuncFromCoeffs(...
          solCoeffTwo, t, x, [obj.fieldBreakpoint, obj.tspan(2)]);
    end

    function solfun = solveTwoFieldsBackward(obj, fieldCoeffsOne, fieldCoeffsTwo)
      % und das ganze nochmal rückwärts
      Lhs = - obj.partTwo.TimeDerivativeMatrix ...
        + obj.coeffLaplacian * obj.partTwo.LaplacianMatrix ...
        + obj.partTwo.ICBackwardMatrix;
      LhsComp = Lhs;
      for cdx = 1:obj.C
        LhsComp = LhsComp + fieldCoeffsTwo(cdx) * obj.partTwo.OmegaMatrices{cdx};
      end
      Rhs = obj.assembly.assembleVectorOnes();
      solCoeffTwo = LhsComp \ Rhs;

      Lhs = - obj.partOne.TimeDerivativeMatrix ...
        + obj.coeffLaplacian * obj.partOne.LaplacianMatrix ...
        + obj.partOne.ICBackwardMatrix;
      LhsComp = Lhs;
      for cdx = 1:obj.C
        LhsComp = LhsComp + fieldCoeffsOne(cdx) * obj.partOne.OmegaMatrices{cdx};
      end
      Rhs = obj.assembly.assembleVectorFromSolutionCoeffs(solCoeffTwo, true);
      solCoeffOne = LhsComp \ Rhs;

      solfun = @(t, x) ...
        (t < obj.fieldBreakpoint) .* obj.assembly.solutionFuncFromCoeffs(...
          solCoeffOne, t, x, [obj.tspan(1), obj.fieldBreakpoint]) ...
        + ~(t < obj.fieldBreakpoint) .* obj.assembly.solutionFuncFromCoeffs(...
          solCoeffTwo, t, x, [obj.fieldBreakpoint, obj.tspan(2)]);
    end

    function plotSolution(obj, solfun)
      gridt = linspace(obj.tspan(1), obj.tspan(2));
      gridx = linspace(obj.xspan(1), obj.xspan(2));
      % gridx = linspace(0, 25, 250);
      [mesht, meshx] = meshgrid(gridt, gridx);

      figure();
      mesh(mesht, meshx, solfun(mesht, meshx));
    end

    function solc = solveTwoFieldsDeprecated(obj, lap, cof1, cof2, xg, tgl, tgr)
      % @deprecated
      N = 25;
      M = 25;
      C = 50;

      assembly = AssemblyFourierLegendre();
      assembly.setNumberOfAnsatzFuncs(N, M);
      assembly.setNumberOfTestFuncsFromAnsatzFuncs();
      assembly.xspan = obj.xspan;

      assembly.coeffLaplacian = lap;
      assembly.coeffOffset = 0;

      % solve for first field
      assembly.tspan = [obj.tspan(1), obj.fieldBreakpoint];
      % assembly.initialData = @(x) ones(size(x, 1), size(x, 2));

      LHS1 = assembly.assembleFieldIndependentMatrix();
      O1 = assembly.assembleFieldDependentMatrixForFourierSeries(C);
      RHS1 = assembly.assembleVectorOnes();

      % compose field
      LHSO1 = LHS1;
      for idx = 1:length(cof1)
        LHSO1 = LHSO1 + cof1(idx) * O1{idx};
      end

      size(LHS1)

      sol1 = LHSO1 \ RHS1;

      % gridt1 = linspace(obj.tspan(1), obj.fieldBreakpoint);
      % gridx1 = linspace(obj.xspan(1), obj.xspan(2));
      [mesht1, meshx1] = meshgrid(tgl(1:end-1), xg);
      tic
      soleval1  = assembly.solutionFuncFromCoeffs(sol1, mesht1, meshx1);
      toc

      % solve for second field
      assembly.tspan = [obj.fieldBreakpoint, obj.tspan(2)];

      LHS2 = assembly.assembleFieldIndependentMatrix();
      O2 = assembly.assembleFieldDependentMatrixForFourierSeries(C);

      cf = zeros(N, 1);
      for idx = 1:N
        cf(idx) = sum(sol1(((idx - 1) * M + 1): (idx * M)));
      end

      RHS2 = assembly.assembleVectorFromSpatialCoeffs(cf);

      % compose field
      LHSO2 = LHS2;
      for idx = 1:length(cof2)
        LHSO2 = LHSO2 + cof2(idx) * O2{idx};
      end

      sol2 = LHSO2 \ RHS2;

      %% visualization stuff
      % gridt2 = linspace(obj.fieldBreakpoint, obj.tspan(2));
      % gridx2 = linspace(obj.xspan(1), obj.xspan(2));
      [mesht2, meshx2] = meshgrid(tgr, xg);

      % evaluate the solution function
      tic
      soleval2  = assembly.solutionFuncFromCoeffs(sol2, mesht2, meshx2);
      toc

      % composite plot
      meshtc = [mesht1, mesht2];
      meshxc = [meshx1, meshx2];
      solc = [soleval1, soleval2];

      % figure()
      % mesh(meshtc, meshxc, solc)
    end

    % Plotting and visualization

    function plotSolutionDeprecated(obj, solfun)
      % @deprecated
      % Plot the solution.
      %
      % Parameters:
      %   solfun: function handle of the solution function @type function_handle

      % create the grids
      gridt = linspace(obj.tspan(1), obj.tspan(2));
      gridx = linspace(obj.xspan(1), obj.xspan(2));
      [mesht, meshx] = meshgrid(gridt, gridx);

      % evaluate the solution function
      solution = solfun(mesht, meshx);

      % visualize
      figure();
      mesh(mesht, meshx, solution);
      title('Solution');
      xlabel('t');
      ylabel('x');
      zlabel('u');
    end

  end

  methods(Access = Private)

  end

end
