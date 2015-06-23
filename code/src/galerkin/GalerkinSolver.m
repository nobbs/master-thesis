classdef GalerkinSolver < handle
  % Galerkin-based Solver.
  %
  % @todo Not yet implemented.

  properties
    % time interval @type vector
    tspan = [0, 1];
    % spatial interval @type vector
    xspan = [0, 10];
    % initial condition
    initialData = 0;
    % source term
    sourceData = 0;
    % first external field
    fieldA = 0;
    % second external field
    fieldB = 0;
    % field-switching point in time @type double
    fieldBreakpoint = 0.5;
    % multiplicative factor for the Laplacian @type double
    coeffLaplacian = 1;
    % additive field-offset `\mu`. @type double
    coeffOffset = 0;
  end

  properties (Dependent)
    % checks if the field-offset is nonzero @type logical
    withOffset;
  end

  methods

    % Constructor

%     function obj = GalerkinSolver()
%       obj;
%     end


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


    % setup of the initial conditions

    function setInitialDataFromFunction(obj, ufun)
      % Set the initial conditions through a function handle.
      %
      % @todo Not yet implemented!
      %
      % Parameters:
      %   ufun: @type function_handle

      error('Not yet implemented!');
    end


    function setInitialDataPointwise(obj)
      % Set the initial conditions through a vector of function values on an
      % equidistant grid.
      %
      % @todo Not yet implemented!
      error('Not yet implemented!');
    end

    function setInitialDataFromFourierCoeffs(obj)
      % Set the initial conditions through a coefficients vector of a Fourier
      % series.
      %
      % @todo Not yet implemented!
      error('Not yet implemented!');
    end

    function setInitialDataFromSineCoeffs(obj)
      % Set the initial conditions through a coefficients vector of a Sine
      % series.
      %
      % @todo Not yet implemented!
      error('Not yet implemented!');
    end


    % setup of the source term

    function setSourceDataToZero(obj)
      % Resets the source data to the default value of zero.
    end

    function setSourceDataPointwise(obj)
      % @todo Not yet implemented!
      error('Not yet implemented!');
    end

    function setSourceDataFromFunc(obj)
      % @todo Not yet implemented!
      error('Not yet implemented!');
    end

    % setup of the fields

    function setFieldFromFunc(obj)
      % @todo Not yet implemented!
      error('Not yet implemented!');
    end

    function setFieldPointwise(obj)
      % @todo Not yet implemented!
      error('Not yet implemented!');
    end

    function setFieldFromFourierCoeffs(obj)
      % @todo Not yet implemented!
      error('Not yet implemented!');
    end

    function setFieldFromSineCoeffs(obj)
      % @todo Not yet implemented!
      error('Not yet implemented!');
    end


    % Solver

    function [solfun] = solve(obj, lap, cof)
      % @todo Not yet implemented!

      % assembly = AssemblySineLegendre();
      assembly = AssemblyFourierLegendre();
      assembly.setNumberOfAnsatzFuncs(40, 40);
      assembly.setNumberOfTestFuncsFromAnsatzFuncs();

      obj.tspan = [0 1];
      assembly.tspan = obj.tspan;
      assembly.xspan = obj.xspan;

      assembly.coeffLaplacian = lap;
      assembly.coeffOffset = 5;

      % assembly.initialData = @(x) sin(pi * 1 * x / obj.xspan(2));

      % assembly.initialData = @(x) ones(size(x, 1), size(x, 2));

      tic
      LHS = assembly.assembleFieldIndependentMatrix();
      toc
      tic
      % O = assembly.assembleFieldDependentMatrixForFourierSeries(5);
      % O = assembly.assembleFieldDependentMatrixForSineSeriesSlow(1);
      toc;
      RHS = assembly.assembleVectorOnes();

      size(LHS)

      % compose field
      LHSO = LHS;
      % for idx = 1:length(cof)
      %   LHSO = LHSO + cof(idx) * O{idx};
      % end

      % nnz(LHS)
      % numel(LHS)
      % nnz(LHS) / numel(LHS)

      % spy(LHS)

      solfun = @(t, x) assembly.solutionFuncFromCoeffs(LHSO \ RHS, t, x);
      obj.plotSolution(solfun);
    end

    function solc = solveTwoFields(obj, lap, cof1, cof2, xg, tgl, tgr)
      N = 50;
      M = 50;
      C = 50;

      assembly = AssemblyFourierLegendre();
      assembly.setNumberOfAnsatzFuncs(N, M);
      assembly.setNumberOfTestFuncsFromAnsatzFuncs();
      assembly.xspan = obj.xspan;

      assembly.coeffLaplacian = lap;
      assembly.coeffOffset = 0;

      % solve for first field
      assembly.tspan = [obj.tspan(1), obj.fieldBreakpoint];
      assembly.initialData = @(x) ones(size(x, 1), size(x, 2));

      LHS1 = assembly.assembleFieldIndependentMatrix();
      O1 = assembly.assembleFieldDependentMatrixForFourierSeries(C);
      RHS1 = assembly.assembleRHS();

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
      RHS2 = assembly.assembleVectorFromCoeffs(sol1);

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

    function plotSolution(obj, solfun)
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
