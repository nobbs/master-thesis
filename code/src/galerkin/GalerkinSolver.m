classdef GalerkinSolver < handle
  % Galerkin-based Solver.
  %
  % @todo Not yet implemented.

  properties
    % time interval @type vector
    tspan = [0, 1];
    % spatial interval @type vector
    xspan = [0, 1];
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
    coeffLaplace = 1;
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

    function [ufun] = solve(obj)
      % @todo Not yet implemented!

      % assembly = AssemblySineLegendre();
      assembly = AssemblyFourierLegendre();
      assembly.setNumberOfAnsatzFuncs(40, 40);
      assembly.setNumberOfTestFuncsFromAnsatzFuncs();

      obj.tspan = [0 obj.fieldBreakpoint];
      assembly.tspan = obj.tspan;
      assembly.xspan = obj.xspan;

      assembly.coeffLaplace = 0.1;
      assembly.coeffOffset = 0;

      % assembly.initialData = @(x) sin(pi * 1 * x / obj.xspan(2));

      assembly.initialData = @(x) ones(size(x, 1), size(x, 2));

      tic
      LHS = assembly.assembleStiffnessMatrixWithoutOmega();
      toc
      tic
      O = assembly.assembleStiffnessMatrixOmegaFromFourier(5);
      % O = assembly.assembleStiffnessMatrixOmegaFromSineSlow(1);
      toc;
      RHS = assembly.assembleRHS();

      size(LHS)

      LHS = LHS + 1 * O{1} - 5 * O{5};

      % nnz(LHS)
      % numel(LHS)
      % nnz(LHS) / numel(LHS)

      % spy(LHS)

      solfun = @(t, x) assembly.solutionFunctionFromCoeffs(LHS \ RHS, t, x);
      obj.plotSolution(solfun);
    end

    function solveTwoConstantFields(obj, fieldA, fieldB)
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
