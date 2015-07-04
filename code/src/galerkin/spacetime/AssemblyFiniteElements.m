classdef AssemblyFiniteElements < SpaceTimeAssemblyAbstract

  properties
    % number of temporal gridpoints
    nK;

    nAnsatzSpatial;
    nTestSpatial;
    nTestSpatialIC;

    % temporal grid (for both ansatz and test space)
    tgrid;
  end

  properties(Dependent)
    nAnsatzDim;
    nTestDim;
  end

  methods
    function dim = get.nAnsatzDim(obj)
      % Calculate the dimension of the ansatz space.

      dim = obj.nAnsatzSpatial * obj.nK;
    end

    function dim = get.nTestDim(obj)
      % Calculate the dimension of the test space.

      dim = obj.nTestSpatial * obj.nK + obj.nTestSpatialIC;
    end

    function val = solutionFuncFromCoeffs(obj, solutionCoeffs, t, x)
      % Construct a function handle of the solution.
      %
      % Use it to define a solution function in (t, x) or evaluate the solution
      % directly for given t and x grids.
      %
      % Parameters:
      %   solutionCoeffs: solution vector of the linear system @type colvec
      %   t: temporal variable @type vector
      %   x: spatial variable @type vector
      %   tspan: custom temporal interval @type vector @default obj.tspan
      %
      % Return values:
      %   solfun: function handle of the solution function @type function_handle

      val = zeros(size(t, 1), size(t, 2));

      % precompute the spatial and temporal basis functions for the given grids
      spatialValues = cell(obj.nAnsatzSpatial, 1);
      for jdx = 1:obj.nAnsatzSpatial
        spatialValues{jdx} = obj.spatialBasisFunc(jdx, x);
      end
      temporalValues = cell(obj.nK, 1);
      for kdx = 1:obj.nK
        temporalValues{kdx} = obj.temporalBasisFunc(kdx, t);
      end

      % we evaluate the solution in two steps. the inner loop adds up all the
      % temporal evaluations that share the same spatial basis function as a
      % factor. the outer loop then multiplies this with the spatial component.
      for jdx = 1:obj.nAnsatzSpatial
        temporalVal = zeros(size(t, 1), size(t, 2));
        for kdx = 1:obj.nK
          % Get the right coefficient
          pos = (kdx - 1) * obj.nK + jdx;

          % evaluate the corresponding basis functions
          temporalVal = temporalVal + solutionCoeffs(pos) * temporalValues{kdx};
        end
        val = val + temporalVal .* spatialValues{jdx};
      end
    end

    function matrices = assembleFieldIndependentMatrix(obj)
      % Assemble the field-independent part of the stiffness matrix.
      %
      % This is done by evaluating the individual summands of the left hand side
      % of the variational problem and by utilization of the orthogonality of
      % the chosen basis functions.
      %
      % Return values:
      %   M1: Time-Derivative part, `\int_{I} \skp{u_{t}(t)}{v_1(t)}{L_2(\Omega)}
      %     \diff t` @type sparsematrix
      %   M2: Laplacian `\int_{I} \skp{\nabla u(t)}{\nabla v_1(t)}{L_2(\Omega)}
      %     \diff t` @type sparsematrix
      %   M3: Offset part `\int_{I} \skp{u(t)}{v_1(t)}{L_2(\Omega)} \diff t`
      %     @type sparsematrix
      %   M4: Initial condition part `\skp{u(0)}{v_2}{L_2(\Omega)}` @type
      %     sparsematrix

      % create the needed cell array
      matrices = {};

      % precompute the needed block matrices
      MtF = TemporalAssemblyLinearConstant.massMatrix(obj.nK, obj.tgrid);
      CtF = TemporalAssemblyLinearConstant.halfStiffnessMatrix(obj.nK);
      eFtF = TemporalAssemblyLinearConstant.forwardInitVector(obj.nK);
      eBtF = TemporalAssemblyLinearConstant.backwardInitVector(obj.nK);

      % spatial matrices
      MxF = SpatialAssemblyFourier.massMatrix(obj.nAnsatzSpatial, obj.nTestSpatial, obj.xspan(2) - obj.xspan(1));
      AxF = SpatialAssemblyFourier.stiffnessMatrix(obj.nAnsatzSpatial, obj.nTestSpatial, obj.xspan(2) - obj.xspan(1));

      % Calculate the first term `\int_{I} \skp{u_t(t)}{v_1(t)}{L_2(\Omega)}
      % \diff t`.
      matrices.TiD = kron(CtF, MxF);
      % resize it to the right dimension
      matrices.TiD(obj.nTestDim, obj.nAnsatzDim) = 0;

      % Calculate the second term `\int_{I} \skp{\grad u(t)}{\grad
      % v_1(t)}{L_2(\Omega)} \diff t`.
      matrices.Lap = kron(MtF, AxF);
      % resize it to the right dimension
      matrices.Lap(obj.nTestDim, obj.nAnsatzDim) = 0;

      % Calculate the third term `\int_{I} \mu \skp{u(t)}{v_1(t)}{L_2(\Omega)}
      % \diff t`.
      matrices.Off = kron(MtF, MxF);
      % resize it to the right dimension
      matrices.Off(obj.nTestDim, obj.nAnsatzDim) = 0;

      % Calculate the fourth term `\skp{u(0)}{v_2}{L_2(\Omega)}`.
      MxF2 = SpatialAssemblyFourier.massMatrix(obj.nAnsatzSpatial, obj.nTestSpatialIC, obj.xspan(2) - obj.xspan(1));
      matrixForward = kron(eFtF, MxF2);
      matrixBackward = kron(eBtF, MxF2);
      % create the needed "big" sparse matrices
      matrices.ICF = sparse(obj.nTestDim, obj.nAnsatzDim);
      matrices.ICB = sparse(obj.nTestDim, obj.nAnsatzDim);
      % place it in the right place
      matrices.ICF((obj.nK * obj.nTestSpatial + 1):end, :) = matrixForward;
      matrices.ICB((obj.nK * obj.nTestSpatial + 1):end, :) = matrixBackward;

    end % assembleFieldIndependentMatrix

    function F = assembleVectorOnes(obj)
      % Assemble the load vector for no source and initial data equal to a
      % constant value of one.
      %
      % Return values:
      %   F: load vector @type colvec

      F = zeros(obj.nTestDim, 1);
      F(obj.nTestSpatial * obj.nK + 1) = obj.xspan(2);

      % normalize if needed
      % if obj.useNormalization
      %   F = obj.TestNormDiag \ F;
      % end
    end % assembleVectorOnes
  end

  methods

    function val = spatialBasisFunc(obj, index, x)
      % Spatial basis functions.
      %
      % Evaluates the spatial basis function, Fourier functions, with the given
      % index for the given values of x. Can be used to define function handles
      % and numerical integration.
      %
      % In this case the spatial basis functions are of the type
      % ``\begin{cases}
      %   1, & i = 1\\
      %   \sin(\pi i x / L), & i~\text{even}\\
      %   \cos(\pi i x / L), & i~\text{odd and}~i > 1
      % \end{cases},`` where `i` corresponds to index and `L` is the width of
      % the spatial interval.
      %
      % Parameters:
      %   index: index of the basis function @type integer
      %   x: values in which the function should be evaluated @type matrix
      %
      % Return values:
      %   val: values of the basis function in x @type matrix

      if index == 1
        val = ones(size(x, 1), size(x, 2));
      elseif mod(index, 2) == 0
        val = sin(pi * index * x / obj.xspan(2));
      else
        val = cos(pi * (index - 1) * x / obj.xspan(2));
      end
    end % spatialBasisFunc


    function val = temporalBasisFunc(obj, index, t)
      % Temporal basis functions.
      %
      % Evaluates the temporal basis function, shifted Legendre polynomials,
      % with the given index for the given values of t. Can be used to define
      % function handles and numerical integration.
      %
      % Warning:
      %   The enumeration of index starts at 1, that means you'll get the
      %   Legendre polynomial of degree (index - 1)!
      %
      % Parameters:
      %   index: index of the basis function @type integer
      %   t: values in which the function should be evaluated @type matrix
      %   tspan: custom temporal interval @type vector @default obj.tspan
      %
      % Return values:
      %   val: values of the basis function in t @type matrix

      if index == 1
        val = (obj.tspan(1) <= t & t < obj.tgrid(2)) .* (obj.tgrid(2) - t) / (obj.tgrid(2) - obj.tgrid(1));
      elseif index == obj.nK
        val = (obj.tgrid(end - 1) <= t & t <= obj.tgrid(end)) .* (t - obj.tgrid(end - 1)) / (obj.tgrid(end) - obj.tgrid(end - 1))
      else
        val = (obj.tgrid(index - 1) <= t & t < obj.tgrid(index)) .* (t - obj.tgrid(index - 1)) / (obj.tgrid(index) - obj.tgrid(index - 1)) + ...
          (obj.tgrid(index) <= t & t < obj.tgrid(index + 1)) .* (obj.tgrid(index + 1) - t) / (obj.tgrid(index + 1) - obj.tgrid(index))
      end
    end % spatialBasisFuncDerivative

  end

end % classdef
