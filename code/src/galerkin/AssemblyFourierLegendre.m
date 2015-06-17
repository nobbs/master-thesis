classdef AssemblyFourierLegendre < AssemblyGlobalAbstract
  % Assembly of the stiffness matrix and load vector for a Galerkin method using
  % global basis functions with periodic boundary conditions.
  %
  % We construct finite-dimensional subspaces for the ansatz space `\mathcal X`
  % and the test space `\mathcal Y` by using global Fourier basis functions
  % (sine and cosine) in space and global Legendre polynomials in time. By
  % orthogonality of the respective functions we gain a sparse stiffness matrix.
  %
  % Let
  %   ``s_j(x) = \begin{cases}
  %       \cos(\frac{2 \pi}{L} (j-1) x), & j = 1 \mod 2\\
  %       \sin(\frac{2 \pi}{L} j x), & j = 0 \mod 2
  %     \end{cases},``
  % then we can choose
  %  ``\mathcal X_{N_s, N_l} = \Set{ s_j(x) L_k(t) \given j = 1,
  %    \dots, N_s,~k = 0, \dots, N_l }``
  % for the finite-dimensional ansatz space `\mathcal X_{N_s, N_l}` and respectively
  % ``\mathcal Y_{M_s, M_l, M_i} = \Set{ (s_l(x) L_m(t), 0 )
  %     \given l = 1, \dots, M_s,~m = 0, \dots, M_l } \cup \Set{ (0,
  %     s_n(x)) \given n = 1, \dots, M_i}``
  % for the test space `\mathcal Y_{M_s, M_l, M_i}` with dimensions `\dim
  % \mathcal X_{N_s, N_l} = N_s  N_l` and `\dim
  % \mathcal Y_{M_s, M_l, M_i} = M_s M_l M_i`.
  %
  % See also:
  %   AssemblySineLegendre
  %
  % @todo Not yet fully implemented.

  properties
    coeffLaplace;
    coeffOffset;
    xspan;
    tspan;

    % Initial condition data
    % @see generateInitialDataStruct
    initialData;

    % External field `\omega \colon \Omega \to \mathbb{R}`
    field;
  end

  properties(Dependent)
    % Checks if the source term is non-zero @type logical
    isSourceNonZero;
  end

  methods

    % Constructor

    function obj = AssemblyFourierLegendre()

    end

    % Custom getters and setters

    function obj = set.xspan(obj, val)
      if size(val, 1) ~= 1 || size(val, 2) ~= 2
        error('Given value is not a row vector of size 2!');
      elseif val(1) ~= 0
        error('The spatial interval must start at 0!');
      elseif val(2) <= val(1)
        error('The given spatial interval is empty or wrongly oriented!');
      end

      obj.xspan = val;
    end

    function val = get.isSourceNonZero(obj)
      val = false;
    end

    function M = assembleStiffnessMatrixWithoutOmega(obj)
      % Assemble the field-independent part of the stiffness matrix.
      %
      % This is done by evaluating the individual summands of the left hand side
      % of the variational problem and by utilization of the orthogonality of
      % the chosen basis functions.
      %
      % Return values:
      %   M: Field-independent part of the stiffness matrix
      %     @type sparsematrix

      % Preparation for the sparse matrix
      Idx = ones(obj.dAnsatz, 1);
      Idy = ones(obj.dAnsatz, 1);
      Val = zeros(obj.dAnsatz, 1);
      ctr = 1;

      % Calculate the first term `\int_{I} \skp{u_t(t)}{v_1(t)}{L_2(\Omega)} \diff t`.
      for jdx = 1:min(obj.nAnsatzSpatial, obj.nTestSpatial)
        for kdx = 1:obj.nAnsatzTemporal
          for mdx = 1:kdx
            if kdx > mdx && mod(kdx + mdx, 2) == 1
              % evaluate spatial integral
              if jdx == 1
                intSpatial = obj.xspan(2);
              else
                intSpatial = obj.xspan(2) / 2;
              end

              % evaluate temporal integral
              intTemporal = 2 / (obj.tspan(2) - obj.tspan(1));

              Idx(ctr) = (jdx - 1) * obj.nAnsatzTemporal + kdx;
              Idy(ctr) = (jdx - 1) * obj.nTestTemporal + mdx;
              Val(ctr) = intSpatial * intTemporal;
              ctr = ctr + 1;
            end
          end
        end
      end

      % Calculate the second term `\int_{I} \skp{\grad u(t)}{\grad v_1(t)}{L_2(\Omega)} \diff t`.
      for jdx = 1:min(obj.nAnsatzSpatial, obj.nTestSpatial)
        for kdx = 1:min(obj.nAnsatzTemporal, obj.nTestTemporal)
          % evaluate spatial integral
          if jdx == 1
            intSpatial = 0;
          elseif mod(jdx, 2) == 0
            intSpatial = (2 / obj.xspan(2)) * (pi * jdx)^2;
          else
            intSpatial = (2 / obj.xspan(2)) * (pi * (jdx - 1))^2;
          end

          % evaluate temporal integral
          intTemporal = (obj.tspan(2) - obj.tspan(1)) / (2 * (kdx - 1) + 1);

          Idx(ctr) = (jdx - 1) * obj.nAnsatzTemporal + kdx;
          Idy(ctr) = (jdx - 1) * obj.nTestTemporal + kdx;
          Val(ctr) = obj.coeffLaplace * intTemporal * intSpatial;
          ctr = ctr + 1;
        end
      end

      % Calculate the third term `\int_{I} \mu \skp{u(t)}{v_1(t)}{L_2(\Omega)} \diff t`.
      for jdx = 1:min(obj.nAnsatzSpatial, obj.nTestSpatial)
        for kdx = 1:min(obj.nAnsatzTemporal, obj.nTestTemporal)
          % evaluate spatial integral
          if jdx == 1
            intSpatial = obj.xspan(2);
          else
            intSpatial = obj.xspan(2) / 2;
          end

          % evaluate temporal integral
          intTemporal = (obj.tspan(2) - obj.tspan(1)) / (2 * (kdx - 1) + 1);

          Idx(ctr) = (jdx - 1) * obj.nAnsatzTemporal + kdx;
          Idy(ctr) = (jdx - 1) * obj.nTestTemporal + kdx;
          Val(ctr) = obj.coeffOffset * intTemporal * intSpatial;;
          ctr = ctr + 1;
        end
      end

      % Calculate the fourth term `\skp{u(0)}{v_2}{L_2(\Omega)}`.
      for jdx = 1:min(obj.nAnsatzSpatial, obj.nTestSpatialIC)
        for kdx = 1:obj.nAnsatzTemporal
          % evaluate spatial integral
          if jdx == 1
            intSpatial = obj.xspan(2);
          else
            intSpatial = obj.xspan(2) / 2;
          end

          % evaluate temporal integral
          intTemporal = (-1)^(kdx - 1);

          Idx(ctr) = (jdx - 1) * obj.nAnsatzTemporal + kdx;
          Idy(ctr) = obj.nTestSpatial * obj.nTestTemporal + jdx;
          Val(ctr) = intTemporal * intSpatial;
          ctr = ctr + 1;
        end
      end

      % Finally assemble the stiffness matrix
      M = sparse(Idy, Idx, Val, obj.dTest, obj.dAnsatz);

    end

    function [O1, O2] = assembleStiffnessMatrixOnlyOmega(obj)
      % Assemble only the field dependet part of the stiffness matrix.
      % @todo Not yet implemented.
      error('Not yet implemented!');
    end

    % function M = assembleOmegaMatrixFrom

    function F = assembleRHS(obj)
      % Assemble the load vector.
      %
      % Return values:
      %   F: RHS of the linear equation system @type colvector
      %
      % @todo Not yet implemented!

      dimY = obj.dTest;
      F = zeros(dimY, 1);

      if obj.isSourceNonZero
        for ldx = 1:obj.nTestSine
          for mdx = 1:obj.nTestTemporal
            val = integral2(@(t, x) obj.sourceFunc(t, x) .* obj.spatialBasisFunc(ldx, x) .* temporalBasisFunc(mdx, t), obj.tspan(1), obj.tspan(2), obj.xspan(1), obj.xspan(2));
            pos = (ldx - 1) * obj.nTestTemporal + mdx;
            F(pos) = val;
          end
        end
      end

      for ndx = 1:obj.nTestSpatialIC
        val = integral(@(x) obj.initialData(x) .* obj.spatialBasisFunc(ndx, x), obj.xspan(1), obj.xspan(2));
        pos = obj.nTestSpatial * obj.nTestTemporal + ndx;
        F(pos) = val;
      end
    end

    function val = spatialBasisFunc(obj, index, x)
      % Spatial basis functions.
      %
      % Evaluates the spatial basis function, Fourier functions, with the given
      % index for the given values of x. Can be used to define function handles
      % and numerical integration.
      %
      % Parameters:
      %   index: index of the basis function
      %   x: values in which the function should be evaluated
      %
      % Return values:
      %   val: values of the basis function in x

      if index == 1
        val = ones(size(x, 1), size(x, 2));
      elseif mod(index, 2) == 0
        val = sin(pi * index * x / obj.xspan(2));
      else
        val = cos(pi * (index - 1) * x / obj.xspan(2));
      end
    end

    function val = temporalBasisFunc(obj, index, t)
      % Temporal basis functions.
      %
      % Evaluates the temporal basis function, shifted Legendre polynomials,
      % with the given index for the given values of x. Can be used to define
      % function handles and numerical integration.
      %
      % Parameters:
      %   index: index of the basis function
      %   x: values in which the function should be evaluated
      %
      % Return values:
      %   val: values of the basis function in x

      val = legendrePolynomial(t, index - 1, obj.tspan);
    end

  end



end
