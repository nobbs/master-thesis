classdef AssemblySineLegendre < AssemblyAbstract
  % Assembly of the stiffness matrix and load vector for a Galerkin method using
  % global basis functions with homogeneous boundary conditions.
  %
  % We construct finite-dimensional subspaces for the ansatz space `\mathcal X`
  % and the test space `\mathcal Y` by using global sine-functions in space and
  % global Legendre polynomials in time. By orthogonality of the respective
  % functions we gain a sparse stiffness matrix.
  %
  % For the finite-dimensional ansatz space `\mathcal X_{N_s, N_l}` we choose
  %  ``\mathcal X_{N_s, N_l} = \Set{ \sin(\frac{\pi j}{L} x) L_k(t) \given j = 1,
  %    \dots, N_s,~k = 0, \dots, N_l }``
  % and respectively for the test space `\mathcal Y_{M_s, M_l, M_i}`
  %  ``\mathcal Y_{M_s, M_l, M_i} = \Set{ (\sin(\frac{\pi l}{L} x) L_m(t), 0 )
  %     \given l = 1, \dots, M_s,~m = 0, \dots, M_l } \cup \Set{ (0,
  %     \sin(\frac{\pi n}{L} x)) \given n = 1, \dots, M_i}.``
  % Obviously it holds that `\dim \mathcal X_{N_s, N_l} = N_s  N_l` and `\dim
  % \mathcal Y_{M_s, M_l, M_i} = M_s M_l M_i`.
  %
  % See also:
  %   AssemblyAbstract
  %
  % @todo Not yet fully implemented.

  properties
    coeffLaplace;
    coeffOffset;
    xspan;
    tspan;

    % Initial condition data
    initialData;

    % Kind of initial condition data. Possible values are
    %   'func': funtion handle
    %   'sineCoeffs': coefficients of a sine series
    %   'vector': values of the function on a equidistant grid
    initialDataKind = 'func';

    % External field `\omega \colon \Omega \to \mathbb{R}`
    field;

    % Number of sine basis functions for the ansatz subspace @type integer
    nAnsatzSine;

    % Number of Legendre basis polynomials for the ansatz subspace @type integer
    nAnsatzLegendre;

    % Number of sine basis functions for the test subspace @type integer
    nTestSine;

    % Number of Legendre basis polynomials for the test subspace @type integer
    nTestLegendre;

    % Number of sine basis functions for the initial condition part of the test
    % subspace @type integer
    nTestSineIC;
  end

  properties(Dependent)
    % Dimension of the ansatz subspace @type integer
    dAnsatz;

    % Dimension of the test subspace @type integer
    dTest;

    % Checks if the source term is non-zero @type logical
    isSourceNonZero;
  end

  methods

    % Constructor

    function obj = AssemblySineLegendre()

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

    function dim = get.dAnsatz(obj)
      dim = obj.nAnsatzSine * obj.nAnsatzLegendre;
    end

    function dim = get.dTest(obj)
      dim = obj.nTestSine * obj.nTestLegendre + obj.nTestSineIC;
    end

    function val = get.isSourceNonZero(obj)
      val = false;
    end

    function obj = setNumberOfAnsatzFuncs(obj, nAnsatzSine, nAnsatzLegendre)
      % Set the number of basis functions for the ansatz subspace.
      %
      % Parameters:
      %   nAnsatzSine: Number of sine basis functions @type integer
      %   nAnsatzLegendre: Number of Legendre basis polynomials @type integer

      obj.nAnsatzSine     = nAnsatzSine;
      obj.nAnsatzLegendre = nAnsatzLegendre;
    end

    function setNumberOfTestFuncs(obj, nTestSine, nTestLegendre, nTestSineIC)
      % Set the number of basis functions for the test subspace.
      %
      % Parameters:
      %   nTestSine: Number of sine basis functions  @type integer
      %   nTestLegendre: Number of Legendre basis polynomials @type integer
      %   nTestSineIC: Number of sine basis functions for the initial condition
      %     @type integer

      obj.nTestSine     = nTestSine;
      obj.nTestLegendre = nTestLegendre;
      obj.nTestSineIC   = nTestSineIC;
    end

    function setNumberOfTestFuncsFromAnsatzFuncs(obj)
      % Set the number of basis functions for the test subspace according to the
      % number of basis functions for the subspace.
      %
      % This guaranties a quadratic system of linear equations.

      obj.nTestSine     = obj.nAnsatzSine;
      obj.nTestLegendre = obj.nAnsatzLegendre - 1;
      obj.nTestSineIC   = obj.nAnsatzSine;
    end

    function StiffnessMatrix = assembleStiffnessMatrixWithoutOmega(obj)
      % Assemble the field-independent part of the stiffness matrix.
      %
      % This is done by evaluating the individual summands of the left hand side
      % of the variational problem and by utilization of the orthogonality of
      % the chosen basis functions.
      %
      % Return values:
      %   StiffnessMatrix: Field-independet part of the stiffness matrix
      %     @type sparsematrix

      % Preparation for the sparse matrix
      Idx = ones(obj.dAnsatz, 1);
      Idy = ones(obj.dAnsatz, 1);
      Val = zeros(obj.dAnsatz, 1);
      ctr = 1;

      % Calculate the first term `\int_{I} \skp{u_t(t)}{v_1(t)}{L_2(\Omega)} \diff t`.
      for jdx = 1:min(obj.nAnsatzSine, obj.nTestSine)
        for kdx = 1:obj.nAnsatzLegendre
          for mdx = 1:kdx
            if kdx > mdx && mod(kdx + mdx, 2) == 1
              Idx(ctr) = (jdx - 1) * obj.nAnsatzLegendre + kdx;
              Idy(ctr) = (jdx - 1) * obj.nTestLegendre + mdx;
              Val(ctr) = (obj.xspan(2) / 2) * 2 / (obj.tspan(2) - obj.tspan(1));
              ctr = ctr + 1;
            end
          end
        end
      end

      % Berechnet den zweiten Summanden der LHS des Variationsproblems.
      %
      % Der zweite Summand (und erster Term der Bilinearform `a`), das
      % heißt ``\int_{I} \int_{\Omega} c \nabla u(t, x) \nabla v_{1}(t,x)
      % dx dt,`` wird ausgewertet.
      %
      % Im homogenen Fall vereinfacht sich dies erneut und kann aufgrund
      % der Orthogonalität der Basisfunktionen ohne numerische Quadratur
      % berechnet werden.
      %
      % Return values:
      %   M: Zweiter Summand der Massematrix @type sparsematrix

      % Calculate the second term `\int_{I} \skp{\grad u(t)}{\grad v_1(t)}{L_2(\Omega)} \diff t`.
      for jdx = 1:min(obj.nAnsatzSine, obj.nTestSine)
        for kdx = 1:min(obj.nAnsatzLegendre, obj.nTestLegendre)
          kk = kdx - 1;

          val = obj.coeffLaplace * (obj.tspan(2) - obj.tspan(1)) / (2 * kk + 1) * (pi * jdx)^2 / (2 * obj.xspan(2));

          Idx(ctr) = (jdx - 1) * obj.nAnsatzLegendre + kdx;
          Idy(ctr) = (jdx - 1) * obj.nTestLegendre + kdx;
          Val(ctr) = val;
          ctr = ctr + 1;
        end
      end

      % Berechnet den dritten Summanden der LHS des Variationsproblems.
      %
      % Der dritte Summand (und zweiter Term der Bilinearform `a`), das
      % heißt ``\int_{I} \int_{\Omega} \mu  u(t, x) v_{1}(t,x) dx dt,``
      % wird ausgewertet.
      %
      % Im homogenen Fall vereinfacht sich dies wieder und kann aufgrund
      % der Orthogonalität der Basisfunktionen ohne numerische Quadratur
      % berechnet werden.
      %
      % Return values:
      %   M: Dritter Summand der Massematrix @type sparsematrix

      % Calculate the third term `\int_{I} \mu \skp{u(t)}{v_1(t)}{L_2(\Omega)} \diff t`.
      for jdx = 1:min(obj.nAnsatzSine, obj.nTestSine)
        for kdx = 1:min(obj.nAnsatzLegendre, obj.nTestLegendre)
          kk = kdx - 1;
          val = obj.coeffOffset * (obj.tspan(2) - obj.tspan(1)) / (2 * kk + 1) * (obj.xspan(2) / 2);

          if val ~= 0
            x_pos = (jdx - 1) * obj.nAnsatzLegendre + kdx;
            y_pos = (jdx - 1) * obj.nTestLegendre + kdx;

            Idx(ctr) = x_pos;
            Idy(ctr) = y_pos;
            Val(ctr) = val;
            ctr = ctr + 1;
          end
        end
      end

      % Berechnet den vierten Summanden der LHS des Variationsproblems.
      %
      % Der vierte Summand, die Anfangsbedingung, das heißt
      % ``\int_{\Omega} c \nabla u(0, x) \nabla v_{2}(x) dx dt,`` wird
      % ausgewertet.
      %
      % Im homogenen Fall vereinfacht sich dies erneut und kann aufgrund
      % der Orthogonalität der Basisfunktionen ohne numerische Quadratur
      % berechnet werden.
      %
      % Return values:
      %   M: Vierter Summand der Massematrix @type sparsematrix

      % Calculate the fourth term `\skp{u(0)}{v_2}{L_2(\Omega)}`.
      for jdx = 1:min(obj.nAnsatzSine, obj.nTestSineIC)
        for kdx = 1:obj.nAnsatzLegendre
          kk = kdx - 1;

          val = (-1)^(kk) * obj.xspan(2) / 2;

          Idx(ctr) = (jdx - 1) * obj.nAnsatzLegendre + kdx;
          Idy(ctr) = obj.nTestSine * obj.nTestLegendre + jdx;
          Val(ctr) = val;
          ctr = ctr + 1;
        end
      end

      % Finally assemble the stiffness matrix
      StiffnessMatrix = sparse(Idy, Idx, Val, obj.dTest, obj.dAnsatz);
    end

    function [O1, O2] = assembleStiffnessMatrixOnlyOmega(obj)
      % Assemble only the field dependet part of the stiffness matrix.
      % @todo Not yet implemented.
      error('Not yet implemented!');
    end

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
          for mdx = 1:obj.nTestLegendre
            mm = mdx - 1;
            val = integral2(@(t, x) obj.sourceFunc(t, x) .* sin(pi * ldx * x / obj.xspan(2)) .* legendrePolynomial(t, mm), obj.tspan(1), obj.tspan(2), obj.xspan(1), obj.xspan(2));
            pos = (ldx - 1) * obj.nTestLegendre + mdx;
            F(pos) = val;
          end
        end
      end

      if obj.initialDataKind == 'func'
        for ndx = 1:obj.nTestSineIC
          val = integral(@(x) obj.initialData(x) .* sin(pi * ndx * x / obj.xspan(2)), obj.xspan(1), obj.xspan(2));
          pos = obj.nTestSine * obj.nTestLegendre + ndx;
          F(pos) = val;
        end
      else
        error('Not yet implemented!');
      end
    end

    function val = solutionFunctionFromCoeffs(obj, solutionCoeffs, t, x)
      % Construct a function handle of the solution.
      %
      % Use it to define a solution function in (t, x) or evaluate the solution
      % directly for given t and x grids.
      %
      % Parameters:
      %   solutionCoeffs: solution vector of the linear system @type colvec
      %   t: time variable
      %   x: spatial variable
      %
      % Return values:
      %   solfun: function handle of the solution function @type function_handle
      %
      % @todo Not yet fully implemented.

      val = zeros(size(t, 1), size(t, 2));

      for jdx = 1:obj.nAnsatzSine
          for kdx = 1:obj.nAnsatzLegendre
            % Get the right coefficient
            pos = (jdx - 1) * obj.nAnsatzLegendre + kdx;
            % normalization constant
            % cnrm = sqrt((1 + (pi * j)^2) / (2 * (2*(k - 1) + 1)) + ...
                      % legendre_dP(1, k - 1));
            % evaluate the corresponding basis functions
            val = val + solutionCoeffs(pos) * sin(pi * jdx * x / obj.xspan(2)) .* legendrePolynomial(t, kdx - 1, obj.tspan);
          end
      end
    end

  end

  methods(Access = 'private')

  end

end
