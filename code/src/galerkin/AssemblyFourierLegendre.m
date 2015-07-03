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

  properties
    % nothing to see here
  end % properties

  methods

    %% Fast implementation of the assembly methods

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
      MtF  = obj.assembleTemporalMassMatrix(obj.nAnsatzTemporal, obj.nTestTemporal);
      CtF  = obj.assembleTemporalHalfStiffnessMatrix(obj.nAnsatzTemporal, obj.nTestTemporal);
      eFtF = obj.assembleTemporalInitForwardVector(obj.nAnsatzTemporal);
      eBtF = obj.assembleTemporalInitBackwardVector(obj.nAnsatzTemporal);

      % spatial matrices
      MxF = obj.assembleSpatialMassMatrix(obj.nAnsatzSpatial, obj.nTestSpatial);
      AxF = obj.assembleSpatialStiffnessMatrix(obj.nAnsatzSpatial, obj.nTestSpatial);

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
      MxF2 = obj.assembleSpatialMassMatrix(obj.nAnsatzSpatial, obj.nTestSpatialIC);
      matrixForward = kron(eFtF, MxF2);
      matrixBackward = kron(eBtF, MxF2);
      % create the needed "big" sparse matrices
      matrices.ICF = sparse(obj.nTestDim, obj.nAnsatzDim);
      matrices.ICB = sparse(obj.nTestDim, obj.nAnsatzDim);
      % place it in the right place
      matrices.ICF((obj.nTestTemporal * obj.nTestSpatial + 1):end, :) = matrixForward;
      matrices.ICB((obj.nTestTemporal * obj.nTestSpatial + 1):end, :) = matrixBackward;

      % normalize if needed
      if obj.useNormalization
        matrices.TiD = (obj.TestNormDiag \ matrices.TiD) / obj.AnsatzNormDiag;
        matrices.Lap = (obj.TestNormDiag \ matrices.Lap) / obj.AnsatzNormDiag;
        matrices.Off = (obj.TestNormDiag \ matrices.Off) / obj.AnsatzNormDiag;
        matrices.ICF = (obj.TestNormDiag \ matrices.ICF) / obj.AnsatzNormDiag;
        matrices.ICB = (obj.TestNormDiag \ matrices.ICB) / obj.AnsatzNormDiag;
      end
    end % assembleFieldIndependentMatrix


    %% Assemble the field-dependent parts of the stiffness matrix based on
    %  different series expansions of the field.

    function O = assembleFieldDependentMatrixForSineSeries(obj, nCoeff)
      % Assemble only the field-dependent components of the stiffness matrix
      % based on a sine series expansion of the field. Optimized version without
      % numerical integration.
      %
      % Integrals over the product of three basis functions of the different
      % spaces are evaluated through formulas which can mostly be derived
      % through the use of trigonometric identities, partial integration and
      % some intuition.
      %
      % Parameters:
      %   nCoeff: number of sine basis functions @type integer
      %
      % Return values:
      %   O: cellarray of the partial stiffness matrices @type struct

      % create the cellarray that will hold the matrices
      O = cell(nCoeff, 1);

      MtF = obj.assembleTemporalMassMatrix(obj.nAnsatzTemporal, obj.nTestTemporal);

      % iterate over the index of the sine basis functions (of the sine series
      % expansion of the field).
      for cdx = 1:nCoeff
        % create needed vectors to assemble the sparse matrix
        Idx = ones(obj.nAnsatzDim, 1);
        Idy = ones(obj.nAnsatzDim, 1);
        Val = zeros(obj.nAnsatzDim, 1);
        ctr = 1;

        % iterate over the indexes of the spatial basis functions of the ansatz
        % and test subspaces
        for jdx = 1:obj.nAnsatzSpatial
          for ldx = 1:obj.nTestSpatial
            % evaluate spatial integral
            % @todo explain optimization and cases
            intSpatial = 0;
            if mod(jdx, 2) == 1 && mod(ldx, 2) == 1 && mod(cdx, 2) == 1
              intSpatial = obj.xspan(2) * (cdx / (pi * (cdx^2 - (jdx - ldx)^2)) + cdx / (pi * (cdx^2 - (jdx + ldx - 2)^2)));
            elseif mod(jdx, 2) == 0 && mod(ldx, 2) == 0 && mod(cdx, 2) == 1
              intSpatial = obj.xspan(2) * (cdx / (pi * (cdx^2 - (jdx - ldx)^2)) - cdx / (pi * (cdx^2 - (jdx + ldx)^2)));
            elseif mod(jdx, 2) == 0 && mod(ldx, 2) == 1 && mod(cdx, 2) == 0
              val = 0;
              if cdx + jdx - ldx + 1 == 0
                val = val - 1;
              end
              if - cdx + jdx + ldx - 1 == 0
                val = val + 1;
              end
              if cdx - jdx + ldx - 1 == 0
                val = val + 1;
              end
              if cdx + jdx + ldx - 1 == 0
                val = val - 1;
              end
              intSpatial = (obj.xspan(2) / 4) * val;
            elseif mod(jdx, 2) == 1 && mod(ldx, 2) == 0 && mod(cdx, 2) == 0
              val = 0;
              if cdx - jdx + ldx + 1 == 0
                val = val - 1;
              end
              if - cdx + jdx + ldx - 1 == 0
                val = val + 1;
              end
              if cdx + jdx - ldx - 1 == 0
                val = val + 1;
              end
              if cdx + jdx + ldx - 1 == 0
                val = val - 1;
              end
              intSpatial = (obj.xspan(2) / 4) * val;
            end

            Idx(ctr) = jdx;
            Idy(ctr) = ldx;
            Val(ctr) = intSpatial;
            ctr = ctr + 1;
          end
        end

        % create the sparse matrix
        O{cdx} = kron(MtF, sparse(Idy, Idx, Val, obj.nTestSpatial, obj.nAnsatzSpatial));
        O{cdx}(obj.nTestDim, obj.nAnsatzDim) = 0;

        % normalize if needed
        if obj.useNormalization
          O{cdx} = (obj.TestNormDiag \ O{cdx}) / obj.AnsatzNormDiag;
        end
      end
    end % assembleFieldDependentMatrixForSineSeries

    function O = assembleFieldDependentMatrixForFourierSeries(obj, nCoeff)
      % Assemble only the field-dependent components of the stiffness matrix
      % based on a sine series expansion of the field. Optimized version without
      % numerical integration.
      %
      % Integrals over the product of three basis functions of the different
      % spaces are evaluated through formulas which can mostly be derived
      % through the use of trigonometric identities, partial integration and
      % some intuition.
      %
      % Parameters:
      %   nCoeff: number of sine basis functions @type integer
      %
      % Return values:
      %   O: cellarray of the partial stiffness matrices @type struct

      % create the cellarray that will hold the matrices
      O = cell(nCoeff, 1);

      MtF = obj.assembleTemporalMassMatrix(obj.nAnsatzTemporal, obj.nTestTemporal);

      % iterate over the index of the Fourier basis function (of the Fourier
      % series expansion of the field).
      % order of functions:
      % 1: cos(2pix), 2: sin(2pix), 3: cos(4pix) ...
      % cdx = i odd: cos((i+1) pi x), cdx = i even: sin(i pi x)
      for cdx = 1:nCoeff
        % create needed vectors to assemble the sparse matrix
        Idx = ones(obj.nAnsatzDim, 1);
        Idy = ones(obj.nAnsatzDim, 1);
        Val = zeros(obj.nAnsatzDim, 1);
        ctr = 1;

        % iterate over the indexes of the spatial basis functions of the ansatz
        % and test subspaces
        for jdx = 1:obj.nAnsatzSpatial
          for ldx = 1:obj.nTestSpatial

            % evaluate the spatial integral. there are several combinations of
            % sine / cosine products we have to consider
            intSpatial = 0;

            if mod(cdx, 2) == 0
              % cdx even: sine function sin(pi * cdx * x)
              if mod(jdx, 2) == 1 && mod(ldx, 2) == 1
                % j odd: cosine, l odd: cosine
                % intSpatial = 0;
              elseif mod(jdx, 2) == 0 && mod(ldx, 2) == 0
                % j even: sine, l even: sine
                % intSpatial = 0;
              elseif mod(jdx, 2) == 0 && mod(ldx, 2) == 1
                % j even: sine, l odd: cosine
                val = 0;
                if cdx + jdx - ldx + 1 == 0
                  val = val - 1;
                end
                if - cdx + jdx + ldx - 1 == 0
                  val = val + 1;
                end
                if cdx - jdx + ldx - 1 == 0
                  val = val + 1;
                end
                if cdx + jdx + ldx - 1 == 0
                  val = val - 1;
                end
                intSpatial = (obj.xspan(2) / 4) * val;
              elseif mod(jdx, 2) == 1 && mod(ldx, 2) == 0
                % j odd: cosine, l even: sine
                val = 0;
                if cdx + jdx - ldx - 1 == 0
                  val = val + 1;
                end
                if - cdx + jdx + ldx - 1 == 0
                  val = val + 1;
                end
                if cdx - jdx + ldx + 1 == 0
                  val = val - 1;
                end
                if cdx + jdx + ldx - 1 == 0
                  val = val - 1;
                end
                intSpatial = (obj.xspan(2) / 4) * val;
              end
            else
              % cdx odd: cosine function cos(pi * (cdx + 1) * x)
              if mod(jdx, 2) == 1 && mod(ldx, 2) == 1
                % j odd: cosine, l odd: cosine
                if jdx == 1 && ldx == 1
                  intSpatial = 0;
                elseif jdx == 1 && ldx > 1 && cdx + 1 == ldx - 1
                  intSpatial = obj.xspan(2) / 2;
                elseif jdx > 1 && ldx == 1 && cdx + 1 == jdx - 1
                  intSpatial = obj.xspan(2) / 2;
                elseif jdx > 1 && ldx > 1
                  val = 0;
                  if cdx + jdx - ldx + 1 == 0
                    val = val + 1;
                  end
                  if - cdx + jdx + ldx - 3 == 0
                    val = val + 1;
                  end
                  if cdx - jdx + ldx + 1 == 0
                    val = val + 1;
                  end
                  if cdx + jdx + ldx - 1 == 0
                    val = val + 1;
                  end
                  intSpatial = (obj.xspan(2) / 4) * val;
                end
              elseif mod(jdx, 2) == 0 && mod(ldx, 2) == 0
                % j even: sine, l even: sine
                val = 0;
                if - cdx + jdx + ldx - 1 == 0
                  val = val - 1;
                end
                if cdx + jdx - ldx + 1 == 0
                  val = val + 1;
                end
                if cdx - jdx + ldx + 1 == 0
                  val = val + 1;
                end
                if cdx + jdx + ldx + 1 == 0
                  val = val - 1;
                end
                intSpatial = (obj.xspan(2) / 4) * val;
              elseif mod(jdx, 2) == 0 && mod(ldx, 2) == 1
                % j even: sine, l odd: cosine
                % intSpatial = 0;
              elseif mod(jdx, 2) == 1 && mod(ldx, 2) == 0
                % j odd: cosine, l even: sine
                % intSpatial = 0;
              end
            end

            Idx(ctr) = jdx;
            Idy(ctr) = ldx;
            Val(ctr) = intSpatial;
            ctr = ctr + 1;
          end
        end

        % create the sparse matrix
        O{cdx} = kron(MtF, sparse(Idy, Idx, Val, obj.nTestSpatial, obj.nAnsatzSpatial));
        O{cdx}(obj.nTestDim, obj.nAnsatzDim) = 0;

        % normalize if needed
        if obj.useNormalization
          O{cdx} = (obj.TestNormDiag \ O{cdx}) / obj.AnsatzNormDiag;
        end
      end
    end % assembleFieldDependentMatrixForFourierSeries


    %% Assembly of the load vector for different kinds of data

    function F = assembleVectorFromSpatialCoeffs(obj, coeffs)
      % Assemble the load vector for no source and initial data given by
      % coefficients corresponding to the spatial basis functions.
      %
      % Parameters:
      %   coeffs: coefficients corresponding to spatial basis functions
      %     @type vector
      %
      % Return values:
      %   F: load vector @type colvec

      F = zeros(obj.nTestDim, 1);

      % iterate over the spatial basis functions corresponding to the initial
      % condition
      for ndx = 1:min(obj.nTestSpatialIC, obj.nAnsatzSpatial)
        if ndx == 1
          % corresponds to the constant basis function
          val = coeffs(ndx) * obj.xspan(2);
        else
          % corresponds to the sine and cosine basis functions
          val = coeffs(ndx) * obj.xspan(2) / 2;
        end

        pos = obj.nTestSpatial * obj.nTestTemporal + ndx;
        F(pos) = val;
      end

      % normalize if needed
      if obj.useNormalization
        F = obj.TestNormDiag \ F;
      end
    end % assembleVectorFromSpatialCoeffs

    function F = assembleVectorFromSolutionCoeffs(obj, solCoeffs, backward)
      % Assemble the load vector for no source and initial data given by
      % coefficients corresponding to the spatial basis and temporal basis
      % functions at the time endpoint.
      %
      % Parameters:
      %   solCoeffs: coefficients corresponding to spatial basis functions
      %     @type vector
      %   backward: toggle between forward or backward propagator @type logical
      %     @default false
      %
      % Return values:
      %   F: load vector @type colvec

      % set default values for optional arguments
      if nargin == 2
        backward = false;
      end

      sumcoeffs = zeros(obj.nAnsatzSpatial, 1);
      Kdx = 1:obj.nAnsatzTemporal;

      % forward and backward propagator have a small difference
      if ~backward
        % just plain sum up the coefficients for the different temporal indices,
        % because the temporal basis function are alwas one at the endpoint.
        for idx = 1:obj.nAnsatzSpatial
          % sumcoeffs(idx) = sum(solCoeffs(((idx - 1) * ...
          %   obj.nAnsatzTemporal) + Kdx));
          sumcoeffs(idx) = sum(solCoeffs(((Kdx - 1) * ...
            obj.nAnsatzSpatial) + idx));
        end
      else
        % sum with alternating sign, because the temporal basis function with
        % index i has the value (-1)^i at the start point.
        for idx = 1:obj.nAnsatzSpatial
          sumcoeffs(idx) = sum((-1).^(idx - 1) * ...
            solCoeffs(((Kdx - 1) * obj.nAnsatzSpatial) + idx));
        end
      end

      % and now: create the right hand side
      F = obj.assembleVectorFromSpatialCoeffs(sumcoeffs);
    end % assembleVectorFromSolutionCoeffs

    function F = assembleVectorOnes(obj)
      % Assemble the load vector for no source and initial data equal to a
      % constant value of one.
      %
      % Return values:
      %   F: load vector @type colvec

      F = zeros(obj.nTestDim, 1);
      F(obj.nTestSpatial * obj.nTestTemporal + 1) = obj.xspan(2);

      % normalize if needed
      if obj.useNormalization
        F = obj.TestNormDiag \ F;
      end
    end % assembleVectorOnes


    %% Assembly methods for the discrete subspace norms

    function M = assembleAnsatzNormMatrix(obj, useNormalization)
      % Assemble the mass matrix of the discrete norm on the ansatz space.
      %
      % Fast implementation that leverages a lot of simplification by hand.
      %
      % Parameters:
      %   useNormalization: toggle if the ansatz functions should be normalized.
      %     you have to set both the class property useNormalization and this
      %     parameter to true to get a normalized matrix. @type logical @default
      %     true
      %
      % Return values:
      %   M: mass matrix @type sparsematrix

      % set default values
      if nargin == 1
        useNormalization = obj.useNormalization;
      end

      % precompute the needed block matrices
      MtF = obj.assembleTemporalMassMatrix(obj.nAnsatzTemporal, obj.nAnsatzTemporal);
      AtF = obj.assembleTemporalStiffnessMatrix(obj.nAnsatzTemporal);

      % spatial matrices
      MxF = obj.assembleSpatialMassMatrix(obj.nAnsatzSpatial, obj.nAnsatzSpatial);
      AxF = obj.assembleSpatialStiffnessMatrix(obj.nAnsatzSpatial, obj.nAnsatzSpatial);

      % first part `\norm{u}_{L_2(I; V)}`: iterate over spatial and temporal
      % basis functions; respectively only in one dimension because boths sets
      % of basis functions are orthogonal in itself.
      mat1 = kron(MtF, MxF + AxF);

      % second part `norm{u_t}_{L_2(I; V')}`: again iterate over spatial and
      % temporal basis functions; this time with two indexes for the temporal
      % component, as the first derivatives of the temporal basis functions are
      % no longer orthogonal
      mat2 = kron(AtF, MxF);

      % Assemble the sparse mass matrix (with normalization, if needed)
      if ~useNormalization || ~obj.useNormalization
        M = mat1 + mat2;
      else
        M = (obj.AnsatzNormDiag \ (mat1 + mat2)) / obj.AnsatzNormDiag;
      end
    end % assembleAnsatzNormMatrix

    function M = assembleTestNormMatrix(obj, useNormalization)
      % Assemble the mass matrix of the discrete norm on the test space.
      %
      % Fast implementation that leverages a lot of simplification by hand.
      %
      % Parameters:
      %   useNormalization: toggle if the test functions should be normalized.
      %     you have to set both the class property useNormalization and this
      %     parameter to true to get a normalized matrix. @type logical @default
      %     true
      %
      % Return values:
      %   M: mass matrix @type sparsematrix

      % set default values
      if nargin == 1
        useNormalization = obj.useNormalization;
      end

      % Preparation for the sparse matrix
      Idx = ones(obj.nTestDim, 1);
      Idy = ones(obj.nTestDim, 1);
      Val = zeros(obj.nTestDim, 1);
      ctr = 1;

      % precompute the needed block matrices
      MtF = obj.assembleTemporalMassMatrix(obj.nTestTemporal, obj.nTestTemporal);
      AtF = obj.assembleTemporalStiffnessMatrix(obj.nTestTemporal);

      % spatial matrices
      MxF = obj.assembleSpatialMassMatrix(obj.nTestSpatial, obj.nTestSpatial);
      AxF = obj.assembleSpatialStiffnessMatrix(obj.nTestSpatial, obj.nTestSpatial);

      % first part; just like the ansatz norm
      mat1 = kron(MtF, MxF + AxF);

      % second part (initial condition)
      mat2 = obj.assembleSpatialMassMatrix(obj.nTestSpatialIC, obj.nTestSpatialIC);

      % Assemble the sparse mass matrix (with normalization, if needed)
      if ~useNormalization || ~obj.useNormalization
        M = sparse(obj.nTestDim, obj.nTestDim);
        M(1:(obj.nTestSpatial * obj.nTestTemporal), 1:(obj.nTestSpatial * obj.nTestTemporal)) = mat1;
        M((obj.nTestSpatial * obj.nTestTemporal + 1):end, (obj.nTestSpatial * obj.nTestTemporal + 1):end) = mat2;
      else
        M = sparse(obj.nTestDim, obj.nTestDim);
        M(1:(obj.nTestSpatial * obj.nTestTemporal), 1:(obj.nTestSpatial * obj.nTestTemporal)) = mat1;
        M((obj.nTestSpatial * obj.nTestTemporal + 1):end, (obj.nTestSpatial * obj.nTestTemporal + 1):end) = mat2;
        M = (obj.TestNormDiag \ M) / obj.TestNormDiag;
      end
    end

  end % methods

  methods % implementing AssemblyTensorAbstract

    function MtF = assembleTemporalMassMatrix(obj, nX, nY)
      % Assemble the temporal mass Matrix, that means we evaluate the integral
      % `\int_{I} \theta_k(t) \xi_m(t) \diff t` for `k = 1 \dots K` and
      % `m = 1 \dots K'`, where the temporal basis functions `\theta_k` and
      % `\xi_m` are defined as shifted Legendre polynomials.
      %
      % The parameters nX and nY correspond to `K` respectively `K'` and are
      % mainly here, because ansatz and test space can have a different number
      % of these basis functions.
      %
      % Parameters:
      %   nX: number of temporal basis functions `K` @type integer
      %   nY: number of temporal basis functions `K'` @type integer
      %
      % Return values:
      %   MtF: temporal mass matrix @type matrix

      tmp = (obj.tspan(2) - obj.tspan(1)) ./ (2 * ((1:min(nX, nY)) - 1) + 1);
      MtF = spdiags(tmp.', 0, nY, nX);
    end

    function CtF = assembleTemporalHalfStiffnessMatrix(obj, nX, nY)
      % Assemble the temporal "half stiffness" matrix, that means we evaluate the
      % integral `\int_{I} \theta'_k(t) \xi_m(t) \diff t` for `k = 1 \dots K` and
      % `m = 1 \dots K'`, where the temporal basis functions `\theta_k` and
      % `\xi_m` are defined as shifted Legendre polynomials.
      %
      % The parameters nX and nY correspond to `K` respectively `K'` and are
      % mainly here, because ansatz and test space can have a different number
      % of these basis functions.
      %
      % Parameters:
      %   nX: number of temporal basis functions `K` @type integer
      %   nY: number of temporal basis functions `K'` @type integer
      %
      % Return values:
      %   CtF: temporal "half stiffness" matrix @type matrix

      CtF = spdiags(2 * ones(nY, ceil(nX / 2)), 1:2:nX, nY, nX);
    end

    function AtF = assembleTemporalStiffnessMatrix(obj, nX)
      % Assemble the temporal stiffness matrix, that means we evaluate the
      % integral `\int_{I} \theta'_{k_1}(t) \theta'_{k_2}(t) \diff t` for
      % `k_1, k_2 = 1 \dots K`, where the temporal basis functions `\theta_k`
      % are defined as shifted Legendre polynomials.
      %
      % As this method is only useful for the assembly of the norm matrices,
      % there's only one parameter nX which corresponds to `K`.
      %
      % Parameters:
      %   nX: number of temporal basis functions `K` @type integer
      %
      % Return values:
      %   AtF: temporal stiffness matrix @type matrix
      %
      % @todo optimize!

      AtF = sparse(nX);
      for kdx1 = 1:nX
        % temporal intregal is zero if kdx1 + kdx2 is odd, so we only iterate
        % over the relevant indexes
        if mod(kdx1, 2) == 0
          startKdx2 = 2;
        else
          startKdx2 = 1;
        end
        for kdx2 = startKdx2:2:nX
          % evaluate the temporal integral;
          if kdx1 >= kdx2
            intTemporal = 2 * kdx2 * (kdx2 - 1) / ...
              (obj.tspan(2) - obj.tspan(1));
          else
            intTemporal = 2 * kdx1 * (kdx1 - 1) / ...
              (obj.tspan(2) - obj.tspan(1));
          end

          % save the evaluated integral
          AtF(kdx1, kdx2) = intTemporal;
        end
      end
    end

    function etF = assembleTemporalInitForwardVector(obj, nX)
      % Assemble the temporal row vector responsible for the propagation of the
      % initial condition in the case of the forward propagator, that means we
      % evaluate `\theta_k(0)` for `k = 1 \dots K` with Legendre polynomials
      % `\theta_k`.
      %
      % Parameters:
      %   nX: corresponds to `K` @type integer
      %
      % Return values:
      %   eFtF: forward propagation vector @type vector

      etF = (-1).^((1:nX) - 1);
    end

    function etF = assembleTemporalInitBackwardVector(obj, nX)
      % Assemble the temporal row vector responsible for the propagation of the
      % initial condition in the case of the forward propagator, that means we
      % evaluate `\theta_k(L)` for `k = 1 \dots K` with Legendre polynomials
      % `\theta_k`.
      %
      % Parameters:
      %   nX: corresponds to `K` @type integer
      %
      % Return values:
      %   eBtF: backward propagation vector @type vector

      etF = ones(1, nX);
    end

    function MxF = assembleSpatialMassMatrix(obj, nX, nY)
      % Assemble the spatial mass matrix, that means we evaluate the integral
      % `\int_{\Omega} \sigma_j(x) \sigma_l(x) \diff x` for some `j` and `l`,
      % which are Fourier basis functions in this case.
      %
      % Parameters:
      %   nX: the limit of iteration for the `j` index @type integer
      %   nY: the limit of iteration for the `l` index @type integer
      %
      % Return values:
      %   MxF: spatial mass matrix @type matrix

      tmp    = (obj.xspan(2) - obj.xspan(1)) * ones(min(nX, nY), 1) / 2;
      tmp(1) = tmp(1) * 2;
      MxF    = spdiags(tmp, 0, nY, nX);
    end

    function AxF = assembleSpatialStiffnessMatrix(obj, nX, nY)
      % Assemble the spatial stiffness matrix, that means we evaluate the
      % integral `\int_{\Omega} \sigma'_j(x) \sigma'_l(x) \diff x` for some `j`
      % and `l`.
      %
      % Parameters:
      %   nX: the limit of iteration for the `j` index @type integer
      %   nY: the limit of iteration for the `l` index @type integer
      %
      % Return values:
      %   AxF: spatial stiffness matrix @type matrix
      %
      % @todo optimize

      tmp = zeros(min(nX, nY), 1);
      for jdx = 2:min(nX, nY)
        if mod(jdx, 2) == 0
          tmp(jdx) = (pi * jdx)^2 / ( 2 * (obj.xspan(2) - obj.xspan(1)));
        else
          tmp(jdx) = (pi * (jdx - 1))^2 / ( 2 * (obj.xspan(2) - obj.xspan(1)));
        end
      end
      AxF = spdiags(tmp, 0, nY, nX);
    end

  end % implementing AssemblyTensorAbstract

  methods % implementing AssemblyGlobalAbstract

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

    function val = spatialBasisFuncDerivative(obj, index, x)
      % First derivatives of spatial basis functions.
      %
      % Evaluates the first derivative of spatial basis function, Fourier
      % functions, with the given index for the given values of x. Can be used
      % to define function handles and for numerical integration.
      %
      % In this case the first derivatives of the spatial basis functions are
      % of the type
      % ``\begin{cases}
      %   0, & i = 1\\
      %   (\pi i / L) \cos(\pi i x / L), & i~\text{even}\\
      %   - (\pi i / L) \sin(\pi i x / L), & i~\text{odd and}~i > 1
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
        val = zeros(size(x, 1), size(x, 2));
      elseif mod(index, 2) == 0
        val = (pi * index / obj.xspan(2)) * cos(pi * index * x / obj.xspan(2));
      else
        val = - (pi * (index - 1) / obj.xspan(2)) * ...
          sin(pi * (index - 1) * x / obj.xspan(2));
      end
    end % spatialBasisFuncDerivative

    function val = temporalBasisFunc(obj, index, t, tspan)
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

      if nargin == 3
        tspan = obj.tspan;
      end

      val = legendrePolynomial(t, index - 1, tspan);
    end % spatialBasisFuncDerivative

    function val = temporalBasisFuncDerivative(obj, index, t, tspan)
      % First derivative of temporal basis functions.
      %
      % Evaluates the first derivative of a temporal basis function, shifted
      % Legendre polynomials, with the given index for the given values of t.
      % Can be used to define function handles and for numerical integration.
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

      if nargin == 3
        tspan = obj.tspan;
      end

      val = legendrePolynomialDerivative(t, index - 1, tspan);
    end % temporalBasisFuncDerivative

    function val = normOfAnsatzFunc(obj, jdx, kdx)
      % Calculate the ansatz space norm of a ansatz basis function.
      %
      % Parameter:
      %   jdx: index of spatial basis function @type integer
      %   kdx: index of temporal basis function @type integer
      %
      % Return values:
      %   val: norm of product of spatial and temporal basis function @type double

      % evaluate the spatial integrals
      if jdx == 1
        % constant
        intSpatialF = obj.xspan(2);
        intSpatialS = obj.xspan(2);
      elseif mod(jdx, 2) == 1
        % cosine
        intSpatialF = obj.xspan(2) / 2 + (pi * (jdx - 1))^2 / (2 * obj.xspan(2));
        intSpatialS = obj.xspan(2) / 2;
      else
        % sine
        intSpatialF = obj.xspan(2) / 2 + (pi * jdx)^2 / (2 * obj.xspan(2));
        intSpatialS = obj.xspan(2) / 2;
      end

      % evaluate the temporal integrals;
      intTemporalF = (obj.tspan(2) - obj.tspan(1)) / (2 * (kdx - 1) + 1);
      intTemporalS = 2 * kdx * (kdx - 1) / (obj.tspan(2) - obj.tspan(1));

      val = sqrt(intTemporalF * intSpatialF + intTemporalS * intSpatialS);
    end % normOfAnsatzFunc

    function val = normOfTestFunc(obj, ldx, mdx, ndx)
      % Calculate the test space norm of a test basis function.
      %
      % Parameter:
      %   jdx: index of spatial basis function @type integer
      %   kdx: index of temporal basis function @type integer
      %
      % Return values:
      %   val: norm of product of spatial and temporal basis function @type double

      tmp1 = 0;
      if ldx > 0 && mdx > 0
        % evaluate the spatial integrals
        if ldx == 1
          % constant
          intSpatial = obj.xspan(2);
        elseif mod(ldx, 2) == 1
          % cosine
          intSpatial = obj.xspan(2) / 2 + (pi * (ldx - 1))^2 / (2 * obj.xspan(2));
        else
          % sine
          intSpatial = obj.xspan(2) / 2 + (pi * ldx)^2 / (2 * obj.xspan(2));
        end

        % evaluate the temporal integral
        intTemporal = (obj.tspan(2) - obj.tspan(1)) / (2 * (mdx - 1) + 1);
        tmp1 = intSpatial * intTemporal;
      end

      tmp2 = 0;
      if ndx > 0
        if ndx == 1
          tmp2 = obj.xspan(2);
        else
          tmp2 = obj.xspan(2) / 2;
        end
      end

      val = sqrt(tmp1 + tmp2);
    end % normOfTestFunc

    function precomputeNormalization(obj)
      % Computes the normalization coefficients if normalization is used.
      %
      % This is done by assembly of the discrete norm matrices for the ansatz
      % and test subspaces. The diagonals of these matrices are the squares of
      % the wanted normalization coefficients.

      if obj.useNormalization
        % we want to use normalization
        obj.AnsatzNormDiag = sqrt(spdiags(...
          spdiags(obj.assembleAnsatzNormMatrix(false), 0), 0, ...
          obj.nAnsatzDim, obj.nAnsatzDim));
        obj.TestNormDiag   = sqrt(spdiags(...
          spdiags(obj.assembleTestNormMatrix(false), 0), 0, ...
          obj.nTestDim, obj.nTestDim));
      else
        % we don't want to use normalization
        obj.AnsatzNormDiag = speye(obj.nAnsatzDim);
        obj.TestNormDiag   = speye(obj.nTestDim);
      end
    end % precomputeNormalization

  end % implementing AssemblyGlobalAbstract

end % classdef
