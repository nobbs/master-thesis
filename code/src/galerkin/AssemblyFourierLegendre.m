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
    %% Normalization and norm related stuff

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


    %% Fast implementation of the assembly methods

    function [M1, M2, M3, M4F, M4B] = assembleFieldIndependentMatrix(obj)
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

      % Preparation for the sparse matrix
      Idx = ones(obj.nAnsatzDim, 1);
      Idy = ones(obj.nAnsatzDim, 1);
      Val = zeros(obj.nAnsatzDim, 1);
      ctr = 1;

      % optimization...
      xwidth = obj.xspan(2);

      % Calculate the first term `\int_{I} \skp{u_t(t)}{v_1(t)}{L_2(\Omega)}
      % \diff t`.
      for jdx = 1:min(obj.nAnsatzSpatial, obj.nTestSpatial)
        for kdx = 1:obj.nAnsatzTemporal
          for mdx = (kdx - 1):-2:1
            % evaluate spatial integral
            if jdx == 1
              intSpatial = xwidth;
            else
              intSpatial = xwidth / 2;
            end

            % evaluate temporal integral
            intTemporal = 2;

            % save the evaluated integrals
            Idx(ctr) = (jdx - 1) * obj.nAnsatzTemporal + kdx;
            Idy(ctr) = (jdx - 1) * obj.nTestTemporal + mdx;
            Val(ctr) = intSpatial * intTemporal;
            ctr = ctr + 1;
          end
        end
      end

      % Assemble the stiffness matrix for this part
      M1 = sparse(Idy, Idx, Val, obj.nTestDim, obj.nAnsatzDim);

      % Reset for the next part
      Idx = ones(obj.nAnsatzDim, 1);
      Idy = ones(obj.nAnsatzDim, 1);
      Val = zeros(obj.nAnsatzDim, 1);
      ctr = 1;

      % Calculate the second term `\int_{I} \skp{\grad u(t)}{\grad
      % v_1(t)}{L_2(\Omega)} \diff t`.
      for jdx = 1:min(obj.nAnsatzSpatial, obj.nTestSpatial)
        for kdx = 1:min(obj.nAnsatzTemporal, obj.nTestTemporal)
          % evaluate spatial integral
          if jdx == 1
            intSpatial = 0;
          elseif mod(jdx, 2) == 0
            intSpatial = (pi * jdx)^2 / ( 2 * xwidth);
          else
            intSpatial = (pi * (jdx - 1))^2 / ( 2 * xwidth);
          end

          % evaluate temporal integral
          intTemporal = (obj.tspan(2) - obj.tspan(1)) / (2 * (kdx - 1) + 1);

          % save the evaluated integrals
          Idx(ctr) = (jdx - 1) * obj.nAnsatzTemporal + kdx;
          Idy(ctr) = (jdx - 1) * obj.nTestTemporal + kdx;
          Val(ctr) = intTemporal * intSpatial;
          ctr = ctr + 1;
        end
      end

      % Assemble the stiffness matrix for this part
      M2 = sparse(Idy, Idx, Val, obj.nTestDim, obj.nAnsatzDim);

      % Reset for the next part
      Idx = ones(obj.nAnsatzDim, 1);
      Idy = ones(obj.nAnsatzDim, 1);
      Val = zeros(obj.nAnsatzDim, 1);
      ctr = 1;

      % Calculate the third term `\int_{I} \mu \skp{u(t)}{v_1(t)}{L_2(\Omega)}
      % \diff t`.
      for jdx = 1:min(obj.nAnsatzSpatial, obj.nTestSpatial)
        for kdx = 1:min(obj.nAnsatzTemporal, obj.nTestTemporal)
          % evaluate spatial integral
          if jdx == 1
            intSpatial = xwidth;
          else
            intSpatial = xwidth / 2;
          end

          % evaluate temporal integral
          intTemporal = (obj.tspan(2) - obj.tspan(1)) / (2 * (kdx - 1) + 1);

          % save the evaluated integrals
          Idx(ctr) = (jdx - 1) * obj.nAnsatzTemporal + kdx;
          Idy(ctr) = (jdx - 1) * obj.nTestTemporal + kdx;
          Val(ctr) = intTemporal * intSpatial;
          ctr = ctr + 1;
        end
      end

      % Assemble the stiffness matrix for this part
      M3 = sparse(Idy, Idx, Val, obj.nTestDim, obj.nAnsatzDim);

      % Reset for the next part
      Idx  = ones(obj.nAnsatzDim, 1);
      Idy  = ones(obj.nAnsatzDim, 1);
      ValF = zeros(obj.nAnsatzDim, 1);
      ValB = zeros(obj.nAnsatzDim, 1);
      ctr  = 1;

      % Calculate the fourth term `\skp{u(0)}{v_2}{L_2(\Omega)}`.
      for jdx = 1:min(obj.nAnsatzSpatial, obj.nTestSpatialIC)
        for kdx = 1:obj.nAnsatzTemporal
          % evaluate spatial integral
          if jdx == 1
            intSpatial = xwidth;
          else
            intSpatial = xwidth / 2;
          end

          % evaluate temporal integral
          intTemporalF  = (-1)^(kdx - 1);
          intTemporalB = 1;

          % save the evaluated integrals
          Idx(ctr)  = (jdx - 1) * obj.nAnsatzTemporal + kdx;
          Idy(ctr)  = obj.nTestSpatial * obj.nTestTemporal + jdx;
          ValF(ctr) = intTemporalF * intSpatial;
          ValB(ctr) = intTemporalB * intSpatial;
          ctr = ctr + 1;
        end
      end

      % Assemble the stiffness matrix for this part
      M4F = sparse(Idy, Idx, ValF, obj.nTestDim, obj.nAnsatzDim);
      M4B = sparse(Idy, Idx, ValB, obj.nTestDim, obj.nAnsatzDim);

      % normalize if needed
      if obj.useNormalization
        M1  = (obj.TestNormDiag \ M1) / obj.AnsatzNormDiag;
        M2  = (obj.TestNormDiag \ M2) / obj.AnsatzNormDiag;
        M3  = (obj.TestNormDiag \ M3) / obj.AnsatzNormDiag;
        M4F = (obj.TestNormDiag \ M4F) / obj.AnsatzNormDiag;
        M4B = (obj.TestNormDiag \ M4B) / obj.AnsatzNormDiag;
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

            % iterate over the index of the temporal basis functions of the
            % ansatz and test subspaces. since both use Legendre polynomials,
            % we can simplify the occurring integrals to the following
            % expression
            for kdx = 1:min(obj.nAnsatzTemporal, obj.nTestTemporal)
              % evaluate temporal integral
              intTemporal = (obj.tspan(2) - obj.tspan(1)) / (2 * (kdx - 1) + 1);

              % save the evaluated integrals
              Idx(ctr) = (jdx - 1) * obj.nAnsatzTemporal + kdx;
              Idy(ctr) = (ldx - 1) * obj.nTestTemporal + kdx;
              Val(ctr) = intTemporal * intSpatial;
              ctr = ctr + 1;
            end
          end
        end

        % create the sparse matrix
        O{cdx} = sparse(Idy, Idx, Val, obj.nTestDim, obj.nAnsatzDim);

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

            % iterate over the index of the temporal basis functions of the
            % ansatz and test subspaces. since both use Legendre polynomials,
            % we can simplify the occurring integrals to the following
            % expression
            for kdx = 1:min(obj.nAnsatzTemporal, obj.nTestTemporal)
              % evaluate temporal integral
              intTemporal = (obj.tspan(2) - obj.tspan(1)) / (2 * (kdx - 1) + 1);

              % save the evaluated integrals
              Idx(ctr) = (jdx - 1) * obj.nAnsatzTemporal + kdx;
              Idy(ctr) = (ldx - 1) * obj.nTestTemporal + kdx;
              Val(ctr) = intTemporal * intSpatial;
              ctr = ctr + 1;
            end
          end
        end

        % create the sparse matrix
        O{cdx} = sparse(Idy, Idx, Val, obj.nTestDim, obj.nAnsatzDim);

        % normalize if needed
        if obj.useNormalization
          O{cdx} = (obj.TestNormDiag \ O{cdx}) / obj.AnsatzNormDiag;
        end
      end
    end % assembleFieldDependentMatrixForFourierSeries


    %% Assembly of the load vector for different kinds of data

    function F = assembleVectorFromSpatialCoeffs(obj, coeffs)
      % Assemble the load vector for no source and initial data given by
      % coeffients corresponding to the spatial basis functions.
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
      % coeffients corresponding to the spatial basis and temporal basis
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
          sumcoeffs(idx) = sum(solCoeffs(((idx - 1) * ...
            obj.nAnsatzTemporal) + Kdx));
        end
      else
        % sum with alternating sign, because the temporal basis function with
        % index i has the value (-1)^i at the start point.
        for idx = 1:obj.nAnsatzSpatial
          sumcoeffs(idx) = sum((-1).^(Kdx - 1) * ...
            solCoeffs(((idx - 1) * obj.nAnsatzTemporal) + Kdx));
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

      % Preparation for the sparse matrix
      Idx = ones(obj.nAnsatzDim, 1);
      Idy = ones(obj.nAnsatzDim, 1);
      Val = zeros(obj.nAnsatzDim, 1);
      ctr = 1;

      % first part `\norm{u}_{L_2(I; V)}`: iterate over spatial and temporal
      % basis functions; respectively only in one dimension because boths sets
      % of basis functions are orthogonal in itself.
      for jdx = 1:obj.nAnsatzSpatial
        for kdx = 1:obj.nAnsatzTemporal
          % handle the different cases of the spatial basis function and its
          % first derivative
          if jdx == 1
            % constant function
            intSpatial = obj.xspan(2);
          elseif mod(jdx, 2) == 1
            % cosine
            intSpatial = obj.xspan(2) / 2 + (pi * (jdx - 1))^2 / ...
              (2 * obj.xspan(2));
          else
            % sine
            intSpatial = obj.xspan(2) / 2 + (pi * jdx)^2 / (2 * obj.xspan(2));
          end

          % evaluate temporal integral
          intTemporal = (obj.tspan(2) - obj.tspan(1)) / (2 * (kdx - 1) + 1);

          % save the evaluated integrals
          Idx(ctr) = (jdx - 1) * obj.nAnsatzTemporal + kdx;
          Idy(ctr) = (jdx - 1) * obj.nAnsatzTemporal + kdx;
          Val(ctr) = intTemporal * intSpatial;
          ctr = ctr + 1;
        end
      end

      % second part `norm{u_t}_{L_2(I; V')}`: again iterate over spatial and
      % temporal basis functions; this time with two indexes for the temporal
      % component, as the first derivatives of the temporal basis functions are
      % no longer orthogonal
      for jdx = 1:obj.nAnsatzSpatial
        for kdx1 = 1:obj.nAnsatzTemporal
          % temporal intregal is zero if kdx1 + kdx2 is odd, so we only iterate
          % over the relevant indexes
          if mod(kdx1, 2) == 0
            startKdx2 = 2;
          else
            startKdx2 = 1;
          end
          for kdx2 = startKdx2:2:obj.nAnsatzTemporal
            % evaluate the spatial integral; handle the different cases of the
            % spatial basis function
            if jdx == 1
              % constant function
              intSpatial = obj.xspan(2);
            else
              % sine and cosine
              intSpatial = obj.xspan(2) / 2;
            end

            % evaluate the temporal integral;
            if kdx1 >= kdx2
              intTemporal = 2 * kdx2 * (kdx2 - 1) / ...
                (obj.tspan(2) - obj.tspan(1));
            else
              intTemporal = 2 * kdx1 * (kdx1 - 1) / ...
                (obj.tspan(2) - obj.tspan(1));
            end

            % save the evaluated integrals
            Idx(ctr) = (jdx - 1) * obj.nAnsatzTemporal + kdx1;
            Idy(ctr) = (jdx - 1) * obj.nAnsatzTemporal + kdx2;
            Val(ctr) = intTemporal * intSpatial;
            ctr = ctr + 1;
          end
        end
      end

      % Assemble the sparse mass matrix (with normalization, if needed)
      if ~useNormalization || ~obj.useNormalization
        M = sparse(Idy, Idx, Val, obj.nAnsatzDim, obj.nAnsatzDim);
      else
        M = (obj.AnsatzNormDiag \ ...
          sparse(Idy, Idx, Val, obj.nAnsatzDim, obj.nAnsatzDim)) / ...
          obj.AnsatzNormDiag;
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

      % first part; iterate over the spatial basis for the first component
      for ldx = 1:obj.nTestSpatial
        % iterate over temporal basis for the first component
        for mdx = 1:obj.nTestTemporal
          % handle the different cases of the spatial basis function
          if ldx == 1
            % constant
            intSpatial = obj.xspan(2);
          elseif mod(ldx, 2) == 1
            % cosine
            intSpatial = obj.xspan(2) / 2 + (pi * (ldx - 1))^2 / ...
              (2 * obj.xspan(2));
          else
            % sine
            intSpatial = obj.xspan(2) / 2 + (pi * ldx)^2 / (2 * obj.xspan(2));
          end

          % evaluate the temporal integral
          intTemporal = (obj.tspan(2) - obj.tspan(1)) / (2 * (mdx - 1) + 1);

          % save the evaluated integrals
          Idx(ctr) = (ldx - 1) * obj.nTestTemporal + mdx;
          Idy(ctr) = (ldx - 1) * obj.nTestTemporal + mdx;
          Val(ctr) = intTemporal * intSpatial;
          ctr = ctr + 1;
        end
      end

      % second part (initial condition); iterate over the spatial basis for the
      % second component
      for ndx = 1:obj.nTestSpatialIC
        % handle the different cases of the spatial basis function
        if ndx == 1
          % constant function
          val = obj.xspan(2);
        else
          % cosine and sine
          val = obj.xspan(2) / 2;
        end

        % save the evaluated integrals
        Idx(ctr) = obj.nTestSpatial * obj.nTestTemporal + ndx;
        Idy(ctr) = obj.nTestSpatial * obj.nTestTemporal + ndx;
        Val(ctr) = val;
        ctr = ctr + 1;
      end

      % Assemble the sparse mass matrix (with normalization, if needed)
      if ~useNormalization || ~obj.useNormalization
        M = sparse(Idy, Idx, Val, obj.nTestDim, obj.nTestDim);
      else
        M = (obj.TestNormDiag \ ...
          sparse(Idy, Idx, Val, obj.nTestDim, obj.nTestDim)) / obj.TestNormDiag;
      end
    end

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


    %% spatial and temporal basis functions as convenient function handles

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
    end

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

  end % methods

end % classdef
