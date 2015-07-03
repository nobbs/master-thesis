classdef AssemblyFourierLegendreSlow < AssemblyFourierLegendre
% Assembly of the stiffness matrix and load vector for a Galerkin method using
% global basis functions with periodic boundary conditions.
%
% This class only contains slow implementations of some crucial methods of
% AssemblyFourierLegendre using numerical integration as the main ingredient to
% assemble the needed matrices and vectors. This should only be used as a
% reference for comparisons and unit tests!
%
% See also:
%   AssemblyFourierLegendre

  methods
    function [M1, M2, M3, M4] = assembleFieldIndependentMatrixSlow(obj)
      % Assemble the field-independent part of the stiffness matrix. Slow
      % implementation, only for comparison!
      %
      % This is done by numerical integration of all the occurring products of
      % spatial and temporal basis functions, so it's kinda slow and not that
      % accurate. Only really useful to test the faster, optimized version.
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

      % Calculate the first term `\int_{I} \skp{u_t(t)}{v_1(t)}{L_2(\Omega)}
      % \diff t`.
      for jdx = 1:obj.nAnsatzSpatial
        for kdx = 1:obj.nAnsatzTemporal
          for ldx = 1:obj.nTestSpatial
            for mdx = 1:obj.nTestTemporal
              val = integral(@(x) obj.spatialBasisFunc(jdx, x) .* ...
                obj.spatialBasisFunc(ldx, x), obj.xspan(1), obj.xspan(2)) * ...
                integral(@(t) obj.temporalBasisFuncDerivative(kdx, t) .* ...
                obj.temporalBasisFunc(mdx, t), obj.tspan(1), obj.tspan(2));

              if obj.useNormalization
                val = val / (obj.normOfAnsatzFuncSlow(jdx, kdx) * ...
                  obj.normOfTestFuncSlow(ldx, mdx, 0));
              end

              % only consider values above a given threshold (since a lot of
              % zero valued integrals won't get exact zero)
              if abs(val) > sqrt(eps)
                Idx(ctr) = (kdx - 1) * obj.nAnsatzSpatial + jdx;
                Idy(ctr) = (mdx - 1) * obj.nTestSpatial + ldx;
                Val(ctr) = val;
                ctr = ctr + 1;
              end
            end
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
      for jdx = 1:obj.nAnsatzSpatial
        for kdx = 1:obj.nAnsatzTemporal
          for ldx = 1:obj.nTestSpatial
            for mdx = 1:obj.nTestTemporal
              intSpatial = integral(@(x) ...
                obj.spatialBasisFuncDerivative(jdx, x) .* ...
                obj.spatialBasisFuncDerivative(ldx, x), ...
                obj.xspan(1), obj.xspan(2));
              intTemporal = integral(@(t) obj.temporalBasisFunc(kdx, t) .* ...
                obj.temporalBasisFunc(mdx, t), obj.tspan(1), obj.tspan(2));

              val = intSpatial * intTemporal;

              if obj.useNormalization
                val = val / (obj.normOfAnsatzFuncSlow(jdx, kdx) * ...
                  obj.normOfTestFuncSlow(ldx, mdx, 0));
              end

              % only consider values above a given threshold (since a lot of
              % zero valued integrals won't get exact zero)
              if abs(val) > sqrt(eps)
                Idx(ctr) = (kdx - 1) * obj.nAnsatzSpatial + jdx;
                Idy(ctr) = (mdx - 1) * obj.nTestSpatial + ldx;
                Val(ctr) = val;
                ctr = ctr + 1;
              end
            end
          end
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
      for jdx = 1:obj.nAnsatzSpatial
        for kdx = 1:obj.nAnsatzTemporal
          for ldx = 1:obj.nTestSpatial
            for mdx = 1:obj.nTestTemporal
              val = integral(@(x) obj.spatialBasisFunc(jdx, x) .* ...
                obj.spatialBasisFunc(ldx, x), obj.xspan(1), obj.xspan(2)) * ...
                integral(@(t) obj.temporalBasisFunc(kdx, t) .* ...
                obj.temporalBasisFunc(mdx, t), obj.tspan(1), obj.tspan(2));

              if obj.useNormalization
                val = val / (obj.normOfAnsatzFuncSlow(jdx, kdx) * ...
                  obj.normOfTestFuncSlow(ldx, mdx, 0));
              end

              % only consider values above a given threshold (since a lot of
              % zero valued integrals won't get exact zero)
              if abs(val) > sqrt(eps)
                Idx(ctr) = (kdx - 1) * obj.nAnsatzSpatial + jdx;
                Idy(ctr) = (mdx - 1) * obj.nTestSpatial + ldx;
                Val(ctr) = val;
                ctr = ctr + 1;
              end
            end
          end
        end
      end

      % Assemble the stiffness matrix for this part
      M3 = sparse(Idy, Idx, Val, obj.nTestDim, obj.nAnsatzDim);

      % Reset for the next part
      Idx = ones(obj.nAnsatzDim, 1);
      Idy = ones(obj.nAnsatzDim, 1);
      Val = zeros(obj.nAnsatzDim, 1);
      ctr = 1;

      % Calculate the fourth term `\skp{u(0)}{v_2}{L_2(\Omega)}`.
      for jdx = 1:obj.nAnsatzSpatial
        for kdx = 1:obj.nAnsatzTemporal
          for ndx = 1:obj.nTestSpatialIC
            val = obj.temporalBasisFunc(kdx, obj.tspan(1)) * ...
              integral(@(x) obj.spatialBasisFunc(jdx, x) .* ...
              obj.spatialBasisFunc(ndx, x), obj.xspan(1), obj.xspan(2));

            if obj.useNormalization
              val = val / (obj.normOfAnsatzFuncSlow(jdx, kdx) * ...
                obj.normOfTestFuncSlow(0, 0, ndx));
            end

            % only consider values above a given threshold (since a lot of zero
            % valued integrals won't get exact zero)
            if abs(val) > sqrt(eps)
              Idx(ctr) = (kdx - 1) * obj.nAnsatzSpatial + jdx;
              Idy(ctr) = obj.nTestTemporal * obj.nTestSpatial + ndx;
              Val(ctr) = val;
              ctr = ctr + 1;
            end
          end
        end
      end

      % Assemble the stiffness matrix for this part
      M4 = sparse(Idy, Idx, Val, obj.nTestDim, obj.nAnsatzDim);
    end % assembleFieldIndependentMatrixSlow

    function O = assembleFieldDependentMatrixForFourierSeriesSlow(obj, nCoeff)
      % Assemble only the field-dependent components of the stiffness matrix
      % based on a fourier series expansion of the field. Slow implementation,
      % only for comparison!
      %
      % This is done by numerical integration of all the occurring products of
      % spatial, temporal and field basis functions, so it's kinda slow and not
      % that accurate. Only really useful to test the faster, optimized version.
      %
      % Parameters:
      %   nCoeff: number of fourier basis functions @type integer
      %
      % Return values:
      %   O: cellarray of the partial stiffness matrices @type struct

      % create the cellarray that will hold the matrices
      O = cell(nCoeff);

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
            if mod(cdx, 2) == 0
              % cdx even: sine function sin(pi * cdx * x / L)
              intSpatial = integral(@(x) sin(pi * cdx * x / obj.xspan(2)) .* ...
                obj.spatialBasisFunc(jdx, x) .* obj.spatialBasisFunc(ldx, x), ...
                obj.xspan(1), obj.xspan(2));
            else
              % cdx odd: cosine function cos(pi * (cdx + 1) * x)
              intSpatial = integral(@(x) ...
                cos(pi * (cdx + 1) * x / obj.xspan(2)) .* ...
                obj.spatialBasisFunc(jdx, x) .* obj.spatialBasisFunc(ldx, x), ...
                obj.xspan(1), obj.xspan(2));
            end

            % only consider values above a given threshold (since a lot of zero
            % valued integrals won't get exact zero)
            if abs(intSpatial) > sqrt(eps)
              % iterate over the index of the temporal basis functions of the
              % ansatz and test subspaces. since both use Legendre polynomials,
              % we can simplify the occurring integrals to the following
              % expression
              for kdx = 1:obj.nAnsatzTemporal
                for mdx = 1:obj.nTestTemporal
                  % evaluate temporal integral
                  intTemporal = integral(@(t) ...
                    obj.temporalBasisFunc(kdx, t) .* ...
                    obj.temporalBasisFunc(mdx, t), obj.tspan(1), obj.tspan(2));

                  val = intSpatial * intTemporal;

                  if obj.useNormalization
                    val = val / (obj.normOfAnsatzFuncSlow(jdx, kdx) * ...
                      obj.normOfTestFuncSlow(ldx, mdx, 0));
                  end

                  if abs(val) > sqrt(eps)
                    % save the evaluated integrals
                    Idx(ctr) = (kdx - 1) * obj.nAnsatzSpatial + jdx;
                    Idy(ctr) = (mdx - 1) * obj.nTestSpatial + ldx;
                    Val(ctr) = val;
                    ctr = ctr + 1;
                  end
                end
              end
            end
          end
        end

        % create the sparse matrix
        O{cdx} = sparse(Idy, Idx, Val, obj.nTestDim, obj.nAnsatzDim);
      end
    end % assembleFieldDependentMatrixForFourierSeriesSlow

    function O = assembleFieldDependentMatrixForSineSeriesSlow(obj, nCoeff)
      % Assemble only the field-dependent components of the stiffness matrix
      % based on a sine series expansion of the field. Slow implementation, only
      % for comparison!
      %
      % This is done by numerical integration of all the occurring products of
      % spatial, temporal and field basis functions, so it's kinda slow and not
      % that accurate. Only really useful to test the faster, optimized version.
      %
      % Parameters:
      %   nCoeff: number of sine basis functions @type integer
      %
      % Return values:
      %   O: cellarray of the partial stiffness matrices @type struct

      % create the cellarray that will hold the matrices
      O = cell(nCoeff);

      % iterate over the index of the sine basis function (of the sine series
      % expansion of the field)
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
            % sine / cosine product we have to consider, but this is all handled
            % by spatialBasisFunc
            intSpatial = integral(@(x) sin(pi * cdx * x / obj.xspan(2)) .* ...
              obj.spatialBasisFunc(jdx, x) .* obj.spatialBasisFunc(ldx, x), ...
              obj.xspan(1), obj.xspan(2));

            % only consider values above a given threshold (since a lot of zero
            % valued integrals won't get exact zero)
            if abs(intSpatial) > sqrt(eps)
              % iterate over the index of the temporal basis functions of the
              % ansatz and test subspaces. since both use Legendre polynomials,
              % we can simplify the occurring integrals to the following
              % expression
              for kdx = 1:obj.nAnsatzTemporal
                for mdx = 1:obj.nTestTemporal
                  % evaluate temporal integral
                  intTemporal = integral(@(t) obj.temporalBasisFunc(kdx, t) .* ...
                    obj.temporalBasisFunc(mdx, t), obj.tspan(1), obj.tspan(2));

                  val = intSpatial * intTemporal;
                  if obj.useNormalization
                    val = val / (obj.normOfAnsatzFuncSlow(jdx, kdx) * ...
                      obj.normOfTestFuncSlow(ldx, mdx, 0));
                  end

                  if abs(val) > sqrt(eps)
                    % save the evaluated integrals
                    Idx(ctr) = (kdx - 1) * obj.nAnsatzSpatial + jdx;
                    Idy(ctr) = (mdx - 1) * obj.nTestSpatial + ldx;
                    Val(ctr) = val;
                    ctr = ctr + 1;
                  end
                end
              end
            end
          end
        end

        % create the sparse matrix
        O{cdx} = sparse(Idy, Idx, Val, obj.nTestDim, obj.nAnsatzDim);
      end
    end % assembleFieldDependentMatrixForSineSeriesSlow

    function val = normOfTestFuncSlow(obj, ldx, mdx, ndx)
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
        tmp1 = integral(@(t) obj.temporalBasisFunc(mdx, t).^2, obj.tspan(1), obj.tspan(2)) * ...
          integral(@(x) obj.spatialBasisFunc(ldx, x).^2 + obj.spatialBasisFuncDerivative(ldx, x).^2, obj.xspan(1), obj.xspan(2));
      end

      tmp2 = 0;
      if ndx > 0
        tmp2 = integral(@(x) obj.spatialBasisFunc(ndx, x).^2, obj.xspan(1), obj.xspan(2));
      end

      val = sqrt(tmp1 + tmp2);
    end % normOfTestFuncSlow

    function val = normOfAnsatzFuncSlow(obj, jdx, kdx)
      % Calculate the ansatz space norm of a ansatz basis function.
      %
      % Parameter:
      %   jdx: index of spatial basis function @type integer
      %   kdx: index of temporal basis function @type integer
      %
      % Return values:
      %   val: norm of product of spatial and temporal basis function @type double

      tmp1 = integral(@(t) obj.temporalBasisFunc(kdx, t).^2, obj.tspan(1), obj.tspan(2)) * ...
        integral(@(x) obj.spatialBasisFunc(jdx, x).^2 + obj.spatialBasisFuncDerivative(jdx, x).^2, obj.xspan(1), obj.xspan(2));

      tmp2 = integral(@(t) obj.temporalBasisFuncDerivative(kdx, t).^2, obj.tspan(1), obj.tspan(2)) * ...
        integral(@(x) obj.spatialBasisFunc(jdx, x).^2, obj.xspan(1), obj.xspan(2));

      val = sqrt(tmp1 + tmp2);
    end % normOfAnsatzFuncSlow

    function M = assembleTestNormMatrixSlow(obj)
      % Assemble the mass matrix of the discrete norm on the test space.
      %
      % Slow implementation using numerical quadrature. Only useful for
      % comparison and testing of faster methods!
      %
      % Return values:
      %   M: mass matrix @type sparsematrix

      % Preparation for the sparse matrix
      Idx = ones(obj.nTestDim, 1);
      Idy = ones(obj.nTestDim, 1);
      Val = zeros(obj.nTestDim, 1);
      ctr = 1;

      % iterate over the spatial basis for the first component
      for ldx1 = 1:obj.nTestSpatial
        for ldx2 = 1:obj.nTestSpatial
          % iterate over temporal basis for the first component
          for mdx1 = 1:obj.nTestTemporal
            for mdx2 = 1:obj.nTestTemporal
              intSpatial  = integral(@(x) obj.spatialBasisFunc(ldx1, x) .* obj.spatialBasisFunc(ldx2, x) + obj.spatialBasisFuncDerivative(ldx1, x) .* obj.spatialBasisFuncDerivative(ldx2, x), obj.xspan(1), obj.xspan(2));
              intTemporal = integral(@(t) obj.temporalBasisFunc(mdx1, t) .* obj.temporalBasisFunc(mdx2, t), obj.tspan(1), obj.tspan(2));
              val = intTemporal * intSpatial;

              if obj.useNormalization
                val = val / (obj.normOfTestFuncSlow(ldx1, mdx1, 0) * ...
                  obj.normOfTestFuncSlow(ldx2, mdx2, 0));
              end

              % only consider values above a given threshold (since a lot of zero
              % valued integrals won't get exact zero)
              if abs(val) > sqrt(eps)
                Idx(ctr) = (mdx1 - 1) * obj.nTestSpatial + ldx1;
                Idy(ctr) = (mdx2 - 1) * obj.nTestSpatial + ldx2;
                Val(ctr) = val;
                ctr = ctr + 1;
              end
            end
          end
        end
      end

      % iterate over the spatial basis for the second component
      for ndx1 = 1:obj.nTestSpatialIC
        for ndx2 = 1:obj.nTestSpatialIC
          val = integral(@(x) obj.spatialBasisFunc(ndx1, x) .* obj.spatialBasisFunc(ndx2, x), obj.xspan(1), obj.xspan(2));

          if obj.useNormalization
            val = val / (obj.normOfTestFuncSlow(0, 0, ndx1) * obj.normOfTestFuncSlow(0, 0, ndx2));
          end

          % only consider values above a given threshold (since a lot of zero
          % valued integrals won't get exact zero)
          if abs(val) > sqrt(eps)
            Idx(ctr) = obj.nTestSpatial * obj.nTestTemporal + ndx1;
            Idy(ctr) = obj.nTestSpatial * obj.nTestTemporal + ndx2;
            Val(ctr) = val;
            ctr = ctr + 1;
          end
        end
      end

      % Assemble the sparse mass matrix
      M = sparse(Idy, Idx, Val, obj.nTestDim, obj.nTestDim);
    end % assembleTestNormMatrixSlow

    function M = assembleAnsatzNormMatrixSlow(obj)
      % Assemble the mass matrix of the discrete norm on the ansatz space.
      %
      % Slow implementation using numerical quadrature. Only useful for
      % comparison and testing of faster methods!
      %
      % Return values:
      %   M: mass matrix @type sparsematrix

      % Preparation for the sparse matrix
      Idx = ones(obj.nAnsatzDim, 1);
      Idy = ones(obj.nAnsatzDim, 1);
      Val = zeros(obj.nAnsatzDim, 1);
      ctr = 1;

      % iterate over spatial component
      for jdx1 = 1:obj.nAnsatzSpatial
        for jdx2 = 1:obj.nAnsatzSpatial
          % iterate over temporal component
          for kdx1 = 1:obj.nAnsatzTemporal
            for kdx2 = 1:obj.nAnsatzTemporal
              % first part `\norm{u}_{L_2(I; V)}`
              intSpatialF  = integral(@(x) obj.spatialBasisFunc(jdx1, x) .* obj.spatialBasisFunc(jdx2, x) + obj.spatialBasisFuncDerivative(jdx1, x) .* obj.spatialBasisFuncDerivative(jdx2, x), obj.xspan(1), obj.xspan(2));
              intTemporalF = integral(@(t) obj.temporalBasisFunc(kdx1, t) .* obj.temporalBasisFunc(kdx2, t), obj.tspan(1), obj.tspan(2));
              intFirstPart = intTemporalF * intSpatialF;

              % second part `norm{u_t}_{L_2(I; V')}`
              intSpatialS  = integral(@(x) obj.spatialBasisFunc(jdx1, x) .* obj.spatialBasisFunc(jdx2, x), obj.xspan(1), obj.xspan(2));
              intTemporalS = integral(@(t) obj.temporalBasisFuncDerivative(kdx1, t) .* obj.temporalBasisFuncDerivative(kdx2, t), obj.tspan(1), obj.tspan(2));
              intSecondPart = intTemporalS * intSpatialS;

              val = intFirstPart + intSecondPart;

              if obj.useNormalization
                val = val / (obj.normOfAnsatzFuncSlow(jdx1, kdx1) * ...
                  obj.normOfAnsatzFuncSlow(jdx2, kdx2));
              end

              % only consider values above a given threshold (since a lot of zero
              % valued integrals won't get exact zero)
              if abs(val) > sqrt(eps)
                Idx(ctr) = (kdx1 - 1) * obj.nAnsatzSpatial + jdx1;
                Idy(ctr) = (kdx2 - 1) * obj.nAnsatzSpatial + jdx2;
                Val(ctr) = val;
                ctr = ctr + 1;
              end
            end
          end
        end
      end

      % Assemble the sparse mass matrix
      M = sparse(Idy, Idx, Val, obj.nAnsatzDim, obj.nAnsatzDim);
    end % assembleAnsatzNormMatrixSlow

  end % methods
end % classdef
