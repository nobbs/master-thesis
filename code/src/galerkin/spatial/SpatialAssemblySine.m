classdef SpatialAssemblySine < SpatialAssemblyAbstract
  % Assemble the spatial mass and stiffness matrices for a sine basis. The
  % resulting spatial discretization has homogenuous boundary conditions.

  methods
    function obj = SpatialAssemblySine(pd, nbasis)
      % Constructor for this assembly class.
      %
      % Parameters:
      %   pd: Reference to the problem data object @type ProblemData
      %   nbasis: Number of basis functions to use @type integer

      % call the superclass constructor
      obj@SpatialAssemblyAbstract(pd, nbasis, nbasis, nbasis);
    end

    function Mx = massMatrix(obj)
      % Assemble the spatial mass matrix, that means we evaluate the integral
      % `\int_{\Omega} \sigma_j(x) \sigma_l(x) \diff x` for some `j` and `l`.
      %
      % Return values:
      %   Mx: spatial mass matrix @type matrix

      Mx = spdiags((obj.pd.xspan(2) - obj.pd.xspan(1)) * ...
                   ones(min(obj.nTrial, obj.nTest), 1) / 2, 0, obj.nTest, ...
                   obj.nTrial);
    end

    function Ax = stiffnessMatrix(obj)
      % Assemble the spatial stiffness matrix, that means we evaluate the
      % integral `\int_{\Omega} \sigma'_j(x) \sigma'_l(x) \diff x` for some `j`
      % and `l`.
      %
      % Return values:
      %   Ax: spatial stiffness matrix @type matrix
      %
      % @todo optimize

      tmp = zeros(min(obj.nTrial, obj.nTest), 1);
      for jdx = 1:min(obj.nTrial, obj.nTest)
        tmp(jdx)  = (pi * jdx)^2 / (2 * (obj.pd.xspan(2) - obj.pd.xspan(1)));
      end
      Ax = spdiags(tmp, 0, obj.nTest, obj.nTrial);
    end

    function val = evaluateFunctional(obj, fun, index)
      % Evaluate the functional given by a basis function.
      %
      % Evaluates the L2 inner product of a fixed basis function for the given
      % function fun.
      %
      % Parameters:
      %   fun: the argument of the functional @type function_handle
      %   index: index of the basis function to use @type integer
      %
      % Return values:
      %   val: value of the evaluated functional @type double

      % use one of the standard matlab numerical quadratures
      val = integral(@(z) fun(z) .* obj.basisFunc(index, z), ...
                     obj.pd.xspan(1), obj.pd.xspan(2));

      % round to zero if we are near enough to maintain sparsity
      if abs(val) < sqrt(eps)
        val = 0;
      end
    end

    function FDx = fieldDependentSine(obj)
      % Assemble only the field-dependet part of the stiffness matrix based on a
      % sine series expansion of the field.
      %
      % Return values:
      %   FDx: cellarray of the spatial matrices @type cell

      % first check if the given data is consistent
      assert(obj.pd.nC == length(obj.pd.seriesIdx));

      % create the cellarray
      FDx = cell(obj.pd.nC, 1);

      % iterate over the index of the sine series expansion functions
      for cdx = 1:obj.pd.nC
        bdx = obj.pd.seriesIdx(cdx);
        % create needed vectors to assemble the sparse matrix
        Idx = ones(obj.nTrial, 1);
        Idy = ones(obj.nTrial, 1);
        Val = zeros(obj.nTrial, 1);
        ctr = 1;

        % iterate over the indexes of the spatial basis functions of the trial
        % and test subspaces
        for jdx = 1:obj.nTrial
          for ldx = 1:obj.nTest
            % evaluate spatial integral
            % @todo explain optimization and cases
            intval = 0;
            % derived condition, if it's not satisfied, the value of the
            % integral is null
            if mod(bdx + ldx + jdx, 2) == 1
              % evaluate the spatial integral
              intval = ((obj.pd.xspan(2) - obj.pd.xspan(1)) / 2 ) * ...
                ( 1 / (pi * (bdx + jdx - ldx)) ...
                + 1 / (pi * (-bdx + jdx + ldx)) ...
                + 1 / (pi * (bdx - jdx + ldx)) ...
                - 1 / (pi * (bdx + jdx + ldx)));
            end

            % save the value
            Idx(ctr) = jdx;
            Idy(ctr) = ldx;
            Val(ctr) = intval;
            ctr = ctr + 1;
          end
        end

        % create the sparse matrix
        FDx{cdx} = sparse(Idy, Idx, Val, obj.nTest, obj.nTrial);
      end
    end

    function FDx = fieldDependentFourier(obj)
      % Assemble only the field-dependet part of the stiffness matrix based on a
      % fourier series expansion of the field (without the constant term).
      %
      % Return values:
      %   FDx: cellarray of the spatial matrices @type cell

      % first check if the given data is consistent
      assert(obj.pd.nC == length(obj.pd.seriesIdx));

      % create the cellarray
      FDx = cell(obj.pd.nC, 1);

      for cdx = 1:obj.pd.nC
        bdx = obj.pd.seriesIdx(cdx);
        % create needed vectors to assemble the sparse matrix
        Idx = ones(obj.nTrial, 1);
        Idy = ones(obj.nTrial, 1);
        Val = zeros(obj.nTrial, 1);
        ctr = 1;

        % iterate over the indexes of the spatial basis functions of the trial
        % and test subspaces
        for jdx = 1:obj.nTrial
          for ldx = 1:obj.nTest

            % evaluate the spatial integral. there are several combinations of
            % sine / cosine products we have to consider
            intval = 0;
            if mod(bdx, 2) == 0
              % bdx even: sine function sin(pi * bdx * x)
              if (mod(jdx, 2) == 0 && mod(ldx, 2) == 1) || ...
                (mod(jdx, 2) == 1 && mod(ldx, 2) == 0)
                intval = ((obj.pd.xspan(2) - obj.pd.xspan(1)) / 2 ) * ...
                  ( 1 / (pi * (bdx + jdx - ldx)) ...
                   + 1 / (pi * (-bdx + jdx + ldx)) ...
                   + 1 / (pi * (bdx - jdx + ldx)) ...
                   - 1 / (pi * (bdx + jdx + ldx)));
              end
            else
              % bdx odd: cosine function cos(pi * (bdx + 1) * x)
              val = 0;
              if - bdx + jdx + ldx - 1 == 0
                val = val - 1;
              end
              if bdx - jdx + ldx + 1 == 0
                val = val + 1;
              end
              if bdx + jdx - ldx + 1 == 0
                val = val + 1;
              end
              if bdx + jdx + ldx + 1 == 0
                val = val - 1;
              end
              intval = ((obj.pd.xspan(2) - obj.pd.xspan(1)) / 4) * val;
            end

            % save the value
            Idx(ctr) = jdx;
            Idy(ctr) = ldx;
            Val(ctr) = intval;
            ctr = ctr + 1;
          end
        end

        % create the sparse matrix
        FDx{cdx} = sparse(Idy, Idx, Val, obj.nTest, obj.nTrial);
      end
    end

    function val = basisFunc(obj, index, x)
      % Evaluate the basis function for a given index and x values.
      %
      % Parameters:
      %   index: index of the basis function @type integer
      %   x: values in which the function should be evaluated @type matrix
      %
      % Return values:
      %   val: values of the basis function in x @type matrix

      val = sin(pi * index * x / (obj.pd.xspan(2) - obj.pd.xspan(1)));
    end
  end % methods
end % classdef
