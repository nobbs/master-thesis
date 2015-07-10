classdef SpatialAssemblyFourier < SpatialAssemblyAbstract
  % Assemble the spatial mass and stiffness matrices for a fourier basis. The
  % resulting spatial discretization has periodic boundary conditions.

  methods
    function obj = SpatialAssemblyFourier(pd, nbasis)
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

      tmp    = (obj.pd.xspan(2) - obj.pd.xspan(1)) * ones(min(obj.nTrial, obj.nTest), 1) / 2;
      tmp(1) = tmp(1) * 2;
      Mx     = spdiags(tmp, 0, obj.nTest, obj.nTrial);
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
      for jdx = 2:min(obj.nTrial, obj.nTest)
        if mod(jdx, 2) == 0
          tmp(jdx) = (pi * jdx)^2 / ( 2 * (obj.pd.xspan(2) - obj.pd.xspan(1)));
        else
          tmp(jdx) = (pi * (jdx - 1))^2 / ( 2 * (obj.pd.xspan(2) - obj.pd.xspan(1)));
        end
      end
      Ax = spdiags(tmp, 0, obj.nTest, obj.nTrial);
    end

    function FDx = fieldDependentSine(obj)
      % Assemble only the field-dependet part of the stiffness matrix based on a
      % sine series expansion of the field.
      %
      % Return values:
      %   FDx: cellarray of the spatial matrices @type cell

      % create the cellarray
      FDx = cell(obj.pd.nC, 1);

      % iterate over the index of the sine series expansion functions
      for cdx = 1:obj.pd.nC
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
            if mod(jdx, 2) == 1 && mod(ldx, 2) == 1 && mod(cdx, 2) == 1
              intval = (obj.pd.xspan(2) - obj.pd.xspan(1)) * ...
                (cdx / (pi * (cdx^2 - (jdx - ldx)^2)) ...
                 + cdx / (pi * (cdx^2 - (jdx + ldx - 2)^2)));
            elseif mod(jdx, 2) == 0 && mod(ldx, 2) == 0 && mod(cdx, 2) == 1
              intval = (obj.pd.xspan(2) - obj.pd.xspan(1)) * ...
                (cdx / (pi * (cdx^2 - (jdx - ldx)^2)) ...
                 - cdx / (pi * (cdx^2 - (jdx + ldx)^2)));
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
              intval = ((obj.pd.xspan(2) - obj.pd.xspan(1)) / 4) * val;
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

    function FDx = fieldDependentFourier(obj)
      % Assemble only the field-dependet part of the stiffness matrix based on a
      % fourier series expansion of the field (without the constant term).
      %
      % Return values:
      %   FDx: cellarray of the spatial matrices @type cell

      % create the cellarray
      FDx = cell(obj.pd.nC, 1);

      for cdx = 1:obj.pd.nC
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

            if mod(cdx, 2) == 0
              % cdx even: sine function sin(pi * cdx * x)
              if mod(jdx, 2) == 1 && mod(ldx, 2) == 1
                % j odd: cosine, l odd: cosine
                % intval = 0;
              elseif mod(jdx, 2) == 0 && mod(ldx, 2) == 0
                % j even: sine, l even: sine
                % intval = 0;
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
                intval = ((obj.pd.xspan(2) - obj.pd.xspan(1)) / 4) * val;
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
                intval = ((obj.pd.xspan(2) - obj.pd.xspan(1)) / 4) * val;
              end
            else
              % cdx odd: cosine function cos(pi * (cdx + 1) * x)
              if mod(jdx, 2) == 1 && mod(ldx, 2) == 1
                % j odd: cosine, l odd: cosine
                if jdx == 1 && ldx == 1
                  intval = 0;
                elseif jdx == 1 && ldx > 1 && cdx + 1 == ldx - 1
                  intval = (obj.pd.xspan(2) - obj.pd.xspan(1)) / 2;
                elseif jdx > 1 && ldx == 1 && cdx + 1 == jdx - 1
                  intval = (obj.pd.xspan(2) - obj.pd.xspan(1)) / 2;
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
                  intval = ((obj.pd.xspan(2) - obj.pd.xspan(1)) / 4) * val;
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
                intval = ((obj.pd.xspan(2) - obj.pd.xspan(1)) / 4) * val;
              elseif mod(jdx, 2) == 0 && mod(ldx, 2) == 1
                % j even: sine, l odd: cosine
                % intval = 0;
              elseif mod(jdx, 2) == 1 && mod(ldx, 2) == 0
                % j odd: cosine, l even: sine
                % intval = 0;
              end
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

      if index == 1
        val = ones(size(x, 1), size(x, 2));
      elseif mod(index, 2) == 0
        val = sin(pi * index * x / (obj.pd.xspan(2) - obj.pd.xspan(1)));
      else
        val = cos(pi * (index - 1) * x / (obj.pd.xspan(2) - obj.pd.xspan(1)));
      end
    end
  end % methods
end % classdef
