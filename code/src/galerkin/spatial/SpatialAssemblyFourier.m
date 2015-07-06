classdef SpatialAssemblyFourier < SpatialAssemblyAbstract
  % Assemble the spatial mass and stiffness matrices for a fourier basis. The
  % resulting spatial discretization has periodic boundary conditions.

  properties
    % witdh of the spatial interval @type double @type double
    xwidth;
  end

  methods

    function obj = SpatialAssemblyFourier(xwidth)
      % Constructor for this assembly class.
      %
      % Parameters:
      %   xwidth: width of the spatial interval @type double

      if nargin > 0
        obj.xwidth = xwidth;
      end
    end

    function Mx = massMatrix(obj, nX, nY)
      % Assemble the spatial mass matrix, that means we evaluate the integral
      % `\int_{\Omega} \sigma_j(x) \sigma_l(x) \diff x` for some `j` and `l`
      % where `\sigma` are Fourier basis functions.
      %
      % Parameters:
      %   nX: the limit of iteration for the `j` index @type integer
      %   nY: the limit of iteration for the `l` index @type integer
      %
      % Return values:
      %   Mx: spatial mass matrix @type matrix

      tmp    = obj.xwidth * ones(min(nX, nY), 1) / 2;
      tmp(1) = tmp(1) * 2;
      Mx     = spdiags(tmp, 0, nY, nX);
    end

    function Ax = stiffnessMatrix(obj, nX, nY)
      % Assemble the spatial stiffness matrix, that means we evaluate the
      % integral `\int_{\Omega} \sigma'_j(x) \sigma'_l(x) \diff x` for some `j`
      % and `l` where `\sigma` are Fourier basis functions.
      %
      % Parameters:
      %   nX: the limit of iteration for the `j` index @type integer
      %   nY: the limit of iteration for the `l` index @type integer
      %
      % Return values:
      %   Ax: spatial stiffness matrix @type matrix
      %
      % @todo optimize

      tmp = zeros(min(nX, nY), 1);
      for jdx = 2:min(nX, nY)
        if mod(jdx, 2) == 0
          tmp(jdx) = (pi * jdx)^2 / ( 2 * obj.xwidth);
        else
          tmp(jdx) = (pi * (jdx - 1))^2 / ( 2 * obj.xwidth);
        end
      end
      Ax = spdiags(tmp, 0, nY, nX);
    end

    function FDx = fieldDependentSine(obj, nX, nY, nC)
      % Assemble only the field-dependet part of the stiffness matrix based on a
      % sine series expansion of the field.
      %
      % Parameters:
      %   nX: the limit of iteration for the `j` index @type integer
      %   nY: the limit of iteration for the `l` index @type integer
      %   nC: number of sine basis functions to use @type integer
      %
      % Return values:
      %   FDx: cellarray of the spatial matrices @type cellarray
      %
      % @todo fixme

      % create the cellarray
      FDx = cell(nC, 1);

      % iterate over the index of the sine series expansion functions
      for cdx = 1:nC
        % create needed vectors to assemble the sparse matrix
        Idx = ones(nX, 1);
        Idy = ones(nX, 1);
        Val = zeros(nX, 1);
        ctr = 1;

        % iterate over the indexes of the spatial basis functions of the trial
        % and test subspaces
        for jdx = 1:nX
          for ldx = 1:nY
            % evaluate spatial integral
            % @todo explain optimization and cases
            intval = 0;
            if mod(jdx, 2) == 1 && mod(ldx, 2) == 1 && mod(cdx, 2) == 1
              intval = obj.xwidth * (cdx / (pi * (cdx^2 - (jdx - ldx)^2)) + cdx / (pi * (cdx^2 - (jdx + ldx - 2)^2)));
            elseif mod(jdx, 2) == 0 && mod(ldx, 2) == 0 && mod(cdx, 2) == 1
              intval = obj.xwidth * (cdx / (pi * (cdx^2 - (jdx - ldx)^2)) - cdx / (pi * (cdx^2 - (jdx + ldx)^2)));
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
              intval = (obj.xwidth / 4) * val;
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
              intval = (obj.xwidth / 4) * val;
            end

            % save the value
            Idx(ctr) = jdx;
            Idy(ctr) = ldx;
            Val(ctr) = intval;
            ctr = ctr + 1;
          end
        end

        % create the sparse matrix
        FDx{cdx} = sparse(Idy, Idx, Val, nY, nX);
      end
    end

    function FDx = fieldDependentFourier(obj, nX, nY, nC)
      % Assemble only the field-dependet part of the stiffness matrix based on a
      % fourier series expansion of the field (without the constant term).
      %
      % Parameters:
      %   nX: the limit of iteration for the `j` index @type integer
      %   nY: the limit of iteration for the `l` index @type integer
      %   nC: number of sine basis functions to use @type integer
      %
      % Return values:
      %   FDx: cellarray of the spatial matrices @type cellarray

      % create the cellarray
      FDx = cell(nC, 1);

      for cdx = 1:nC
        % create needed vectors to assemble the sparse matrix
        Idx = ones(nX, 1);
        Idy = ones(nX, 1);
        Val = zeros(nX, 1);
        ctr = 1;

        % iterate over the indexes of the spatial basis functions of the trial
        % and test subspaces
        for jdx = 1:nY
          for ldx = 1:nX

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
                intval = (obj.xwidth / 4) * val;
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
                intval = (obj.xwidth / 4) * val;
              end
            else
              % cdx odd: cosine function cos(pi * (cdx + 1) * x)
              if mod(jdx, 2) == 1 && mod(ldx, 2) == 1
                % j odd: cosine, l odd: cosine
                if jdx == 1 && ldx == 1
                  intval = 0;
                elseif jdx == 1 && ldx > 1 && cdx + 1 == ldx - 1
                  intval = obj.xwidth / 2;
                elseif jdx > 1 && ldx == 1 && cdx + 1 == jdx - 1
                  intval = obj.xwidth / 2;
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
                  intval = (obj.xwidth / 4) * val;
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
                intval = (obj.xwidth / 4) * val;
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
        FDx{cdx} = sparse(Idy, Idx, Val, nY, nX);
      end
    end

  end % methods

end % classdef
