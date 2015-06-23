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
    % nothing to see here
  end
  
  methods
    
    %% Fast implementation of the assembly methods
    
    function M = assembleFieldIndependentMatrix(obj)
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
              intTemporal = 2;% / (obj.tspan(2) - obj.tspan(1));
              
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
            intSpatial = (pi * jdx)^2 / ( 2 * obj.xspan(2));
          else
            intSpatial = (pi * (jdx - 1))^2 / ( 2 * obj.xspan(2));
          end
          
          % evaluate temporal integral
          intTemporal = (obj.tspan(2) - obj.tspan(1)) / (2 * (kdx - 1) + 1);
          
          Idx(ctr) = (jdx - 1) * obj.nAnsatzTemporal + kdx;
          Idy(ctr) = (jdx - 1) * obj.nTestTemporal + kdx;
          Val(ctr) = obj.coeffLaplacian * intTemporal * intSpatial;
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
      O = {};
      
      % iterate over the index of the sine basis functions (of the sine series
      % expansion of the field).
      for cdx = 1:nCoeff
        % create needed vectors to assemble the sparse matrix
        Idx = ones(obj.dAnsatz, 1);
        Idy = ones(obj.dAnsatz, 1);
        Val = zeros(obj.dAnsatz, 1);
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
              Val(ctr) = intTemporal * intSpatial;;
              ctr = ctr + 1;
            end
          end
        end
        
        % create the sparse matrix
        O{cdx} = sparse(Idy, Idx, Val, obj.dTest, obj.dAnsatz);;
      end
    end
    
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
      O = {};
      
      % iterate over the index of the Fourier basis function (of the Fourier
      % series expansion of the field).
      % order of functions:
      % 1: cos(2pix), 2: sin(2pix), 3: cos(4pix) ...
      % cdx = i odd: cos((i+1) pi x), cdx = i even: sin(i pi x)
      for cdx = 1:nCoeff
        % create needed vectors to assemble the sparse matrix
        Idx = ones(obj.dAnsatz, 1);
        Idy = ones(obj.dAnsatz, 1);
        Val = zeros(obj.dAnsatz, 1);
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
              Val(ctr) = intTemporal * intSpatial;;
              ctr = ctr + 1;
            end
          end
        end
        
        % create the sparse matrix
        O{cdx} = sparse(Idy, Idx, Val, obj.dTest, obj.dAnsatz);
      end
    end
    
    
    %% Assembly of the load vector for different kinds of data
    
    function F = assembleVectorFromCoeffs(obj, coeffs)
      % Assemble the load vector for no source and initial data given by
      % coeffients corresponding to the spatial basis functions.
      %
      % Parameters:
      %   coeffs: coefficients corresponding to spatial basis functions @type vector
      %
      % Return values:
      %   F: load vector @type colvec
      
      F = zeros(obj.dTest, 1);
      
      % iterate over the spatial basis functions corresponding to the initial
      % condition
      for ndx = 1:obj.nTestSpatialIC
        
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
    end
    
    function F = assembleVectorOnes(obj)
      % Assemble the load vector for no source and initial data equal to a
      % constant value of one.
      %
      % Return values:
      %   F: load vector @type colvec
      
      F = zeros(obj.dTest, 1);
      F(obj.nTestSpatial * obj.nTestTemporal + 1) = obj.xspan(2);
    end
    
    
    %% Slow implementations of the assembly methods
    
    function M = assembleFieldIndependentMatrixSlow(obj)
      % Assemble the field-independent part of the stiffness matrix. Slow
      % implementation, only for comparison!
      %
      % This is done by numerical integration of all the occurring products of
      % spatial and temporal basis functions, so it's kinda slow and not that
      % accurate. Only really useful to test the faster, optimized version.
      %
      % Return values:
      %   M: Field-independent part of the stiffness matrix
      %     @type sparsematrix
      
      % Preparation for the sparse matrix
      Idx = ones(obj.dAnsatz, 1);
      Idy = ones(obj.dAnsatz, 1);
      Val = zeros(obj.dAnsatz, 1);
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
              
              % only consider values above a given threshold (since a lot of
              % zero valued integrals won't get exact zero)
              if abs(val) > sqrt(eps)
                Idx(ctr) = (jdx - 1) * obj.nAnsatzTemporal + kdx;
                Idy(ctr) = (ldx - 1) * obj.nTestTemporal + mdx;
                Val(ctr) = val;
                ctr = ctr + 1;
              end
            end
          end
        end
      end
      
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
              
              % only consider values above a given threshold (since a lot of
              % zero valued integrals won't get exact zero)
              if abs(intSpatial * intTemporal) > sqrt(eps)
                Idx(ctr) = (jdx - 1) * obj.nAnsatzTemporal + kdx;
                Idy(ctr) = (ldx - 1) * obj.nTestTemporal + mdx;
                Val(ctr) = obj.coeffLaplacian * intSpatial * intTemporal;
                ctr = ctr + 1;
              end
            end
          end
        end
      end
      
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
              
              % only consider values above a given threshold (since a lot of
              % zero valued integrals won't get exact zero)
              if abs(val) > sqrt(eps)
                Idx(ctr) = (jdx - 1) * obj.nAnsatzTemporal + kdx;
                Idy(ctr) = (ldx - 1) * obj.nTestTemporal + mdx;
                Val(ctr) = obj.coeffOffset * val;
                ctr = ctr + 1;
              end
            end
          end
        end
      end
      
      % Calculate the fourth term `\skp{u(0)}{v_2}{L_2(\Omega)}`.
      for jdx = 1:obj.nAnsatzSpatial
        for kdx = 1:obj.nAnsatzTemporal
          for ndx = 1:obj.nTestSpatialIC
            val = obj.temporalBasisFunc(kdx, obj.tspan(1)) * ...
              integral(@(x) obj.spatialBasisFunc(jdx, x) .* ...
              obj.spatialBasisFunc(ndx, x), obj.xspan(1), obj.xspan(2));
            
            % only consider values above a given threshold (since a lot of zero
            % valued integrals won't get exact zero)
            if abs(val) > sqrt(eps)
              Idx(ctr) = (jdx - 1) * obj.nAnsatzTemporal + kdx;
              Idy(ctr) = obj.nTestSpatial * obj.nTestTemporal + jdx;
              Val(ctr) = val;
              ctr = ctr + 1;
            end
          end
        end
      end
      
      % Finally assemble the stiffness matrix
      M = sparse(Idy, Idx, Val, obj.dTest, obj.dAnsatz);
    end
    
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
        Idx = ones(obj.dAnsatz, 1);
        Idy = ones(obj.dAnsatz, 1);
        Val = zeros(obj.dAnsatz, 1);
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
                  
                  if abs(intTemporal) > sqrt(eps)
                    % save the evaluated integrals
                    Idx(ctr) = (jdx - 1) * obj.nAnsatzTemporal + kdx;
                    Idy(ctr) = (ldx - 1) * obj.nTestTemporal + mdx;
                    Val(ctr) = intTemporal * intSpatial;;
                    ctr = ctr + 1;
                  end
                end
              end
            end
          end
        end
        
        % create the sparse matrix
        O{cdx} = sparse(Idy, Idx, Val, obj.dTest, obj.dAnsatz);
      end
    end
    
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
        Idx = ones(obj.dAnsatz, 1);
        Idy = ones(obj.dAnsatz, 1);
        Val = zeros(obj.dAnsatz, 1);
        ctr = 1;
        
        % iterate over the indexes of the spatial basis functions of the ansatz
        % and test subspaces
        for jdx = 1:obj.nAnsatzSpatial
          for ldx = 1:obj.nTestSpatial
            % evaluate the spatial integral. there are several combinations of
            % sine / cosine products we have to consider
            intSpatial = 0;
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
                  
                  if abs(intTemporal) > sqrt(eps)
                    % save the evaluated integrals
                    Idx(ctr) = (jdx - 1) * obj.nAnsatzTemporal + kdx;
                    Idy(ctr) = (ldx - 1) * obj.nTestTemporal + mdx;
                    Val(ctr) = intTemporal * intSpatial;;
                    ctr = ctr + 1;
                  end
                end
              end
            end
          end
        end
        
        % create the sparse matrix
        O{cdx} = sparse(Idy, Idx, Val, obj.dTest, obj.dAnsatz);
      end
    end
    
    
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
    end
    
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
      %
      % Return values:
      %   val: values of the basis function in t @type matrix
      
      val = legendrePolynomial(t, index - 1, obj.tspan);
    end
    
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
    
    function val = temporalBasisFuncDerivative(obj, index, t)
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
      %
      % Return values:
      %   val: values of the basis function in t @type matrix
      
      val = legendrePolynomialDerivative(t, index - 1, obj.tspan);
    end
    
  end
  
end
