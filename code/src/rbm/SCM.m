classdef SCM < handle
  % Successive constraint method.
  %
  % Allows the computation of a lower bound for the inf-sup-condition of a
  % variational problem.
  %
  % This class is an implementation of the successive constraint method
  % described in @cite Huynh2007 for the computation of lower and upper bounds
  % of the inf-sup-constant for a given parameter in an efficient offline /
  % online decomposition.
  
  properties
    % Verbosity toggle. Set it to true to get some messages about what scm is
    % currently doing. @type logical
    verbose = true;
    
    % Max error in each iteration @type vector
    errors;
  end
  
  properties(Access = 'protected')
    % Reference to the reduced basis solver object from which this object was
    % created @type RBMSolverAbstract
    rbm;
    % Holds the affine representation of the Y-norm of the supremizers @type
    % cell
    affineT;
    % Discrete norm of the truth trial space @type matrix
    normX;
    % storage of data generated in the offline stage @type struct
    offlineData;
    % Total number of addends in the affine representation of the supremizers
    % @type integer
    nQt;
    % Total number of addends in the affine representation of the underlying
    % variational problem @type integer
    nQb;
  end
  
  methods
    function obj = SCM(rbm)
      % Constructor for this class.
      %
      % Parameters:
      %   rbm: Reference to the reduced basis solver that created this scm
      %     object. @type RBMSolverAbstract
      
      % first set the reference
      obj.rbm = rbm;
      
      % and now prepare the scm offline stage
      obj.prepare();
    end
    
    % Main methods of the successive constraint method algorithm
    
    function exflag = offlinePhase(obj, paramTest, paramTrain, tol, maxIter)
      % Perform the offline greedy training stage.
      %
      % This method is the main offline stage of the successive constraint
      % method. The main part is a greedy style algorithm that looks for
      % the parameter with the largest gap between lower and upper bound
      % returned from the online solver with the current data set, then
      % solving the generalized eigenvalue problem for this parameter to
      % obtain the respective inf-sup-constant. The generated data of this
      % greedy algorithm can then be used in the online stage to compute
      % lower and upper bounds for any parameter with a certified error
      % bound.
      %
      % Possible values of exflag are:
      %   1: the largest remaining error is smaller than the given
      %      tolerance
      %   2: maximal number of iterations exceeded
      %   3: the training set is empty
      %
      % Parameters:
      %   paramTest: test parameter set which is used for the positivity
      %     constraints. this set should be really fine, like factor 10-100
      %     larger than paramTest. @type matrix
      %   paramTrain: training parameter set for the greedy algorithm. @type
      %     matrix
      %   tol: tolerance for the largest gap between certified upper and lower
      %     bound over all parameters in the training set @type double
      %     @default 1e-10
      %   maxIter: max number of greedy iterations @type integer @default 1000
      %
      % Return values:
      %   exflag: exit flag of the greedy algorithm. can be used to determine
      %     why the algorithm stopped. @type integer
      
      if obj.verbose
        fprintf('# Starting SCM offline stage\n');
        tic
      end
      
      % set defaults
      if nargin == 3
        tol = 1e-10;
        maxIter = 1000;
      elseif nargin == 4
        maxIter = 1000;
      end
      
      % defaults for the online phase
      Malpha = 15;
      Mplus  = 150;
      
      % create the structure that will hold the offline data
      obj.offlineData = struct();
      
      % calculate bounds for the variables in the linprog
      [lb, ub] = obj.getBounds();
      obj.offlineData.lowerBounds = lb;
      obj.offlineData.upperBounds = ub;
      
      % storage of the test parameters
      obj.offlineData.paramTest = paramTest;
      % convert and store the test parameters
      obj.offlineData.paramTestT = zeros(obj.nQt, size(paramTest, 2));
      for idx = 1:size(paramTest, 2)
        obj.offlineData.paramTestT(:, idx) = obj.mapParam(paramTest(:, idx));
      end
      
      % create the storage for the greedy-selected parameters and alle the
      % values we have to compute with them
      obj.offlineData.paramCk = [];
      obj.offlineData.paramCkT = [];
      obj.offlineData.ystarCk = [];
      obj.offlineData.alphaCk = [];
      
      % the first parameter is randomly selected from the training set
      idxTrain = randi(size(paramTrain, 2), 1, 1);
      curTrain = paramTrain(:, idxTrain);
      
      % final setup for the greedy algorithm
      isDone = false;
      exflag = 0;
      
      if obj.verbose
        fprintf('# Starting greedy loop ');
      end
      
      % start the greedy loop
      while ~isDone
        if obj.verbose
          % progress indicator, every dot represents one pass of this loop
          fprintf('.');
        end
        
        % assemble the system for the chosen parameter and calculate the
        % smallest eigenvalue and the corresponding eigenvector
        M = obj.assembleAffineT(curTrain);
        [mi, mivec] = obj.computeEV(M, obj.normX, 0);
        % if the eigenvalue mi is positive, then we got the largest eigenvalue
        % and have to do it again with a shift
        if mi >= 0
          [mi, mivec] = obj.computeEV(M, obj.normX, mi);
        end
        
        % calculate the values of ystar for the current parameter. this is done
        % by iterating over all the addends of the affine decomposition and
        % multiplication with the eigenvector from above from both sides
        ystar = zeros(obj.nQt, 1);
        for idx = 1:obj.nQt
          ystar(idx) = mivec.' * obj.affineT{idx} * mivec;
        end
        
        % save the newly generated data
        obj.offlineData.paramCk(:, end + 1)  = curTrain;
        obj.offlineData.alphaCk(end + 1)     = mi;
        obj.offlineData.ystarCk(:, end + 1)  = ystar;
        obj.offlineData.paramCkT(:, end + 1) = obj.mapParam(curTrain);
        
        % now we compute the relative gaps of upper and lower bound for all the
        % parameters given in the training set to find the best choice for the
        % next loop cycle
        relgap = zeros(size(paramTrain, 2), 1);
        for idx = 1:size(paramTrain, 2)
          [lb, ub] = obj.onlineSolve(paramTrain(:, idx), Malpha, Mplus, true);
          relgap(idx) = abs((ub - lb) / ub);
        end
        
        % get the parameter with the largest relgap
        [maxgap, maxdx]     = max(relgap);
        obj.errors(end + 1) = maxgap;
        
        % check the break conditions, maybe we are already done
        if maxgap < tol
          % we stop after successfully getting a largest remaining relative gap
          % smaller than the desired tolerance
          exflag = 1;
          isDone = true;
        elseif size(obj.offlineData.paramCk, 2) > maxIter
          % we stop because we performed the maxIter number of greedy cycles.
          % it's now maybe a good idea to check maxgap and decide whether the
          % offline stage should be restarted with a larger maxIter
          exflag = 2;
          isDone = true;
        elseif isempty(paramTrain) == 1
          % we stop because the training set is empty. this really shouldn't
          % happen if we want good bounds...
          exflag = 3;
          isDone = true;
        else
          % we aren't done yet, so let's choose the best parameter for the next
          % cycle and continue
          isDone = false;
          curTrain = paramTrain(:, maxdx);
          paramTrain(:, maxdx) = [];
        end
      end
      
      if obj.verbose
        t = toc;
        fprintf('done\n# SCM offline stage is done with exit flag %d after %f seconds.', exflag, t);
      end
    end
    
    function [lb, ub] = onlineSolve(obj, param, Malpha, Mplus, internal)
      % Online computation of the upper and lower bound.
      %
      % Calculates bounds for the coercivity constant of the modified
      % supremizers for a given parameter, which corresponds to the square of
      % the inf-sup-constant of the underlying variational problem.
      %
      % Attention:
      %   The parameters Malpha and Mplus control how many constraints are used
      %   for the linear program. It's likely that higher values ensure better
      %   bounds, so you may have to experiment to get a good tradeoff between
      %   accuracy and efficiency.
      %
      % Parameters:
      %   param: parameter for which we want the inf-sup-constant @type colvec
      %   Malpha: number of stability constraints to use @type integer
      %   Mplus: number of positivity constraints to use @type integer
      %   internal: flag, if the method is called internally from the offline
      %     stage @type logical @default false
      %
      % Return values:
      %   lb: lower bound for the coercivity or inf-sup-constant @type double
      %   ub: upper bound for the coercivity or inf-sup-constant @type double
      
      if nargin == 4
        internal = false;
      end
      
      % get the nearest parameters for stability and positivity constraints
      [~, Idxalpha] = obj.getNeighbors(Malpha, param, obj.offlineData.paramCk);
      [~, Idxplus]  = obj.getNeighbors(Mplus, param, obj.offlineData.paramTest);
      
      %| @todo as we are starting with the upper bounds as our initial value
      %|   it's maybe a good idea to perform a feasibility check
      
      % set up the linear problem:
      % start value are the upper bounds
      x0 = obj.offlineData.upperBounds;
      % objective vector
      f = obj.mapParam(param);
      % rhs and lhs of the inequalities
      b = [obj.offlineData.upperBounds;
        -obj.offlineData.lowerBounds;
        -obj.offlineData.alphaCk(Idxalpha).';
        zeros(size(Idxplus, 2), 1)];
      A = [eye(size(obj.offlineData.upperBounds, 1));
        -eye(size(obj.offlineData.lowerBounds, 1));
        -obj.offlineData.paramCkT(:, Idxalpha).';
        -obj.offlineData.paramTestT(:, Idxplus).'];
      
      % set some options for linprog
      opts = struct('LargeScale','off', ...
        'Algorithm', 'active-set', ...
        'Display', 'off');
      
      % disable the warning that active-set will be removed. may suppress
      % some other options-bad warnings...
      warning('off', 'optim:linprog:AlgOptsWillError');
      
      % now solve the linear problem
      [~, fval, flag, ~] = linprog(f, A, b, [],  [],  [],  [], x0, opts);
      
      % check the exit flag of linprog and break if it looks fishy
      if flag <= 0
        warning(['The exit flag of linprog is ', num2str(flag), '. ', ...
          'Looks like somethings wrong!', sprintf('\n'), ...
          'Falling back to console. Type return to continue.']);
        keyboard;
      end
      
      % the upper bound is simply the minimum over all objective values of y*
      if internal
        ub = min(obj.offlineData.ystarCk.' * f);
        lb = fval;
      else
        ub = sqrt(min(obj.offlineData.ystarCk.' * f));
        lb = sqrt(fval);
      end
    end
    
    function [lb, ub] = getBounds(obj)
      % Calculate the lower and upper bounds.
      %
      % This is done by solving for the largest and smallest eigenvalue of the
      % respective affine addend.
      %
      % Return values:
      %   lb: lower bounds @type colvec
      %   ub: upper bounds @type colvec
      
      if obj.verbose
        fprintf('# Computing bounds ');
      end
      
      lb = zeros(obj.nQt, 1);
      ub = zeros(obj.nQt, 1);
      
      % solve the respective generalized eigenvalue problem to obtain the
      % smallest and largest eigenvalue.
      for idx = 1:obj.nQt
        if obj.verbose
          fprintf('.');
        end
        
        [mi, ma] = obj.computeMinMaxEV(obj.affineT{idx}, obj.normX);
        ub(idx) = ma;
        lb(idx) = mi;
      end
      
      if obj.verbose
        fprintf(' done\n');
      end
    end
    
    function [Pm, Idx] = getNeighbors(~, M, param, C)
      % Return the neighbors of a given parameter.
      %
      % This method returns the M closest parameters to param in C according to
      % the euclidian distance.
      %
      % Parameters:
      %   M: number of neighbors to return @type integer
      %   param: parameter whose neighbors we want @type colvec
      %   C: set of parameters (each column is one parameter) @type matrix
      %
      % Return values:
      %   Pm: neighbors of the given parameter (again: each column is one
      %     parameter) @type matrix
      %   Idx: indexes of the returned paremeters in terms of columns of C
      %     @type vector
      
      % create the structures
      Pm = [];
      Idx = [];
      
      % if M is zero or the set C is empty, then there are no neighbors
      if M ~= 0 && isempty(C) ~= 1
        % total number of parameters in C
        nc = size(C, 2);
        
        % Calculate the euclidian distance to param for every parameter in C
        Dist  = sqrt(sum((repmat(param, 1, nc) - C).^2, 1));
        
        % return the M nearest neighbors
        [~, order] = sort(Dist);
        Idx        = order(1:min(nc, M));
        Pm         = C(:, Idx);
      end
    end
    
  end
  
  methods(Access = 'private')
    
    % Helper methods
    
    function newparam = mapParam(obj, param)
      % Remaps the parameter for the affine representation.
      %
      % This converts a parameter from the variational problem to a modified
      % parameter needed for the final assembly of the Y-norm of the supremizer.
      % For a motivation of this, see my master thesis.
      %
      % Parameters:
      %   param: parameter from the variational problem @type colvec
      %
      % Return values:
      %   newparam: parameter for the supremizer @type colvec
      
      % create the needed vector
      newparam = zeros(obj.nQt, 1);
      
      % pad the parameter from the variational problem with 1 for the field
      % independent part
      param = [1; shiftdim(param)];
      
      % and now we handle the assembly of the new parameter the same way we
      % handled the assembly of the addends in SCM@prepare.
      
      idx = 1;
      % first kind is handled first
      for qdx = 1:obj.rbm.nQb
        newparam(idx) = param(qdx) * (param(qdx) - sum(param(1:end ~= qdx)));
        idx = idx + 1;
      end
      
      % and now the simpler second kind
      for qdx = 1:obj.rbm.nQb
        for pdx = (qdx + 1):obj.rbm.nQb
          newparam(idx) = param(qdx) * param(pdx);
          idx = idx + 1;
        end
      end
    end
    
    function T = assembleAffineT(obj, param)
      % Assemble the supremizer for the given parameter.
      %
      % Parameters:
      %   param: parameter from the variational problem @type colvec
      %
      % Return values:
      %   T: supremizer for the given parameter @type matrix
      
      % remap the parameter
      param = obj.mapParam(param);
      
      % and sum up the addends of the affine representation
      T = sparse(obj.rbm.nTrialTruth, obj.rbm.nTrialTruth);
      for idx = 1:obj.nQt
        T = T + param(idx) * obj.affineT{idx};
      end
    end
    
    function [mi, ma] = computeMinMaxEV(obj, A, B)
      % Computes the minimal and maximal generalized eigenvalues.
      %
      % Parameters:
      %   A: matrix on the left side @type matrix
      %   B: matrix on the right side @type matrix
      %
      % Return values:
      %   mi: minimal eigenvalue @type double
      %   ma: maximal eigenvalue @type double
      
      % first we calculate the eigenvalue with the largest absolute value
      lm = obj.computeEV(A, B, 0);
      % if it's positive, then it's the maximal eigenvalue, else the
      % minimal eigenvalue
      if lm >= 0
        ma = lm;
        % now shift to get the minimal eigenvalue
        mi = obj.computeEV(A, B, lm);
      else
        mi = lm;
        % now shift to get the maximal eigenvalue
        ma = obj.computeEV(A, B, lm);
      end
    end
    
    function [lm, lmvec] = computeEV(~, A, B, shift)
      % Compute the largest generalized eigenvalue.
      %
      % This method computes the largest eigenvalue and the respective
      % eigenvector of the generalized eigenvalue problem `(\mat{A} +
      % \text{shift} * \mat{B})\vec{x} = \lambda \mat{B} \vec{x}`.
      %
      % Attention:
      %   There is some error handling included for the case that the matlab
      %   function eigs can't get the eigenvalue with default settings. So it
      %   may take some time, but it should get the eigenvalue most of the time,
      %   but possibly with a bad accuracy.
      %
      % Parameters:
      %   A: matrix on the left side of the eigenvalue problem @type matrix
      %   B: matrix on the right side of the eigenvalue problem @type matrix
      %   shift: shifts A on the left side by shift * B. this is needed if you
      %     want to compute the smallest eigenvalue; in this case you have to
      %     set shift to the value of the largest eigenvalue. @type double
      %
      % Return values:
      %   lm: largest eigenvalue @type double
      %   lmvec: eigenvector of the largest eigenvalue @type colvec
      %
      % @todo improve the error handling / accuracy decrease options
      
      % change the `eigenvalue won't converge` warning to an error, so that we
      % can catch it.
      origstate = warning;
      warning('error', 'MATLAB:eigs:NoEigsConverged'); %#ok<CTPCT>
      
      % setup the options structure
      opts = struct();
      opts.tol = eps;
      
      try
        % we first try to solve the eigenvalue problem with default settings.
        [lmvec, lm] = eigs(A - shift * B, B, 1, 'lm');
        % revert the shift
        lm = lm + shift;
      catch err
        % if the the default settings won't work, we have to try again with some
        % modifications
        if strcmp(err, 'MATLAB:eigs:NoEigsConverged') || ...
            strcmp(err.identifier, 'MATLAB:eigs:ARPACKroutineErrorMinus14')
          % flag to check whether we found a eigenvalue
          hasEV = false;
          % decrease the accuracy
          opts.tol = opts.tol * 100;
          % set the number of used Lanczos vectors
          opts.p = 50;
          while ~hasEV
            try
              % and now try again
              [lmvec, lm] = eigs(A - shift * B, B, 1, 'lm', opts);
              lm = lm  + shift;
              hasEV = true;
            catch err2
              hasEV = false;
              if strcmp(err, 'MATLAB:eigs:NoEigsConverged') || ...
                  strcmp(err2.identifier, 'MATLAB:eigs:ARPACKroutineErrorMinus14')
                % decrease the accuracy
                opts.tol = opts.tol * 100;
                % set the number of used Lanczos vectors
                opts.p = 50;
              else
                % other errors should be passed back
                rethrow(err2);
              end
            end
          end
        else
          % other errors should be passed back
          rethrow(err);
        end
      end
      
      % if we found an eigenvalue that looks like it's zero, it probably is
      % zero, so ...
      if abs(lm) < sqrt(eps)
        lm = 0;
      end
      
      % revert the warning state changes
      warning(origstate);
    end
    
    function prepare(obj)
      % Prepare the SCM object for the offline phase.
      %
      % This mainly consists of the creation of the addends for the affine
      % representation of the Y-norm of the supremizers.
      
      % Create local copies of the truth solver structures.
      obj.normX   = obj.rbm.solver.TrNorm;
      % more truth solver structures, this time only needed for the preparation
      normY   = obj.rbm.solver.TeNorm;
      affineB = obj.rbm.solver.Lhs;
      
      % set the total number of needed addends for the affine representation
      obj.nQt = obj.rbm.nQb * (obj.rbm.nQb + 1) / 2;
      
      % create the needed structure
      obj.affineT = cell(obj.nQt, 1);
      
      % iterate over the different addends of the affine representation. as we
      % are using an alternative expansion of the field independent operators to
      % ensure that they are all symmetric and positiv semi-definite, we have to
      % deal with the two following types of addends:
      %  1. <T_q u, T_q u>_Y for q = 1..Qb,
      %  2. <T_q u + T_p u, T_q u + T_p u>_Y for q = 1..Qb and p = (q+1)..Qb.
      
      % index for the affineT structure (a mapping of the index pair (q, p) to a
      % useful 1d index would be nice, but... meh).
      idx = 1;
      
      % assemble the affine addends of the first kind
      for qdx = 1:obj.rbm.nQb
        obj.affineT{idx} = affineB{qdx}.' * (normY \ affineB{qdx});
        idx = idx + 1;
      end
      
      % and now assemble the affine addends of the second kind
      for qdx = 1:obj.rbm.nQb
        for pdx = (qdx + 1):obj.rbm.nQb
          obj.affineT{idx} = (affineB{qdx} + affineB{pdx}).' * ...
            (normY \ (affineB{qdx} + affineB{pdx}));
          idx = idx + 1;
        end
      end
    end
  end
end
