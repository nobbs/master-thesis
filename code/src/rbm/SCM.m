classdef SCM < handle
  % Successive constraint method.
  %
  % Allows the computation of a lower bound for the inf-sup-condition of a
  % variational problem.
  %
  % This class is an implementation of the successive constraint method
  % described in @cite Huynh2007 for the computation of lower and upper bounds
  % of the inf-sup-constant for a given parameter in an efficient offline /
  % online decomposition with some improvements from @cite Chen2009.
  %
  % Warning:
  %   The Access tags are commented so that it's easier to check the data the
  %   scm generated, but you still shouldn't mess with them.

  properties
    % Max error in each iteration @type vector
    errors;
  end

  properties%(Access = 'protected')
    % Reference to the "truth" solver object @type SolverAbstract
    solver;
    % Holds the affine representation of the Y-norm of the supremizers @type
    % cell
    affineT;
    % Discrete norm of the truth trial space @type matrix
    normX;
    % Total number of addends in the affine representation of the supremizers
    % @type integer
    nQt;
    % Number of stability constraints to use @type integer
    Malpha;
    % Number of positivity constraints to use @type integer
    Mplus;
    % Number of training parameters @type integer
    nMuTrain;

    % Helper vector @type vector
    MuIndexes;
    % Logger handle @type log4m
    L;
  end

  properties%(Access = 'protected') % stop and go support
    % Holds the current upper bounds for every parameter in the training set.
    % @type vector
    muTrainUB;
    % Holds the current lower bounds for every parameter in the training set.
    % @type vector
    muTrainLB;
    % Index of the parameter chosen in the previous to last greedy loop. @type
    % integer
    muPrevIdx;
    % Index of the parameter chosen in the last greedy loop. @type integer
    muCurIdx;
    % Total number of cycles done by the greedy loop. @type integer
    totalLoopCtr;
  end

  properties%(Access = 'protected') % Offline Data
    % Lower bounds for the linear program @type colvec
    offLowerBounds;
    % Upper bounds for the linear program @type colvec
    offUpperBounds;

    % parameter training set @type matrix
    offMuTrain;
    % parameter training set remapped for the non-coercive case @type matrix
    offMuTrainMapped;

    % parameter set Ck @type matrix
    offMuCk;
    % parameter set Ck remapped @type matrix
    offMuCkMapped;
    % parameter set Ck indexes @type vector
    offMuCkIndex;
    % parameter set Ck y* values @type matrix
    offYstarCk;
    % parameter set Ck alpha values @type vector
    offAlphaCk;

    % lower bounds for the current alpha values @type vector
    offMuLB;
    % boolean that marks parameters of the training set in Ck @type vector
    offMuChosen;
  end

  properties % for debugging and analysis of the scm
    save_scm_bounds = false;
    scm_lower_bounds;
    scm_upper_bounds;
  end

  methods
    function obj = SCM(solver, muTrain, Malpha, Mplus)
      % Constructor for this class.
      %
      % Parameters:
      %   solver: Reference to the "truth" solver object @type SolverAbstract
      %   muTrain: training parameter set for the greedy algorithm. @type
      %     matrix
      %   Malpha: number of stability constraints to use in the linear program
      %     @type integer @default 30
      %   Mplus: number of positivity constraints to use in the linear program
      %     @type integer @default 10

      % get a logger reference
      obj.L = log4m.getLogger('scftrbm.log');

      % default values
      if nargin < 2
        error('');
      else
        if ~exist('Malpha', 'var'),  Malpha  = 30; end;
        if ~exist('Mplus', 'var'),   Mplus   = 10; end;
      end
      obj.Malpha = Malpha;
      obj.Mplus  = Mplus;

      % first set the reference
      obj.solver = solver;

      % and now prepare the scm offline stage
      obj.prepare();

      % save the training set
      obj.offMuTrain = muTrain;

      if obj.save_scm_bounds
        obj.scm_lower_bounds = {};
        obj.scm_upper_bounds = {};
      end
    end

    % Main methods of the successive constraint method algorithm

    function [lb, ub] = onlineQuery(obj, param, internal)
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
      %   internal: flag, if the method is called internally from the offline
      %     stage @type logical @default false
      %
      % Return values:
      %   lb: lower bound for the coercivity or inf-sup-constant @type double
      %   ub: upper bound for the coercivity or inf-sup-constant @type double

      if nargin == 2
        internal = false;
      end

      % get the neighbors from C_k
      [~, Idxalpha] = obj.getNeighbors(obj.Malpha, param, obj.offMuCk);

      % get the neighbors from Xi \setminus C_k,
      % first we have to remap the indexes back to the training set
      paramsNotInCk                   = true(obj.nMuTrain, 1);
      paramsNotInCk(obj.offMuCkIndex) = false;
      paramsNotInCkIndexes = obj.MuIndexes(paramsNotInCk);
      [~, IdxplusXi]       = obj.getNeighbors(obj.Mplus, param, obj.offMuTrain(:, paramsNotInCk));
      Idxplus              = paramsNotInCkIndexes(IdxplusXi);

      % set up the linear problem:
      % start value are the upper bounds
      x0 = obj.offUpperBounds;
      % objective vector
      f = obj.mapParam(param);
      % rhs and lhs of the inequalities
      b = [obj.offUpperBounds;
        -obj.offLowerBounds;
        -obj.offAlphaCk(Idxalpha).';
        -obj.offMuLB(Idxplus)];
      A = [eye(size(obj.offUpperBounds, 1));
        -eye(size(obj.offLowerBounds, 1));
        -obj.offMuCkMapped(:, Idxalpha).';
        -obj.offMuTrainMapped(:, Idxplus).'];

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

      % we compute the square of the wanted inf-sup-constant, but don't take the
      % square root if we are in the offline stage. the upper bound is simply
      % the minimum over all objective values of y*
      if internal
        ub = min(obj.offYstarCk.' * f);
        lb = fval;
      else
        ub = sqrt(min(obj.offYstarCk.' * f));
        lb = sqrt(fval);
      end
    end

    function exflag = startOfflineStage(obj, tol, maxIter)
      % Perform the offline greedy training stage.
      %
      % This method sets up the greedy loop by performing the calculations for
      % the first randomly selected parameter and consequent creation of all the
      % needed structures. It then calls the greedyLoop method.
      %
      % Possible values of exflag are:
      %   1: the largest remaining error is smaller than the given
      %      tolerance
      %   2: maximal number of iterations exceeded
      %   3: the training set is empty
      %
      % Parameters:
      %   tol: tolerance for the largest gap between certified upper and lower
      %     bound over all parameters in the training set @type double
      %     @default 1e-3
      %   maxIter: max number of greedy iterations @type integer @default 1e5
      %
      % Return values:
      %   exflag: exit flag of the greedy algorithm. can be used to determine
      %     why the algorithm stopped. @type integer

      % default values
      if nargin < 1
        error('');
      else
        if ~exist('tol', 'var'),     tol     = 1e-2; end;
        if ~exist('maxIter', 'var'), maxIter = 1e5; end;
      end

      % logging
      obj.L.info('SCM', ...
        sprintf(' Starting offline stage with %d max iterations and tolerance %.2e', ...
          maxIter, tol));

      % get some needed data from the training set
      obj.nMuTrain  = size(obj.offMuTrain, 2);
      obj.MuIndexes = 1:obj.nMuTrain;

      % convert and store the training parameters
      obj.offMuTrainMapped = zeros(obj.nQt, obj.nMuTrain);
      for idx = 1:obj.nMuTrain
        obj.offMuTrainMapped(:, idx) = obj.mapParam(obj.offMuTrain(:, idx));
      end

      % calculate bounds for the variables in the linear program
      [lb, ub] = obj.getBounds();
      obj.offLowerBounds = lb;
      obj.offUpperBounds = ub;

      % create the storage for the greedy-selected parameters and alle the
      % values we have to compute with them
      obj.offMuCk       = [];
      obj.offMuCkMapped = [];
      obj.offMuCkIndex  = [];
      obj.offYstarCk    = [];
      obj.offAlphaCk    = [];
      obj.offMuChosen   = false(obj.nMuTrain, 1);
      obj.offMuLB       = zeros(obj.nMuTrain, 1);

      % these two vectors will hold the current upper and lower bounds for
      % every parameter in the training set
      obj.muTrainUB = zeros(obj.nMuTrain, 1);
      obj.muTrainLB = zeros(obj.nMuTrain, 1);

      % the first parameter is randomly selected from the training set
      obj.muPrevIdx = -1;
      obj.muCurIdx  = randi(obj.nMuTrain, 1, 1);

      % calculate stuff for the first parameter. assemble the system for the
      % chosen parameter and calculate the smallest eigenvalue and the
      % corresponding eigenvector
      M = obj.assembleAffineT(obj.offMuTrainMapped(:, obj.muCurIdx));
      [mi, mivec] = obj.computeEV(M, obj.normX, 0);
      % if the eigenvalue mi is positive, then we got the largest eigenvalue
      % and have to do it again with a shift
      if mi >= 0
        [mi, mivec] = obj.computeEV(M, obj.normX, mi);
      end

      % the eigenvalue is lower and upper bound for this parameter mu
      obj.muTrainLB(obj.muCurIdx) = mi;
      obj.muTrainUB(obj.muCurIdx) = mi;

      % calculate the values of ystar for the current parameter. this is done
      % by iterating over all the addends of the affine decomposition and
      % multiplication with the eigenvector from above from both sides
      ystar = zeros(obj.nQt, 1);
      for idx = 1:obj.nQt
        ystar(idx) = mivec.' * obj.affineT{idx} * mivec;
      end

      % save the newly generated data
      obj.offMuCk(:, 1)         = obj.offMuTrain(:, obj.muCurIdx);
      obj.offMuCkMapped(:, 1)   = obj.offMuTrainMapped(:, obj.muCurIdx);
      obj.offMuCkIndex(1)       = obj.muCurIdx;
      obj.offYstarCk(:, 1)      = ystar;
      obj.offAlphaCk(1)         = mi;
      obj.offMuChosen(obj.muCurIdx) = true;

      % start the greedy loop
      exflag = obj.greedyLoop(tol, maxIter);
    end

    function exflag = continueOfflineStage(obj, tol, maxIter)
      % Continue with the offline stage.
      %
      % This method allows to "restart" the greedy loop without deletion of the
      % already obtained data.
      %
      % Possible values of exflag are:
      %   1: the largest remaining error is smaller than the given
      %      tolerance
      %   2: maximal number of iterations exceeded
      %   3: the training set is empty
      %
      % Parameters:
      %   tol: tolerance for the largest gap between certified upper and lower
      %     bound over all parameters in the training set @type double
      %     @default 1e-10
      %   maxIter: max number of greedy iterations @type integer @default 1000
      %
      % Return values:
      %   exflag: exit flag of the greedy algorithm. can be used to determine
      %     why the algorithm stopped. @type integer

      % default values
      if nargin < 1
        error('');
      else
        if ~exist('tol', 'var'),     tol     = 1e-2; end;
        if ~exist('maxIter', 'var'), maxIter = 1e5; end;
      end

      % logging
      obj.L.info('SCM', ...
        sprintf(' Continuation of offline stage with %d max iterations and tolerance %.2e (already %d cycles done)', ...
          maxIter, tol, obj.totalLoopCtr));

      exflag = obj.greedyLoop(tol, maxIter);
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

      % logging
      obj.L.info('SCM', ...
        sprintf(' Computing %d pairs of upper / lower bounds', obj.nQt));

      lb = zeros(obj.nQt, 1);
      ub = zeros(obj.nQt, 1);

      % solve the respective generalized eigenvalue problem to obtain the
      % smallest and largest eigenvalue.
      for idx = 1:obj.nQt
        [mi, ma, retries] = obj.computeMinMaxEV(obj.affineT{idx}, obj.normX);
        ub(idx) = ma;
        lb(idx) = mi;

        % logging
        obj.L.info('SCM', ...
          sprintf(' solved eigenvalue problem %3d of %3d (after %2d retries): [%f, %f]', idx, obj.nQt, retries, mi, ma));
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

    function exflag = greedyLoop(obj, tol, maxIter)
      % Start the greedy loop.
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
      %   tol: tolerance for the largest gap between certified upper and lower
      %     bound over all parameters in the training set @type double
      %     @default 1e-10
      %   maxIter: max number of greedy iterations @type integer @default 1000
      %
      % Return values:
      %   exflag: exit flag of the greedy algorithm. can be used to determine
      %     why the algorithm stopped. @type integer

      % final setup for the greedy algorithm
      isDone  = false;
      exflag  = 0;
      curLoopCtr = 1;

      obj.L.info('SCM', ...
        sprintf(' Starting greedy loop with %d max iterations and tolerance %.2e', ...
          maxIter, tol));

      % Now we can finally start the greedy loop
      while ~isDone
        scmtic = tic;
        % now we compute the relative gaps of upper and lower bound for all the
        % parameters given in the training set to find the best choice for the
        % next loop cycle
        reLBStability  = zeros(obj.nMuTrain, 1);
        reLBPositivity = zeros(obj.nMuTrain, 1);

        disp('checking parameters: ');
        muTrainLB = obj.muTrainLB;
        muTrainUB = obj.muTrainUB;
        parfor idx = 1:obj.nMuTrain
          %textprogressbar(idx * 100 / obj.nMuTrain);
          if ~obj.offMuChosen(idx)
            param = obj.offMuTrain(:, idx);

            % set up the needed sets for the positivity constraint calculation
            paramsNotInCkPrev = true(obj.nMuTrain, 1);
            paramsNotInCkPrev(obj.offMuCkIndex(1:(end-1))) = false;
            paramsNotInCkPrevIndexes = obj.MuIndexes(paramsNotInCkPrev);

            % get the nearest parameters for stability and positivity constraints
            [~, Idxalpha]  = obj.getNeighbors(obj.Malpha, param, obj.offMuCk);
            [~, IdxplusXi] = obj.getNeighbors(obj.Mplus, param, obj.offMuTrain(:, paramsNotInCkPrev));
            Idxplus        = paramsNotInCkPrevIndexes(IdxplusXi);

            % check whether the lower bound could change
            reLBStability(idx)  = any(obj.offMuCkIndex(Idxalpha) == obj.muCurIdx);
            reLBPositivity(idx) = any(Idxplus == obj.muCurIdx);

            % if one of the above conditions holds, then we compute new
            % (hopefully better) lower bounds, if not, there is no improvement
            % and we take the value from the last loop
            if reLBStability(idx) || reLBPositivity(idx)
              [lb, ub] = obj.onlineQuery(param, true);
              muTrainLB(idx) = lb;
              muTrainUB(idx) = ub;
            else
              muTrainLB(idx) = obj.offMuLB(idx);
              muTrainUB(idx) = min(obj.offYstarCk.' * obj.offMuTrainMapped(:, idx));
            end
          end
        end
        obj.muTrainLB = muTrainLB;
        obj.muTrainUB = muTrainUB;
        disp(' done!');

        % save the newly calculated lower bounds for online use
        obj.offMuLB = obj.muTrainLB;

        if obj.save_scm_bounds
          obj.scm_lower_bounds{obj.totalLoopCtr + 1} = sqrt(obj.muTrainLB);
          obj.scm_upper_bounds{obj.totalLoopCtr + 1} = sqrt(obj.muTrainUB);
        end

        % Look for the parameter with the greatest gap / error
        obj.muPrevIdx          = obj.muCurIdx;
        [maxGap, obj.muCurIdx] = max((obj.muTrainUB - obj.muTrainLB) ./ obj.muTrainUB);
        if ~isempty(maxGap)
          obj.errors(end + 1) = maxGap;
        end

        % first check whether we are done here, check if the given tolerance
        % is obtained
        if maxGap < tol
          exflag = 1;
          isDone = true;
          obj.L.info('SCM', sprintf(' stopped greedy loop, max gap %f is smaller than tolerance %f!', maxGap, tol));
          break;
        end

        % looks like we're not done yet, so let's solve the eigenvalue
        % problem. assemble the system for the chosen parameter and calculate
        % the smallest eigenvalue and the corresponding eigenvector
        M = obj.assembleAffineT(obj.offMuTrainMapped(:, obj.muCurIdx));
        [mi, mivec] = obj.computeEV(M, obj.normX, 0);
        % if the eigenvalue mi is positive, then we got the largest eigenvalue
        % and have to do it again with a shift
        if mi >= 0
          [mi, mivec] = obj.computeEV(M, obj.normX, mi);
        end

        % the new eigenvalue is lower and upper bound for this parameter mu
        obj.muTrainLB(obj.muCurIdx) = mi;
        obj.muTrainUB(obj.muCurIdx) = mi;

        % calculate the values of ystar for the current parameter. this is done
        % by iterating over all the addends of the affine decomposition and
        % multiplication with the eigenvector from above from both sides
        ystar = zeros(obj.nQt, 1);
        for idx = 1:obj.nQt
          ystar(idx) = mivec.' * obj.affineT{idx} * mivec;
        end

        % save the newly generated data for the online stage
        obj.offMuCk(:, end + 1)       = obj.offMuTrain(:, obj.muCurIdx);
        obj.offMuCkMapped(:, end + 1) = obj.offMuTrainMapped(:, obj.muCurIdx);
        obj.offMuCkIndex(end + 1)     = obj.muCurIdx;
        obj.offYstarCk(:, end + 1)    = ystar;
        obj.offAlphaCk(end + 1)       = mi;
        obj.offMuChosen(obj.muCurIdx) = true;

        % logging
        obj.L.info('SCM', sprintf(' cycle %7.d of max %7.d done. max gap %g', ...
                   curLoopCtr, maxIter, maxGap));

        % update the current and total loop count
        curLoopCtr       = curLoopCtr + 1;
        obj.totalLoopCtr = obj.totalLoopCtr + 1;

        % and now we check some more breaking conditions
        if curLoopCtr > maxIter
          % we stop because we performed the maxIter number of greedy cycles.
          % it's now maybe a good idea to check maxGap and decide whether the
          % offline stage should be restarted with a larger maxIter
          exflag = 2;
          isDone = true;
          obj.L.info('SCM', sprintf(' stopped greedy loop, number of max iterations %d reached!', maxIter));
          break;
        elseif size(obj.offMuCk, 2) == obj.nMuTrain
          % we stop because the training set is empty. this really shouldn't
          % happen if we want good bounds...
          exflag = 3;
          isDone = true;
          obj.L.info('SCM', ' stopped greedy loop, training set is already empty!');
          break;
        end
        toc(scmtic)
      end
    end

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
      for qdx = 1:obj.solver.nQb
        newparam(idx) = param(qdx) * (param(qdx) - sum(param(1:end ~= qdx)));
        idx = idx + 1;
      end

      % and now the simpler second kind
      for qdx = 1:obj.solver.nQb
        for pdx = (qdx + 1):obj.solver.nQb
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

      % sum the addends of the affine representation
      T = sparse(obj.solver.nTrialDim, obj.solver.nTrialDim);
      for idx = 1:obj.nQt
        T = T + param(idx) * obj.affineT{idx};
      end
    end

    function [mi, ma, retries] = computeMinMaxEV(obj, A, B)
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
      [lm, ~, retries1] = obj.computeEV(A, B, 0);
      % if it's positive, then it's the maximal eigenvalue, else the
      % minimal eigenvalue
      if lm >= 0
        ma = lm;
        % now shift to get the minimal eigenvalue
        [mi, ~, retries2] = obj.computeEV(A, B, lm);
      else
        mi = lm;
        % now shift to get the maximal eigenvalue
        [ma, ~, retries2] = obj.computeEV(A, B, lm);
      end

      retries = retries1 + retries2;
    end

    function [lm, lmvec, retries] = computeEV(~, A, B, shift)
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

      retries = 0;

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
        if strcmp(err.identifier, 'MATLAB:eigs:NoEigsConverged') || strcmp(err.identifier, 'MATLAB:eigs:ARPACKroutineErrorMinus14')
          % flag to check whether we found a eigenvalue
          hasEV = false;
          % decrease the accuracy
          opts.tol = opts.tol * 100;
          % set the number of used Lanczos vectors
          opts.p = 50;
          while ~hasEV
            try
              retries = retries + 1;
              % and now try again
              [lmvec, lm] = eigs(A - shift * B, B, 1, 'lm', opts);
              lm = lm  + shift;
              hasEV = true;
            catch err2
              hasEV = false;
              if strcmp(err2.identifier, 'MATLAB:eigs:NoEigsConverged') || strcmp(err2.identifier, 'MATLAB:eigs:ARPACKroutineErrorMinus14')
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
      % Prepare the SCM object for the offline stage.
      %
      % This mainly consists of the creation of the addends for the affine
      % representation of the Y-norm of the supremizers.

      % Create local copies of the truth solver structures.
      obj.normX   = obj.solver.TrNorm;
      % more truth solver structures, this time only needed for the preparation
      normY   = obj.solver.TeNorm;
      affineB = obj.solver.Lhs;

      % set the total number of needed addends for the affine representation
      obj.nQt = obj.solver.nQb * (obj.solver.nQb + 1) / 2;

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
      for qdx = 1:obj.solver.nQb
        obj.affineT{idx} = affineB{qdx}.' * (normY \ affineB{qdx});
        idx = idx + 1;
      end

      % and now assemble the affine addends of the second kind
      for qdx = 1:obj.solver.nQb
        for pdx = (qdx + 1):obj.solver.nQb
          obj.affineT{idx} = (affineB{qdx} + affineB{pdx}).' * ...
            (normY \ (affineB{qdx} + affineB{pdx}));
          idx = idx + 1;
        end
      end

      % reset the greedy loop counter
      obj.totalLoopCtr = 0;
    end
  end
end


