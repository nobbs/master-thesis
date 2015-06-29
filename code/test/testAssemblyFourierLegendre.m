%% same number of temporal and spatial basis functions
% test the assembly of the stiffness matrix for chosen numbers of spatial and
% temporal basis functions resulting in a quadratic system of linear equations
X = [2, 5];

assembly = AssemblyFourierLegendre();
assembly.tspan = [0 1];
assembly.xspan = [0 1];

assembly.coeffLaplacian = 0.1;
assembly.coeffOffset  = 2;

% without normalization
for idx = 1:length(X)
  assembly.setNumberOfAnsatzFuncs(X(idx), X(idx));
  assembly.setNumberOfTestFuncsFromAnsatzFuncs();

  assembly.precomputeNormalization();
  [LRef, ~]  = assembly.assembleFieldIndependentMatrixSlow();
  [LTest, ~] = assembly.assembleFieldIndependentMatrix();

  assert(max(max(abs(LTest - LRef))) < 1e-8);
end

% with normalization
assembly.useNormalization = true;
for idx = 1:length(X)
  assembly.setNumberOfAnsatzFuncs(X(idx), X(idx));
  assembly.setNumberOfTestFuncsFromAnsatzFuncs();

  assembly.precomputeNormalization();
  [LRef, ~]  = assembly.assembleFieldIndependentMatrixSlow();
  [LTest, ~] = assembly.assembleFieldIndependentMatrix();

  assert(max(max(abs(LTest - LRef))) < 1e-8);
end


%% non-standard tspan
% test the assembly of the stiffness matrix for a temporal interval different
% from [0 1]
X = 5;
span1 = [0 1/3];
span2 = [1/3 2];

assembly = AssemblyFourierLegendre();
assembly.xspan = [0 1];
assembly.coeffLaplacian = 0.1;
assembly.coeffOffset  = 2;
assembly.setNumberOfAnsatzFuncs(X, X);
assembly.setNumberOfTestFuncsFromAnsatzFuncs();

% without normalization
assembly.tspan = span1;
assembly.precomputeNormalization();
[LRef, ~]  = assembly.assembleFieldIndependentMatrixSlow();
[LTest, ~] = assembly.assembleFieldIndependentMatrix();
assert(max(max(abs(LTest - LRef))) < 1e-8);

assembly.tspan = span2;
assembly.precomputeNormalization();
[LRef, ~]  = assembly.assembleFieldIndependentMatrixSlow();
[LTest, ~] = assembly.assembleFieldIndependentMatrix();
assert(max(max(abs(LTest - LRef))) < 1e-8);

% with normalization
assembly.useNormalization = true;
assembly.tspan = span1;
assembly.precomputeNormalization();
[LRef, ~]  = assembly.assembleFieldIndependentMatrixSlow();
[LTest, ~] = assembly.assembleFieldIndependentMatrix();
assert(max(max(abs(LTest - LRef))) < 1e-8);

assembly.tspan = span2;
assembly.precomputeNormalization();
[LRef, ~]  = assembly.assembleFieldIndependentMatrixSlow();
[LTest, ~] = assembly.assembleFieldIndependentMatrix();
assert(max(max(abs(LTest - LRef))) < 1e-8);


%% non-standard xspan
% test the assembly of the stiffness matrix for a spatial interval different
% from [0 1]
X = 5;
span1 = [0 10];

assembly = AssemblyFourierLegendre();
assembly.tspan = [0 1];
assembly.coeffLaplacian = 0.1;
assembly.coeffOffset  = 2;
assembly.setNumberOfAnsatzFuncs(X, X);
assembly.setNumberOfTestFuncsFromAnsatzFuncs();

% without normalization
assembly.xspan = span1;
assembly.precomputeNormalization();
[LRef, ~]  = assembly.assembleFieldIndependentMatrixSlow();
[LTest, ~] = assembly.assembleFieldIndependentMatrix();
assert(max(max(abs(LTest - LRef))) < 1e-8);

% with normalization
assembly.useNormalization = true;
assembly.xspan = span1;
assembly.precomputeNormalization();
[LRef, ~]  = assembly.assembleFieldIndependentMatrixSlow();
[LTest, ~] = assembly.assembleFieldIndependentMatrix();
assert(max(max(abs(LTest - LRef))) < 1e-8);


%% non-quadratic system
% test whether the assembly works correctly for the cases where a non-quadratic
% system of linear equation is constructed.
XN = [2, 4];
XM = [3, 5];

YN = [3, 5];
YM = [4, 7];
YO = [2, 6];

assembly = AssemblyFourierLegendre();

assembly.tspan = [0 1];
assembly.xspan = [0 1];

assembly.coeffLaplacian = 0.1;
assembly.coeffOffset  = 2;

% without normalization
for idx = 1:length(XN)
  assembly.setNumberOfAnsatzFuncs(XN(idx), XM(idx));
  assembly.setNumberOfTestFuncs(YN(idx), YM(idx), YO(idx));

  assembly.precomputeNormalization();
  [LRef, ~]  = assembly.assembleFieldIndependentMatrixSlow();
  [LTest, ~] = assembly.assembleFieldIndependentMatrix();

  assert(max(max(abs(LTest - LRef))) < 1e-8);
end

% with normalization
assembly.useNormalization = true;
for idx = 1:length(XN)
  assembly.setNumberOfAnsatzFuncs(XN(idx), XM(idx));
  assembly.setNumberOfTestFuncs(YN(idx), YM(idx), YO(idx));

  assembly.precomputeNormalization();
  [LRef, ~]  = assembly.assembleFieldIndependentMatrixSlow();
  [LTest, ~] = assembly.assembleFieldIndependentMatrix();

  assert(max(max(abs(LTest - LRef))) < 1e-8);
end


%% Omega from sine series expansion
% check if the field-dependent matrices are correct for the sine series
% expansion case
XN = [2, 4];
XM = [3, 5];

YN = [3, 5];
YM = [4, 7];
YO = [2, 6];

N = 15;

assembly = AssemblyFourierLegendre();

assembly.tspan = [1/2 5/2];
assembly.xspan = [0 10];

assembly.coeffLaplacian = 0.1;
assembly.coeffOffset  = 2;

% without normalization
for idx = 1:length(XN)
  assembly.setNumberOfAnsatzFuncs(XN(idx), XM(idx));
  assembly.setNumberOfTestFuncs(YN(idx), YM(idx), YO(idx));

  assembly.precomputeNormalization();
  ORef = assembly.assembleFieldDependentMatrixForSineSeriesSlow(N);
  OTest = assembly.assembleFieldDependentMatrixForSineSeries(N);

  for jdx = 1:N
    assert(max(max(abs(OTest{jdx} - ORef{jdx}))) < 1e-8);
  end
end

% with normalization
assembly.useNormalization = true;
for idx = 1:length(XN)
  assembly.setNumberOfAnsatzFuncs(XN(idx), XM(idx));
  assembly.setNumberOfTestFuncs(YN(idx), YM(idx), YO(idx));

  assembly.precomputeNormalization();
  ORef = assembly.assembleFieldDependentMatrixForSineSeriesSlow(N);
  OTest = assembly.assembleFieldDependentMatrixForSineSeries(N);

  for jdx = 1:N
    assert(max(max(abs(OTest{jdx} - ORef{jdx}))) < 1e-8);
  end
end


%% Omega from fourier series expansion
% check if the field-dependent matrices are correct for the fourier series
% expansion case
XN = [2, 4];
XM = [3, 5];

YN = [3, 5];
YM = [4, 7];
YO = [2, 6];

N = 15;

assembly = AssemblyFourierLegendre();

assembly.tspan = [1/2 5/2];
assembly.xspan = [0 10];

assembly.coeffLaplacian = 0.1;
assembly.coeffOffset  = 2;


% without normalization
for idx = 1:length(XN)
  assembly.setNumberOfAnsatzFuncs(XN(idx), XM(idx));
  assembly.setNumberOfTestFuncs(YN(idx), YM(idx), YO(idx));

  assembly.precomputeNormalization();
  ORef = assembly.assembleFieldDependentMatrixForFourierSeriesSlow(N);
  OTest = assembly.assembleFieldDependentMatrixForFourierSeries(N);

  for jdx = 1:N
    assert(max(max(abs(OTest{jdx} - ORef{jdx}))) < 1e-8);
  end
end

% with normalization
assembly.useNormalization = true;
for idx = 1:length(XN)
  assembly.setNumberOfAnsatzFuncs(XN(idx), XM(idx));
  assembly.setNumberOfTestFuncs(YN(idx), YM(idx), YO(idx));

  assembly.precomputeNormalization();
  ORef = assembly.assembleFieldDependentMatrixForFourierSeriesSlow(N);
  OTest = assembly.assembleFieldDependentMatrixForFourierSeries(N);

  for jdx = 1:N
    assert(max(max(abs(OTest{jdx} - ORef{jdx}))) < 1e-8);
  end
end



%% ansatz and test norm matrix assembly
% test the assembly of the discrete ansatz and test subspace norms by comparison
% with the slow numerical quadrature assembly methods
XN = [2, 4];
XM = [3, 5];

YN = [3, 5];
YM = [4, 7];
YO = [2, 6];

assembly = AssemblyFourierLegendre();

assembly.tspan = [0 1];
assembly.xspan = [0 1];

assembly.coeffLaplacian = 0.1;
assembly.coeffOffset  = 2;

for idx = 1:length(XN)
  assembly.setNumberOfAnsatzFuncs(XN(idx), XM(idx));
  assembly.setNumberOfTestFuncs(YN(idx), YM(idx), YO(idx));

  assembly.precomputeNormalization();

  MARef  = assembly.assembleAnsatzNormMatrixSlow();
  MATest = assembly.assembleAnsatzNormMatrix();

  MTRef  = assembly.assembleTestNormMatrixSlow();
  MTTest = assembly.assembleTestNormMatrix();

  assert(max(max(abs(MATest - MARef))) < 1e-8);
  assert(max(max(abs(MTTest - MTRef))) < 1e-8);
end


%% Quick test without normalization
% Assemble the whole system of linear equations, solve it in time till an
% intermediate endpoint, set the value of this solution as the initial condition
% of a second run and compare that with a solution from start to end.

X = 15;
Y = 15;
N = 10;
rcoeffs = randn(N, 1);
breakpoint = rand(1);
gt1 = linspace(0, breakpoint);
gt1 = gt1(1:end-1);
gt2 = linspace(breakpoint, 1);
gtw = [gt1, gt2];
[mt1, mx1] = meshgrid(gt1, linspace(0, 5));
[mt2, mx2] = meshgrid(gt2, linspace(0, 5));
mtw = [mt1, mt2];
mxw = [mx1, mx2];

assembly = AssemblyFourierLegendre();
assembly.coeffLaplacian = 0.1;
assembly.coeffOffset = 0;
assembly.setNumberOfAnsatzFuncs(X, Y);
assembly.setNumberOfTestFuncsFromAnsatzFuncs();

% compute solution for the whole temporal interval
assembly.tspan = [0, 1];
assembly.xspan = [0, 5];
assembly.precomputeNormalization();
[LhsFI, ~] = assembly.assembleFieldIndependentMatrix();
LhsFD = assembly.assembleFieldDependentMatrixForFourierSeries(N);

Lhs = LhsFI;
for idx = 1:N
  Lhs = Lhs + rcoeffs(idx) * LhsFD{idx};
end

RhsOnes = assembly.assembleVectorOnes();
solComplete = Lhs \ RhsOnes;
solCompleteEval = assembly.solutionFuncFromCoeffs(solComplete, mtw, mxw);

% compute solution for the first part of the temporal interval
assembly.tspan = [0, breakpoint];
assembly.precomputeNormalization();
[LhsFIOne, ~] = assembly.assembleFieldIndependentMatrix();
LhsFDOne = assembly.assembleFieldDependentMatrixForFourierSeries(N);
LhsOne = LhsFIOne;
for idx = 1:N
  LhsOne = LhsOne + rcoeffs(idx) * LhsFDOne{idx};
end
RhsOnes = assembly.assembleVectorOnes();
solPartOne = LhsOne \ RhsOnes;
solPartOneEval = assembly.solutionFuncFromCoeffs(solPartOne, mt1, mx1);

% compute solution for the second part of the temporal interval
assembly.tspan = [breakpoint, 1];
assembly.precomputeNormalization();
[LhsFITwo, ~] = assembly.assembleFieldIndependentMatrix();
LhsFDTwo = assembly.assembleFieldDependentMatrixForFourierSeries(N);
LhsTwo = LhsFITwo;
for idx = 1:N
  LhsTwo = LhsTwo + rcoeffs(idx) * LhsFDTwo{idx};
end
RhsMid = assembly.assembleVectorFromSolutionCoeffs(solPartOne);
solPartTwo = LhsTwo \ RhsMid;
solPartTwoEval = assembly.solutionFuncFromCoeffs(solPartTwo, mt2, mx2);

% @todo remove visual debugging stuff
% figure(1)
% mesh(mtw, mxw, [solPartOneEval, solPartTwoEval]);
% figure(2)
% mesh(mtw, mxw, solCompleteEval);
% max(max(abs([solPartOneEval, solPartTwoEval] - solCompleteEval)))

assert(max(max(abs([solPartOneEval, solPartTwoEval] - solCompleteEval))) < 1e-4);


%% Quick test with normalization
% Assemble the whole system of linear equations, solve it in time till an
% intermediate endpoint, set the value of this solution as the initial condition
% of a second run and compare that with a solution from start to end.

X = 15;
Y = 15;
N = 10;
rcoeffs = randn(N, 1);
breakpoint = rand(1);
gt1 = linspace(0, breakpoint);
gt1 = gt1(1:end-1);
gt2 = linspace(breakpoint, 1);
gtw = [gt1, gt2];
[mt1, mx1] = meshgrid(gt1, linspace(0, 5));
[mt2, mx2] = meshgrid(gt2, linspace(0, 5));
mtw = [mt1, mt2];
mxw = [mx1, mx2];

assembly = AssemblyFourierLegendre();
assembly.coeffLaplacian = 0.1;
assembly.coeffOffset = 0;
assembly.setNumberOfAnsatzFuncs(X, Y);
assembly.setNumberOfTestFuncsFromAnsatzFuncs();
assembly.useNormalization = true;

% compute solution for the whole temporal interval
assembly.tspan = [0, 1];
assembly.xspan = [0, 5];
assembly.precomputeNormalization();
[LhsFI, ~] = assembly.assembleFieldIndependentMatrix();
LhsFD = assembly.assembleFieldDependentMatrixForFourierSeries(N);

Lhs = LhsFI;
for idx = 1:N
  Lhs = Lhs + rcoeffs(idx) * LhsFD{idx};
end

RhsOnes = assembly.assembleVectorOnes();
solComplete = Lhs \ RhsOnes;
solCompleteEval = assembly.solutionFuncFromCoeffs(solComplete, mtw, mxw);

% compute solution for the first part of the temporal interval
assembly.tspan = [0, breakpoint];
assembly.precomputeNormalization();
[LhsFIOne, ~] = assembly.assembleFieldIndependentMatrix();
LhsFDOne = assembly.assembleFieldDependentMatrixForFourierSeries(N);
LhsOne = LhsFIOne;
for idx = 1:N
  LhsOne = LhsOne + rcoeffs(idx) * LhsFDOne{idx};
end
RhsOnes = assembly.assembleVectorOnes();
solPartOne = LhsOne \ RhsOnes;
solPartOneEval = assembly.solutionFuncFromCoeffs(solPartOne, mt1, mx1);

TmpMatrix = assembly.AnsatzNormDiag;

% compute solution for the second part of the temporal interval
assembly.tspan = [breakpoint, 1];
assembly.precomputeNormalization();
[LhsFITwo, ~] = assembly.assembleFieldIndependentMatrix();
LhsFDTwo = assembly.assembleFieldDependentMatrixForFourierSeries(N);
LhsTwo = LhsFITwo;
for idx = 1:N
  LhsTwo = LhsTwo + rcoeffs(idx) * LhsFDTwo{idx};
end
RhsMid = assembly.assembleVectorFromSolutionCoeffs(solPartOne ./ sqrt(diag(TmpMatrix)));
solPartTwo = LhsTwo \ RhsMid;
solPartTwoEval = assembly.solutionFuncFromCoeffs(solPartTwo, mt2, mx2);

% @todo remove visual debugging stuff
% figure(1)
% mesh(mtw, mxw, [solPartOneEval, solPartTwoEval]);
% figure(2)
% mesh(mtw, mxw, solCompleteEval);
% max(max(abs([solPartOneEval, solPartTwoEval] - solCompleteEval)))

assert(max(max(abs([solPartOneEval, solPartTwoEval] - solCompleteEval))) < 1e-4);
