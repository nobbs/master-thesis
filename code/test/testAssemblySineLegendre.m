%% same number of temporal and spatial basis functions
% test the assembly of the stiffness matrix for chosen numbers of spatial and
% temporal basis functions resulting in a quadratic system of linear equations
X = [1, 5, 10];

assembly = AssemblySineLegendre();

assembly.tspan = [0 1];
assembly.xspan = [0 1];

assembly.coeffLaplacian = 0.1;
assembly.coeffOffset  = 2;

for idx = 1:length(X)
  assembly.setNumberOfAnsatzFuncs(X(idx), X(idx));
  assembly.setNumberOfTestFuncsFromAnsatzFuncs();

  LRef  = assembly.assembleFieldIndependentMatrixSlow();
  LTest = assembly.assembleFieldIndependentMatrix();

  assert(max(max(abs(LTest - LRef))) < 1e-8);
end


%% non-standard tspan
% test the assembly of the stiffness matrix for a temporal interval different
% from [0 1]
X = 5;
span1 = [0 1/3];
span2 = [1/3 2];

assembly = AssemblySineLegendre();
assembly.xspan = [0 1];
assembly.coeffLaplacian = 0.1;
assembly.coeffOffset  = 2;
assembly.setNumberOfAnsatzFuncs(X, X);
assembly.setNumberOfTestFuncsFromAnsatzFuncs();

assembly.tspan = span1;
LRef  = assembly.assembleFieldIndependentMatrixSlow();
LTest = assembly.assembleFieldIndependentMatrix();
assert(max(max(abs(LTest - LRef))) < 1e-8);

assembly.tspan = span2;
LRef  = assembly.assembleFieldIndependentMatrixSlow();
LTest = assembly.assembleFieldIndependentMatrix();
assert(max(max(abs(LTest - LRef))) < 1e-8);


%% non-standard xspan
% test the assembly of the stiffness matrix for a spatial interval different
% from [0 1]
X = 5;
span1 = [0 10];

assembly = AssemblySineLegendre();
assembly.tspan = [0 1];
assembly.coeffLaplacian = 0.1;
assembly.coeffOffset  = 2;
assembly.setNumberOfAnsatzFuncs(X, X);
assembly.setNumberOfTestFuncsFromAnsatzFuncs();

assembly.xspan = span1;
LRef  = assembly.assembleFieldIndependentMatrixSlow();
LTest = assembly.assembleFieldIndependentMatrix();
assert(max(max(abs(LTest - LRef))) < 1e-8);


%% non-quadratic system
% test whether the assembly works correctly for the cases where a non-quadratic
% system of linear equation is constructed.
XN = [3, 6];
XM = [4, 8];

YN = [5, 10];
YM = [4, 8];
YO = [3, 9];

assembly = AssemblySineLegendre();

assembly.tspan = [0 1];
assembly.xspan = [0 1];

assembly.coeffLaplacian = 0.1;
assembly.coeffOffset  = 2;

for idx = 1:length(XN)
  assembly.setNumberOfAnsatzFuncs(XN(idx), XM(idx));
  assembly.setNumberOfTestFuncs(YN(idx), YM(idx), YO(idx));

  LRef  = assembly.assembleFieldIndependentMatrixSlow();
  LTest = assembly.assembleFieldIndependentMatrix();

  assert(max(max(abs(LTest - LRef))) < 1e-8);
end


%% Omega from sine series expansion
% check if the field-dependent matrices are correct for the sine series
% expansion case
XN = [3, 6];
XM = [4, 8];

YN = [5, 10];
YM = [4, 8];
YO = [3, 9];

N = 15;

assembly = AssemblySineLegendre();

assembly.tspan = [1/2 5/2];
assembly.xspan = [0 10];

assembly.coeffLaplacian = 0.1;
assembly.coeffOffset  = 2;

for idx = 1:length(XN)
  assembly.setNumberOfAnsatzFuncs(XN(idx), XM(idx));
  assembly.setNumberOfTestFuncs(YN(idx), YM(idx), YO(idx));

  ORef = assembly.assembleFieldDependentMatrixForSineSeriesSlow(N);
  OTest = assembly.assembleFieldDependentMatrixForSineSeries(N);

  for jdx = 1:N
    assert(max(max(abs(OTest{jdx} - ORef{jdx}))) < 1e-8);
  end
end


%% Omega from fourier series expansion
% check if the field-dependent matrices are correct for the fourier series
% expansion case
XN = [3, 6];
XM = [4, 8];

YN = [5, 10];
YM = [4, 8];
YO = [3, 9];

N = 15;

assembly = AssemblySineLegendre();

assembly.tspan = [1/2 5/2];
assembly.xspan = [0 10];

assembly.coeffLaplacian = 0.1;
assembly.coeffOffset  = 2;

for idx = 1:length(XN)
  assembly.setNumberOfAnsatzFuncs(XN(idx), XM(idx));
  assembly.setNumberOfTestFuncs(YN(idx), YM(idx), YO(idx));

  ORef = assembly.assembleFieldDependentMatrixForFourierSeriesSlow(N);
  OTest = assembly.assembleFieldDependentMatrixForFourierSeries(N);

  for jdx = 1:N
    assert(max(max(abs(OTest{jdx} - ORef{jdx}))) < 1e-8);
  end
end
