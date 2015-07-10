problem = ProblemData;
problem.laplacian = 1;
problem.offset    = 0;
problem.xspan     = [0 1];
problem.xgrid     = linspace(0, 1);
problem.tspan     = [0 1];
problem.tgrid     = linspace(0, 1);
problem.f         = [];
problem.nC        = 1;
problem.seriesIdx = [1];
problem.nF        = 1;

problem
