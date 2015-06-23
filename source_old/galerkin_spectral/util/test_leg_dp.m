LP = {};
LP{1} = @(x) ones(size(x, 1), size(x, 2));
LP{2} = @(x) 2 * x - 1;
LP{3} = @(x) 6*x.^2 - 6 * x + 1;
LP{4} = @(x) 20*x.^3 -30*x.^2 + 12 * x -1;
LP{5} = @(x) 70*x.^4-140*x.^3 + 90*x.^2 - 20 *x +1;

dLP = {};
dLP{1} = @(x) zeros(size(x, 1), size(x, 2));
dLP{2} = @(x) 2 * ones(size(x, 1), size(x, 2));
dLP{3} = @(x) 12*x - 6;
dLP{4} = @(x) 60*x.^2 -60*x + 12;
dLP{5} = @(x) 280*x.^3-420*x.^2 + 180*x - 20;
