% This file loads the precomputed self consistent field theory example data and
% creates the plots as seen in Chapter 1.

% load the data
load('ch1_scft_example_data.mat');

% first, plot the fields
f1 = figure(1);
plot(x.grid, omegaA, x.grid, omegaB);

% set the number of coefficients to plot
numCoeffs = 64;

% calculate sine / cosine coefficients
coeffsA = Utility.valuesToSineCosine(omegaA);
coeffsB = Utility.valuesToSineCosine(omegaB);

% second, plot coefficients
grid = 1:numCoeffs;
f2 = figure(2);
semilogy(grid, abs(coeffsA(grid)), 'x', grid, abs(coeffsB(grid)), 'o');

% and now convert figures to tikz code
% set(f1, 'Position', [300, 300, 400, 300])
% set(f2, 'Position', [300, 300, 400, 300])
% figure(1);
% matlab2tikz('scft_example_fields.tikz');
% figure(2);
% matlab2tikz('scft_example_fourier_coeffs.tikz');
