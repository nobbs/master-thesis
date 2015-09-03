% This file loads the precomputed self consistent field theory example data and
% creates the plots as seen in Chapter 1.

% handle the tikz generation toggle
if ~exist('generate_tikz', 'var')
    generate_tikz = false;
end

% load the data
load('ch1_scft_example_data.mat');
% load('ch1_scft_example_data2.mat');

%% first, plot the fields
figure();
plot(x.grid, omegaA, x.grid, omegaB);

if generate_tikz
    matlab2tikz('scft_example_fields.tikz');
end

%% fourier coefficients
% set the number of coefficients to plot
numCoeffs = 64;

% calculate sine / cosine coefficients
coeffsA = Utility.valuesToSineCosine(omegaA);
coeffsB = Utility.valuesToSineCosine(omegaB);

% second, plot coefficients
grid = 1:numCoeffs;
= figure();
semilogy(grid, abs(coeffsA(grid)), 'x', grid, abs(coeffsB(grid)), 'o');

if generate_tikz
    matlab2tikz('scft_example_fourier_coeffs.tikz');
end
