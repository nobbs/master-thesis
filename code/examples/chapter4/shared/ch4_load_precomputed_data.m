% handle the tikz generation toggle
if ~exist('generate_tikz', 'var')
    generate_tikz = false;
end

%% example 1
load('ch4ex1_precomputed.mat');

%% and now, let's plot that stuff
% first figure is the exact inf-sup-constant `beta_{\mathcal N}`
figure();
semilogy(tvalues, exact, '-o')
ylim([5e-4, 1]);
xlim([00, 310]);
legend('5', '10', '15', '20', '25', 'Location', 'southeast');
legend('boxoff');
xlabel('Anzahl Zeitgitterpunkte');
drawnow

if generate_tikz
    matlab2tikz('stability_sine_dataset1_fig_1.tikz');
end

figure();
semilogy(tvalues, bounds, '-o')
ylim([5e-4, 1]);
xlim([00, 310]);
legend('5', '10', '15', '20', '25', 'Location', 'southeast');
legend('boxoff');
xlabel('Anzahl Zeitgitterpunkte');
drawnow

if generate_tikz
    matlab2tikz('stability_sine_dataset1_fig_2.tikz');
end

% clear;

%% example 2
load('ch4ex2_precomputed.mat');

%% and now, let's plot that stuff
% first figure is the exact inf-sup-constant `beta_{\mathcal N}`
figure();
semilogy(tvalues, exact, '-o')
ylim([5e-4, 1]);
xlim([00, 310]);
legend('5', '10', '15', '20', '25', 'Location', 'southeast');
legend('boxoff');
xlabel('Anzahl Zeitgitterpunkte');
drawnow

if generate_tikz
    matlab2tikz('stability_sine_dataset2_fig_1.tikz');
end

figure();
semilogy(tvalues, bounds, '-o')
ylim([5e-4, 1]);
xlim([00, 310]);
legend('5', '10', '15', '20', '25', 'Location', 'southeast');
legend('boxoff');
xlabel('Anzahl Zeitgitterpunkte');
drawnow

if generate_tikz
    matlab2tikz('stability_sine_dataset2_fig_2.tikz');
end

% clear;

%% example 3
load('ch4ex3_precomputed.mat');

%% and now, let's plot that stuff
% first figure is the exact inf-sup-constant `beta_{\mathcal N}`
figure();
semilogy(tvalues, exact, '-o')
ylim([5e-4, 1]);
xlim([00, 310]);
legend('5', '10', '15', '20', '25', 'Location', 'southeast');
legend('boxoff');
xlabel('Anzahl Zeitgitterpunkte');
drawnow

if generate_tikz
    matlab2tikz('stability_fourier_dataset1_fig_1.tikz');
end

figure();
semilogy(tvalues, bounds, '-o')
ylim([5e-4, 1]);
xlim([00, 310]);
legend('5', '10', '15', '20', '25', 'Location', 'southeast');
legend('boxoff');
xlabel('Anzahl Zeitgitterpunkte');
drawnow

if generate_tikz
    matlab2tikz('stability_fourier_dataset1_fig_2.tikz');
end

% clear;
