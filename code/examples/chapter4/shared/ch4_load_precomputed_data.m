%% example 1
load('ch4ex1_precomputed.mat');

%% and now, let's plot that stuff
% first figure is the exact inf-sup-constant `beta_{\mathcal N}`
figure(1);
semilogy(tvalues, exact, '-o')
ylim([5e-4, 1]);
xlim([00, 310]);
legend('5', '10', '15', '20', '25', 'Location', 'southeast');
legend('boxoff');
xlabel('Anzahl Zeitgitterpunkte');

drawnow
matlab2tikz('stability_sine_dataset1_fig_1.tikz');

figure(2);
semilogy(tvalues, bounds, '-o')
ylim([5e-4, 1]);
xlim([00, 310]);
legend('5', '10', '15', '20', '25', 'Location', 'southeast');
legend('boxoff');
xlabel('Anzahl Zeitgitterpunkte');

drawnow
matlab2tikz('stability_sine_dataset1_fig_2.tikz');
clear;

%% example 2
load('ch4ex2_precomputed.mat');

%% and now, let's plot that stuff
% first figure is the exact inf-sup-constant `beta_{\mathcal N}`
figure(3);
semilogy(tvalues, exact, '-o')
ylim([5e-4, 1]);
xlim([00, 310]);
legend('5', '10', '15', '20', '25', 'Location', 'southeast');
legend('boxoff');
xlabel('Anzahl Zeitgitterpunkte');

drawnow
matlab2tikz('stability_sine_dataset2_fig_1.tikz');

figure(4);
semilogy(tvalues, bounds, '-o')
ylim([5e-4, 1]);
xlim([00, 310]);
legend('5', '10', '15', '20', '25', 'Location', 'southeast');
legend('boxoff');
xlabel('Anzahl Zeitgitterpunkte');

drawnow
matlab2tikz('stability_sine_dataset2_fig_2.tikz');
clear;

%% example 3
load('ch4ex3_precomputed.mat');

%% and now, let's plot that stuff
% first figure is the exact inf-sup-constant `beta_{\mathcal N}`
figure(5);
semilogy(tvalues, exact, '-o')
ylim([5e-4, 1]);
xlim([00, 310]);
legend('5', '10', '15', '20', '25', 'Location', 'southeast');
legend('boxoff');
xlabel('Anzahl Zeitgitterpunkte');

drawnow
matlab2tikz('stability_fourier_dataset1_fig_1.tikz');

figure(6);
semilogy(tvalues, bounds, '-o')
ylim([5e-4, 1]);
xlim([00, 310]);
legend('5', '10', '15', '20', '25', 'Location', 'southeast');
legend('boxoff');
xlabel('Anzahl Zeitgitterpunkte');

drawnow
matlab2tikz('stability_fourier_dataset1_fig_2.tikz');
clear;
