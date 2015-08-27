% example 1
tic;
ch4ex1_run;
toc
save('ch4ex1_precomputed.mat', 'cfls', 'bounds', 'exact', 'tvalues', 'svalues');
clear

% example 2
tic;
ch4ex2_run;
toc
save('ch4ex2_precomputed.mat', 'cfls', 'bounds', 'exact', 'tvalues', 'svalues');
clear

% example 4
tic;
ch4ex3_run;
toc
save('ch4ex3_precomputed.mat', 'cfls', 'bounds', 'exact', 'tvalues', 'svalues');
clear
