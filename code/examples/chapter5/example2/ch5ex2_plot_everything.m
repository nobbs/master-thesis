% handle the tikz generation toggle
if ~exist('generate_tikz', 'var')
    generate_tikz = false;
end

%% data source: you can execute the following line. as this could take some time
% (several hours+), it's recommended to load the already precomputed data
% ch5ex2_run;
load('ch5ex2_precomputed.mat');

% As this example has a four dimensional parameter space the only useful plot is
% the maximal rb error plotted against the dimension of the reduced basis.

% create the plot
figure();
semilogy(1:rbm.nTrialRB, rbm.errors, '-*');
xlabel('Anzahl N der Parameter-Samples');
ylabel('Maximum des a posteriori-Fehlerschätzers');
xlim([1, rbm.nTrialRB]);
drawnow;

if generate_tikz
    matlab2tikz(['ch5ex2_rbm_error.tikz']);
end

%% And now: table generation!
% rbm error
table1 = [1:length(rbm.errors); rbm.errors].';
print_table(table1, {'%d', '%e'}, {'N', 'err'}, 'printMode', 'latex')
