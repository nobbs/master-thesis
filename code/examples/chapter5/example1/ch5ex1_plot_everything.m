% handle the tikz generation toggle
if ~exist('generate_tikz', 'var')
  generate_tikz = false;
end

%% data source: you can execute the following line. as this could take some time
% (several minutes+), it's recommended to load the already precomputed data
% ch5ex1_run; ch5ex1_run_analysis;
load('ch5ex1_precomputed.mat');

%% First: Plot data from the scm greedy construction that shows how the error
%gets smaller with each iteration
for idx = 1:length(scm.scm_lower_bounds)
  figure();
  plot(ptrain, ex_values.', '-', ...
    ptrain, scm.scm_lower_bounds{idx}, '--', ...
    ptrain, scm.scm_upper_bounds{idx}, '-.', ...
    scm.offMuCk(1:idx), sqrt(scm.offAlphaCk(1:idx)), '*');
  drawnow;

  if generate_tikz
    ylim([0.0920, 0.1020]);
    drawnow;
    matlab2tikz(['ch5ex1_scm_plot_', num2str(idx), '.tikz']);
  end
end

%% Second: Same for the greedy stage of the rbm.
% plotting
for idx = 1:length(rbm.rbm_exact_error)
  figure();
  semilogy(ptrain, rbm.rbm_exact_error{idx}, ...
    ptrain, rbm.rbm_error_bounds{idx});
  drawnow;

  if generate_tikz
    ylim([1e-15, 1]);
    drawnow;
    matlab2tikz(['ch5ex1_rbmerr_plot_', num2str(idx), '.tikz']);
  end
end

%% Third: Plot the effectivity
eff = {};
for idx = 1:length(rbm.rbm_exact_error)
  figure();
  eff{idx} = rbm.rbm_error_bounds{idx} ./ rbm.rbm_exact_error{idx}.';
  semilogy(ptrain, eff{idx});
  drawnow;

  if generate_tikz
    matlab2tikz(['ch5ex1_eff_plot_', num2str(idx), '.tikz']);
  end
end

%% And now: table generation!
% scm error
table1 = [1:length(scm.errors); scm.errors].';
print_table(table1, {'%d', '%e'}, {'N', 'err'}, 'printMode', 'latex')

% rbm error // average, median, max effectivity
ave_eff = zeros(1, length(rbm.errors));
med_eff = zeros(1, length(rbm.errors));
max_eff = zeros(1, length(rbm.errors));
for idx = 1:length(rbm.errors)
  ave_eff(idx) = mean(eff{idx});
  med_eff(idx) = median(eff{idx});
  max_eff(idx) = max(eff{idx});
end

table2 = [(1:length(rbm.errors)); rbm.errors; ave_eff; med_eff; max_eff].';
print_table(table2, {'%d', '%e', '%e', '%e', '%e'}, {'N', 'err', 'average eff', 'median eff', 'max eff'}, 'printMode', 'latex')

% runtimes
disp(sprintf('scm runtime: %f', runTimes.scm));
disp(sprintf('rbm runtime: %f', runTimes.rbm));
