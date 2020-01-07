
p_tstats = nan(1,3);

% part 1:
h1 = process_group;
S1 = plot_results(h1, 'g1');
% p_tstats(1) = S1.stats.t{2};

% part 2:
h2 = process_group_v3;
S2 = plot_results(h2, 'g2');
% p_tstats(2) = S2.stats.t{2};

% part 3:
h3 = process_group_v4;
S3 = plot_results(h3, 'g3');
% p_tstats(3) = S3.stats.t{2};

% Experiment E6:
h6 = process_group_E6;
S6 = plot_results(h6, 'g4');

% part 1:
% S1 = plot_results_learning(process_group, 'g1');
% p_tstats(1) = S1.stats.t{2};

% part 2:
% S2 = plot_results_learning(process_group_v3, 'g2');
% p_tstats(2) = S2.stats.t{2};

% part 3:
% S3 = plot_results_learning(process_group_v4, 'g3');
% p_tstats(3) = S3.stats.t{2};