
p_tstats = nan(1,3);

% part 1:
S1 = plot_results(process_group, 'g1');
% p_tstats(1) = S1.stats.t{2};

% part 2:
S2 = plot_results(process_group_v3, 'g2');
% p_tstats(2) = S2.stats.t{2};

% part 3:
S3 = plot_results(process_group_v4, 'g3');
% p_tstats(3) = S3.stats.t{2};