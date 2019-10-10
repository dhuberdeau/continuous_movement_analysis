function stats_results = plot_results(H_results)
% function plot_results(H_results)
%
% Plots results generated from process_group* function

summary_stats_freq_range = [4 6];

H_h_nc_2 = H_results.high_pt.nc;
H_l_nc_2 = H_results.low_pt.nc;
H_h_np_2 = H_results.high_pt.np;
H_l_np_2 = H_results.low_pt.np;
H_h_p_2 = H_results.high_pt.p;
H_l_p_2 = H_results.low_pt.p;
f = H_results.freq;

figure; 
subplot(2,1,1);
hold on
plot(f, nanmean(H_h_nc_2(:,:,end), 2), 'k-', 'Linewidth',2)
plot(f, nanmean(H_l_nc_2(:,:,end), 2), 'r-', 'Linewidth',2)
title('No-cue');
legend('High-pt', 'Low-pt');
subplot(2,1,2); hold on;
plot(f, nanmean(H_h_np_2(:,:,end) - H_l_np_2(:,:,end), 2),...
    'r-', 'Linewidth',2)
plot([0 30], [0 0], '-', 'Color', [.5 .5 .5]);
title('No-cue');
legend('Diff: High - Low');
axis([0 10 -10 4])

figure; 
subplot(2,2,1); hold on
plot(f, nanmean(H_l_nc_2(:,:,end), 2), 'b-', 'Linewidth',2)
plot(f, nanmean(H_h_p_2(:,:,1), 2), 'k.-', 'Linewidth',1)
plot(f, nanmean(H_h_p_2(:,:,end), 2), 'k-', 'Linewidth',2)
plot(f, nanmean(H_l_p_2(:,:,1), 2), 'r.-', 'Linewidth',1)
plot(f, nanmean(H_l_p_2(:,:,end), 2), 'r-', 'Linewidth',2)
title('Predictive')
legend('low-pt / no-cue', 'high-pt / block1', 'high-pt / known',...
    'low-pt / block1', 'low-pt / known')
axis([0 10 -110 -80])

subplot(2,2,2); hold on
plot(f, nanmean(H_l_nc_2(:,:,end), 2), 'b-', 'Linewidth',2)
plot(f, nanmean(H_h_np_2(:,:,1), 2), 'k.-', 'Linewidth',1)
plot(f, nanmean(H_h_np_2(:,:,end), 2), 'k-', 'Linewidth',2)
plot(f, nanmean(H_l_np_2(:,:,1), 2), 'r.-', 'Linewidth',1)
plot(f, nanmean(H_l_np_2(:,:,end), 2), 'r-', 'Linewidth',2)
title('Non-predictive')
legend('low-pt / no-cue', 'high-pt / block1', 'high-pt / known',...
    'low-pt / block1', 'low-pt / known')
axis([0 10 -110 -80])

subplot(2,2,3); hold on;
plot(f, nanmean(H_h_p_2(:,:,end) - H_l_nc_2(:,:,end), 2), 'b-', 'Linewidth',1)
plot(f, nanmean(H_h_p_2(:,:,end) - H_h_p_2(:,:,1), 2), 'k.-', 'Linewidth',1)
plot(f, nanmean(H_h_p_2(:,:,end) - H_l_p_2(:,:,1), 2), 'r.-', 'Linewidth',1)
plot(f, nanmean(H_h_p_2(:,:,end) - H_l_p_2(:,:,end), 2), 'r-', 'Linewidth',2)
title('Diff: predictive')
plot([0 30], [0 0], '-', 'Color', [.5 .5 .5]);
axis([0 10 -10 4])

subplot(2,2,4); hold on;
plot(f, nanmean(H_h_np_2(:,:,end) - H_l_nc_2(:,:,end), 2), 'b-', 'Linewidth',1)
plot(f, nanmean(H_h_np_2(:,:,end) - H_h_np_2(:,:,1), 2), 'k.-', 'Linewidth',1)

plot(f, nanmean(H_h_np_2(:,:,end) - H_l_np_2(:,:,1), 2), 'r.-', 'Linewidth',1)
plot(f, nanmean(H_h_np_2(:,:, end) - H_l_np_2(:,:,end), 2), 'r-', 'Linewidth',2)
title('Diff: Non-predictive')
plot([0 30], [0 0], '-', 'Color', [.5 .5 .5]);
axis([0 10 -10 4])

%% summary statistics:

% find indicies for desired range of frequencies:
[~, k_f_start] = min(abs(summary_stats_freq_range(1) - f));
[~, k_f_end] = min(abs(summary_stats_freq_range(2) - f));

% plot execution error for each trial type across blocks:
if size(H_results.low_pt.nc, 3) > 1
    % first and last block had NC
    figure; hold on
    
    plot([1, 6], reshape(nanmean(nanmean(...
        H_results.low_pt.nc(k_f_start:k_f_end, :,:),1), 2),...
        1,size(H_results.low_pt.nc,3)), 'rs', 'MarkerSize', 14)
    plot((1:size(H_results.low_pt.p,3))+1,...
        reshape(nanmean(nanmean(...
        H_results.low_pt.p(k_f_start:k_f_end, :,:),1), 2),...
        1,size(H_results.low_pt.p,3)), 'gs', 'MarkerSize', 14)
    plot((1:size(H_results.low_pt.np,3))+1,...
        reshape(nanmean(nanmean(...
        H_results.low_pt.np(k_f_start:k_f_end, :,:),1), 2),...
        1,size(H_results.low_pt.np,3)), 'ms', 'MarkerSize', 14)
    axis([1 6 -99 -90])

    figure; hold on
    plot([1, 6], reshape(nanmean(nanmean(...
        H_results.high_pt.nc(k_f_start:k_f_end, :,:),1), 2),...
        1,size(H_results.high_pt.nc,3)), 'rs', 'MarkerSize', 14)
    plot((1:size(H_results.low_pt.p,3))+1, reshape(nanmean(nanmean(...
        H_results.high_pt.p(k_f_start:k_f_end, :,:),1), 2),...
        1,size(H_results.high_pt.p,3)), 'gs', 'MarkerSize', 14)
    plot((1:size(H_results.low_pt.np,3))+1, reshape(nanmean(nanmean(...
        H_results.high_pt.np(k_f_start:k_f_end, :,:),1), 2),...
        1,size(H_results.high_pt.np,3)), 'ms', 'MarkerSize', 14)
    axis([1 6 -99 -90])
else
    % only first did
    figure; hold on
    plot(1, nanmean(nanmean(H_results.low_pt.nc(k_f_start:k_f_end, :),1), 2),...
        'rs', 'MarkerSize', 14)
    plot(2:6, reshape(nanmean(nanmean(...
        H_results.low_pt.p(k_f_start:k_f_end, :,:),1), 2),...
        1,size(H_results.low_pt.p,3)), 'gs', 'MarkerSize', 14)
    plot(2:6, reshape(nanmean(nanmean(...
        H_results.low_pt.np(k_f_start:k_f_end, :,:),1), 2),...
        1,size(H_results.low_pt.np,3)), 'ms', 'MarkerSize', 14)
    axis([1 6 -99 -90])

    figure; hold on
    plot(1, nanmean(nanmean(H_results.high_pt.nc(k_f_start:k_f_end, :),1), 2),...
        'rs', 'MarkerSize', 14)
    plot(2:6, reshape(nanmean(nanmean(...
        H_results.high_pt.p(k_f_start:k_f_end, :,:),1), 2),...
        1,size(H_results.high_pt.p,3)), 'gs', 'MarkerSize', 14)
    plot(2:6, reshape(nanmean(nanmean(...
        H_results.high_pt.np(k_f_start:k_f_end, :,:),1), 2),...
        1,size(H_results.high_pt.np,3)), 'ms', 'MarkerSize', 14)
    axis([1 6 -99 -90])
end

%% Plot kinematic measures:

temp = H_results.kin_all;
hpt = H_results.metric_all(:, 3);
pred = H_results.metric_all(:, 2);
nocue = H_results.metric_all(:, 4);

hpt(isnan(hpt)) = 0;
pred(isnan(pred)) = 0;
nocue(isnan(nocue)) = 0;
hpt = boolean(hpt);
pred = boolean(pred);
nocue = boolean(nocue);

figure;
subplot(411); hold on
t_HPN = temp(hpt & pred & ~nocue, :, 1);
y_HPN = temp(hpt & pred & ~nocue, :, 2);
[t_, a_] = timeaverage(t_HPN', y_HPN');
plot(t_, a_);
xlim([-.5 .5])

subplot(412); hold on
vy_HPN = temp(hpt & pred & ~nocue, :, 3);
[t_, a_] = timeaverage(t_HPN', vy_HPN');
plot(t_, a_);
xlim([-.5 .5])

subplot(413); hold on
ay_HPN = temp(hpt & pred & ~nocue, :, 4);
[t_, a_] = timeaverage(t_HPN', ay_HPN');
plot(t_, a_);
xlim([-.5 .5])

subplot(414); hold on
am_HPN = temp(hpt & pred & ~nocue, :, 5);
[t_, a_] = timeaverage(t_HPN', am_HPN');
plot(t_, a_);
xlim([-.5 .5])

figure;
subplot(411); hold on
t_HPN = temp(~hpt & ~pred & ~nocue, :, 1);
y_HPN = temp(~hpt & ~pred & ~nocue, :, 2);
[t_, a_] = timeaverage(t_HPN', y_HPN');
plot(t_, a_);
xlim([-.5 .5])

subplot(412); hold on
vy_HPN = temp(~hpt & ~pred & ~nocue, :, 3);
[t_, a_] = timeaverage(t_HPN', vy_HPN');
plot(t_, a_);
xlim([-.5 .5])

subplot(413); hold on
ay_HPN = temp(~hpt & ~pred & ~nocue, :, 4);
[t_, a_] = timeaverage(t_HPN', ay_HPN');
plot(t_, a_);
xlim([-.5 .5])

subplot(414); hold on
am_HPN = temp(~hpt & ~pred & ~nocue, :, 5);
[t_, a_] = timeaverage(t_HPN', am_HPN');
plot(t_, a_);
xlim([-.5 .5])


%%
% compute difference in power in desired freq. range for each trial type
% and for the first and the last blocks. Difference is taken relative to 
% Long-PT predictive trials in the last block (theoretically the most 
% proficient trial type). (note: for experiment conditions where there was
% no trial type in the last block, the first block will basically get 
% copied to the last block):
% % S_l_nc_1 = nanmean(H_h_p_2(k_f_start:k_f_end, :, end),1) - nanmean(H_l_nc_2(k_f_start:k_f_end, :, 1), 1);
% % S_l_np_1 = nanmean(H_h_p_2(k_f_start:k_f_end, :, end),1) - nanmean(H_l_np_2(k_f_start:k_f_end, :, 1), 1);
% % S_l_p_1 = nanmean(H_h_p_2(k_f_start:k_f_end, :, end),1) - nanmean(H_l_p_2(k_f_start:k_f_end, :, 1), 1);
% % 
% % S_l_nc_2 = nanmean(H_h_p_2(k_f_start:k_f_end, :, end),1) - nanmean(H_l_nc_2(k_f_start:k_f_end, :, end), 1);
% % S_l_np_2 = nanmean(H_h_p_2(k_f_start:k_f_end, :, end),1) - nanmean(H_l_np_2(k_f_start:k_f_end, :, end), 1);
% % S_l_p_2 = nanmean(H_h_p_2(k_f_start:k_f_end, :, end),1) - nanmean(H_l_p_2(k_f_start:k_f_end, :, end), 1);
% % 
% % S_h_nc_1 = nanmean(H_h_p_2(k_f_start:k_f_end, :, end),1) - nanmean(H_h_nc_2(k_f_start:k_f_end, :, 1), 1);
% % S_h_np_1 = nanmean(H_h_p_2(k_f_start:k_f_end, :, end),1) - nanmean(H_h_np_2(k_f_start:k_f_end, :, 1), 1);
% % % S_h_p_1 = nanmean(H_h_p_2(k_f_start:k_f_end, :, 1), 1);
% % 
% % S_h_nc_2 = nanmean(H_h_p_2(k_f_start:k_f_end, :, end),1) - nanmean(H_h_nc_2(k_f_start:k_f_end, :, end), 1);
% % S_h_np_2 = nanmean(H_h_p_2(k_f_start:k_f_end, :, end),1) - nanmean(H_h_np_2(k_f_start:k_f_end, :, end), 1);
% % % S_h_p_2 = nanmean(H_h_p_2(k_f_start:k_f_end, :, end), 1);
% % 
% % stats_results.high_pt.nc = [S_h_nc_1; S_h_nc_2]';
% % stats_results.high_pt.np = [S_h_np_1; S_h_np_2]';
% % % stats_results.high_pt.p = [S_h_p_1; S_h_p_2]';
% % stats_results.low_pt.nc = [S_l_nc_1; S_l_nc_2]';
% % stats_results.low_pt.np = [S_l_np_1; S_l_np_2]';
% % stats_results.low_pt.p = [S_l_p_1; S_l_p_2]';
% % 
% % % compute average difference for each trial type for first and last blocks:
% % m_mat = -[nanmean(stats_results.high_pt.nc);...
% %     nanmean(stats_results.high_pt.np);...
% %     nanmean(stats_results.low_pt.nc);...
% %     nanmean(stats_results.low_pt.np);...
% %     nanmean(stats_results.low_pt.p)];
% % s_mat = -[sqrt(nanvar(stats_results.high_pt.nc)./20);...
% %     sqrt(nanvar(stats_results.high_pt.np)./20);...
% %     sqrt(nanvar(stats_results.low_pt.nc)./20);...
% %     sqrt(nanvar(stats_results.low_pt.np)./20);...
% %     sqrt(nanvar(stats_results.low_pt.p)./20)];
% % 
% % % plot bars of the average difference for the last block:
% % figure; 
% % bar(m_mat(:,2), 'k'); hold on;
% % errorbar(1:length(m_mat), m_mat(:,2), s_mat(:,2), 'k.');
% % axis([0 7 -4 10])
% % a = gca;
% % a.XTickLabel = {'H-NC', 'H-NP', 'L-NC', 'L-NP', 'L-P'};
% % % fill in any subject values that are nan with the group mean:
% % % (this might be a bit dubious.. consider removing them instead, or filling
% % % in missing values based on a more proper missing value analysis)
% % temp_h_np_2 = S_h_np_2;
% % temp_h_np_2(isnan(temp_h_np_2)) = nanmean(temp_h_np_2);
% % 
% % % b = S_h_p_2;
% % % b(isnan(b)) = nanmean(b);
% % 
% % temp_l_np_2 = S_l_np_2;
% % temp_l_np_2(isnan(temp_l_np_2)) = nanmean(temp_l_np_2);
% % 
% % temp_l_p_2 = S_l_p_2;
% % temp_l_p_2(isnan(temp_l_p_2)) = nanmean(temp_l_p_2);
% % 
% % [h_aov1, d_aov1, p_aov1] = anova1([temp_h_np_2', temp_l_np_2', temp_l_p_2']);
% % 
% % [h_t, d_t, c_t, p_t] = ttest(temp_l_np_2 - temp_l_p_2);
% % 
% % stats_results.stats.aov = {h_aov1, d_aov1, p_aov1};
% % stats_results.stats.t = {h_t, d_t, c_t, p_t};

stats_results = [];

