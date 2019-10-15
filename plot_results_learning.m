function stats_results = plot_results(H_results, varargin)
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
% metric_all has columns: 1)metric-value; 2)pred; 3)hpt; 4)nocue; 5)total_trial; 6)block
hpt = H_results.metric_all(:, 3);
pred = H_results.metric_all(:, 2);
nocue = H_results.metric_all(:, 4);
block_ind = H_results.metric_all(:, 6);
block_fin = block_ind == nanmax(unique(block_ind)); % pick out last block trials;
block_1 = block_ind < 3;
block_2 = block_ind >= 3 & ~block_fin;


hpt(isnan(hpt)) = 0;
pred(isnan(pred)) = 0;
nocue(isnan(nocue)) = 0;
block_fin(isnan(block_fin)) = 0;
block_1(isnan(block_1)) = 0;
block_2(isnan(block_2)) = 0;
hpt = boolean(hpt);
pred = boolean(pred);
nocue = boolean(nocue); % these are a combination that will never be analyzed (because it didn't happen in the exp.)
% block = boolean(block);

HPC1_cat = hpt & pred & ~nocue & block_1;
HPC2_cat = hpt & pred & ~nocue & block_2;
LPC1_cat = ~hpt & pred & ~nocue & block_1;
LPC2_cat = ~hpt & pred & ~nocue & block_2;

HUC1_cat = hpt & ~pred & ~nocue & block_1;
HUC2_cat = hpt & ~pred & ~nocue & block_2;
LUC1_cat = ~hpt & ~pred & ~nocue & block_1;
LUC2_cat = ~hpt & ~pred & ~nocue & block_2;

HPC_fin_cat = hpt & pred & ~nocue & block_fin;
LPC_fin_cat = ~hpt & pred & ~nocue & block_fin;

HUC_fin_cat = hpt & ~pred & ~nocue & block_fin;
LUC_fin_cat = ~hpt & ~pred & ~nocue & block_fin;

t_HPC = temp(HPC1_cat, :, 1);
y_HPC = temp(HPC1_cat, :, 2);
vy_HPC = temp(HPC1_cat, :, 3);
ay_HPC = temp(HPC1_cat, :, 4);
am_HPC = temp(HPC1_cat, :, 5);

t_HPC2 = temp(HPC2_cat, :, 1);
y_HPC2 = temp(HPC2_cat, :, 2);
vy_HPC2 = temp(HPC2_cat, :, 3);
ay_HPC2 = temp(HPC2_cat, :, 4);
am_HPC2 = temp(HPC2_cat, :, 5);

t_LPC = temp(LPC1_cat, :, 1);
y_LPC = temp(LPC1_cat, :, 2);
vy_LPC = temp(LPC1_cat, :, 3);
ay_LPC = temp(LPC1_cat, :, 4);
am_LPC = temp(LPC1_cat, :, 5);

t_LPC2 = temp(LPC2_cat, :, 1);
y_LPC2 = temp(LPC2_cat, :, 2);
vy_LPC2 = temp(LPC2_cat, :, 3);
ay_LPC2 = temp(LPC2_cat, :, 4);
am_LPC2 = temp(LPC2_cat, :, 5);

t_HUC = temp(HUC1_cat, :, 1);
y_HUC = temp(HUC1_cat, :, 2);
vy_HUC = temp(HUC1_cat, :, 3);
ay_HUC = temp(HUC1_cat, :, 4);
am_HUC = temp(HUC1_cat, :, 5);

t_HUC2 = temp(HUC2_cat, :, 1);
y_HUC2 = temp(HUC2_cat, :, 2);
vy_HUC2 = temp(HUC2_cat, :, 3);
ay_HUC2 = temp(HUC2_cat, :, 4);
am_HUC2 = temp(HUC2_cat, :, 5);

t_LUC = temp(LUC1_cat, :, 1);
y_LUC = temp(LUC1_cat, :, 2);
vy_LUC = temp(LUC1_cat, :, 3);
ay_LUC = temp(LUC1_cat, :, 4);
am_LUC = temp(LUC1_cat, :, 5);

t_LUC2 = temp(LUC2_cat, :, 1);
y_LUC2 = temp(LUC2_cat, :, 2);
vy_LUC2 = temp(LUC2_cat, :, 3);
ay_LUC2 = temp(LUC2_cat, :, 4);
am_LUC2 = temp(LUC2_cat, :, 5);

% t_HN = temp(hpt & nocue, :, 1);
% y_HN = temp(hpt & nocue, :, 2);
% vy_HN = temp(hpt & nocue, :, 3);
% ay_HN = temp(hpt & nocue, :, 4);
% am_HN = temp(hpt & nocue, :, 5);
% 
% t_LN = temp(~hpt & nocue, :, 1);
% y_LN = temp(~hpt & nocue, :, 2);
% vy_LN = temp(~hpt & nocue, :, 3);
% ay_LN = temp(~hpt & nocue, :, 4);
% am_LN = temp(~hpt & nocue, :, 5);

t_HPC_fin = temp(HPC_fin_cat, :, 1);
y_HPC_fin = temp(HPC_fin_cat, :, 2);
vy_HPC_fin = temp(HPC_fin_cat, :, 3);
ay_HPC_fin = temp(HPC_fin_cat, :, 4);
am_HPC_fin = temp(HPC_fin_cat, :, 5);

t_LPC_fin = temp(LPC_fin_cat, :, 1);
y_LPC_fin = temp(LPC_fin_cat, :, 2);
vy_LPC_fin = temp(LPC_fin_cat, :, 3);
ay_LPC_fin = temp(LPC_fin_cat, :, 4);
am_LPC_fin = temp(LPC_fin_cat, :, 5);

t_HUC_fin = temp(HUC_fin_cat, :, 1);
y_HUC_fin = temp(HUC_fin_cat, :, 2);
vy_HUC_fin = temp(HUC_fin_cat, :, 3);
ay_HUC_fin = temp(HUC_fin_cat, :, 4);
am_HUC_fin = temp(HUC_fin_cat, :, 5);

t_LUC_fin = temp(LUC_fin_cat, :, 1);
y_LUC_fin = temp(LUC_fin_cat, :, 2);
vy_LUC_fin = temp(LUC_fin_cat, :, 3);
ay_LUC_fin = temp(LUC_fin_cat, :, 4);
am_LUC_fin = temp(LUC_fin_cat, :, 5);

% t_LUN = temp(~hpt & ~pred & ~nocue, :, 1);
% y_LUN = temp(~hpt & ~pred & ~nocue, :, 2);
% vy_LUN = temp(~hpt & ~pred & ~nocue, :, 3);
% ay_LUN = temp(~hpt & ~pred & ~nocue, :, 4);
% am_LUN = temp(~hpt & ~pred & ~nocue, :, 5);

% t_N = [t_HPC; t_LUC; t_HN; t_LN; t_HPC_fin];
% y_N = [y_HPC; y_LUC; y_HN; y_LN; y_HPC_fin];
% vy_N = [vy_HPC; vy_LUC; vy_HN; vy_LN; vy_HPC_fin];
% ay_N = [ay_HPC; ay_LUC; ay_HN; ay_LN; ay_HPC_fin];
% am_N = [am_HPC; am_LUC; am_HN; am_LN; am_HPC_fin];
% cat = [zeros(size(t_HPC,1),1); zeros(size(t_LUC,1),1);...
%     zeros(size(t_HN,1),1); zeros(size(t_LN,1),1); ones(size(t_HPC_fin,1),1)];
% [t_main, y_interp, v_interp, a_interp, m_interp, cat_out] = ...
%     fit_transform_movement_direction(t_N, y_N, vy_N, ay_N, am_N, cat);

t_N = [t_HPC; t_HPC2; t_LPC; t_LPC2; t_HUC; t_HUC2; t_LUC; t_LUC2; t_HUC_fin; t_LUC_fin; t_LPC_fin; t_HPC_fin];
y_N = [y_HPC; y_HPC2; y_LPC; y_LPC2; y_HUC; y_HUC2; y_LUC; y_LUC2; y_HUC_fin; y_LUC_fin; y_LPC_fin; y_HPC_fin];
vy_N = [vy_HPC; vy_HPC2; vy_LPC; vy_LPC2; vy_HUC; vy_HUC2; vy_LUC; vy_LUC2; vy_HUC_fin; vy_LUC_fin; vy_LPC_fin; vy_HPC_fin];
ay_N = [ay_HPC; ay_HPC2; ay_LPC; ay_LPC2; ay_HUC; ay_HUC2; ay_LUC; ay_LUC2; ay_HUC_fin; ay_LUC_fin; ay_LPC_fin; ay_HPC_fin];
am_N = [am_HPC; am_HPC2; am_LPC; am_LPC2; am_HUC; am_HUC2; am_LUC; am_LUC2; am_HUC_fin; am_LUC_fin; am_LPC_fin; am_HPC_fin];
cat = [zeros(size(t_HPC,1),1); zeros(size(t_HPC2,1),1); zeros(size(t_LPC,1),1); zeros(size(t_LPC2,1),1);...
    zeros(size(t_HUC,1),1); zeros(size(t_HUC2,1),1); zeros(size(t_LUC,1),1); zeros(size(t_LUC2,1),1); ...
    zeros(size(t_HUC_fin,1),1); zeros(size(t_LUC_fin,1),1); zeros(size(t_LPC_fin,1),1); ones(size(t_HPC_fin,1),1)];

[t_main, y_interp, v_interp, a_interp, m_interp, cat_out] = ...
    fit_transform_movement_direction(t_N, y_N, vy_N, ay_N, am_N, cat);

y_0 = y_interp(cat_out(:,1), :);
v_0 = v_interp(cat_out(:,1), :);
a_0 = a_interp(cat_out(:,1), :);
ma_0 = m_interp(cat_out(:,1), :);
kin_0 = {y_0, v_0, a_0, ma_0};

y_HPC1_up = y_interp(cat_out(:,2) & HPC1_cat, :);
v_HPC1_up = v_interp(cat_out(:,2) & HPC1_cat, :);
a_HPC1_up = a_interp(cat_out(:,2) & HPC1_cat, :);
ma_HPC1_up = m_interp(cat_out(:,2) & HPC1_cat, :);
kin_HPC1_up = {y_HPC1_up, v_HPC1_up, a_HPC1_up, ma_HPC1_up};

y_HPC1_down = y_interp(cat_out(:,3) & HPC1_cat, :);
v_HPC1_down = v_interp(cat_out(:,3) & HPC1_cat, :);
a_HPC1_down = a_interp(cat_out(:,3) & HPC1_cat, :);
ma_HPC1_down = m_interp(cat_out(:,3) & HPC1_cat, :);
kin_HPC1_down = {y_HPC1_down, v_HPC1_down, a_HPC1_down, ma_HPC1_down};

y_HPC1_int = y_interp(cat_out(:,4) & HPC1_cat, :);
v_HPC1_int = v_interp(cat_out(:,4) & HPC1_cat, :);
a_HPC1_int = a_interp(cat_out(:,4) & HPC1_cat, :);
ma_HPC1_int = m_interp(cat_out(:,4) & HPC1_cat, :);
kin_HPC1_int = {y_HPC1_int, v_HPC1_int, a_HPC1_int, ma_HPC1_int};

y_HPC2_up = y_interp(cat_out(:,2) & HPC2_cat, :);
v_HPC2_up = v_interp(cat_out(:,2) & HPC2_cat, :);
a_HPC2_up = a_interp(cat_out(:,2) & HPC2_cat, :);
ma_HPC2_up = m_interp(cat_out(:,2) & HPC2_cat, :);
kin_HPC2_up = {y_HPC2_up, v_HPC2_up, a_HPC2_up, ma_HPC2_up};

y_HPC2_down = y_interp(cat_out(:,3) & HPC2_cat, :);
v_HPC2_down = v_interp(cat_out(:,3) & HPC2_cat, :);
a_HPC2_down = a_interp(cat_out(:,3) & HPC2_cat, :);
ma_HPC2_down = m_interp(cat_out(:,3) & HPC2_cat, :);
kin_HPC2_down = {y_HPC2_down, v_HPC2_down, a_HPC2_down, ma_HPC2_down};

y_HPC2_int = y_interp(cat_out(:,4) & HPC2_cat, :);
v_HPC2_int = v_interp(cat_out(:,4) & HPC2_cat, :);
a_HPC2_int = a_interp(cat_out(:,4) & HPC2_cat, :);
ma_HPC2_int = m_interp(cat_out(:,4) & HPC2_cat, :);
kin_HPC2_int = {y_HPC2_int, v_HPC2_int, a_HPC2_int, ma_HPC2_int};

y_LPC1_up = y_interp(cat_out(:,2) & LPC1_cat, :);
v_LPC1_up = v_interp(cat_out(:,2) & LPC1_cat, :);
a_LPC1_up = a_interp(cat_out(:,2) & LPC1_cat, :);
ma_LPC1_up = m_interp(cat_out(:,2) & LPC1_cat, :);
kin_LPC1_up = {y_LPC1_up, v_LPC1_up, a_LPC1_up, ma_LPC1_up};

y_LPC1_down = y_interp(cat_out(:,3) & LPC1_cat, :);
v_LPC1_down = v_interp(cat_out(:,3) & LPC1_cat, :);
a_LPC1_down = a_interp(cat_out(:,3) & LPC1_cat, :);
ma_LPC1_down = m_interp(cat_out(:,3) & LPC1_cat, :);
kin_LPC1_down = {y_LPC1_down, v_LPC1_down, a_LPC1_down, ma_LPC1_down};

y_LPC1_int = y_interp(cat_out(:,4) & LPC1_cat, :);
v_LPC1_int = v_interp(cat_out(:,4) & LPC1_cat, :);
a_LPC1_int = a_interp(cat_out(:,4) & LPC1_cat, :);
ma_LPC1_int = m_interp(cat_out(:,4) & LPC1_cat, :);
kin_LPC1_int = {y_LPC1_int, v_LPC1_int, a_LPC1_int, ma_LPC1_int};

y_LPC2_up = y_interp(cat_out(:,2) & LPC2_cat, :);
v_LPC2_up = v_interp(cat_out(:,2) & LPC2_cat, :);
a_LPC2_up = a_interp(cat_out(:,2) & LPC2_cat, :);
ma_LPC2_up = m_interp(cat_out(:,2) & LPC2_cat, :);
kin_LPC2_up = {y_LPC2_up, v_LPC2_up, a_LPC2_up, ma_LPC2_up};

y_LPC2_down = y_interp(cat_out(:,3) & LPC2_cat, :);
v_LPC2_down = v_interp(cat_out(:,3) & LPC2_cat, :);
a_LPC2_down = a_interp(cat_out(:,3) & LPC2_cat, :);
ma_LPC2_down = m_interp(cat_out(:,3) & LPC2_cat, :);
kin_LPC2_down = {y_LPC2_down, v_LPC2_down, a_LPC2_down, ma_LPC2_down};

y_LPC2_int = y_interp(cat_out(:,4) & LPC2_cat, :);
v_LPC2_int = v_interp(cat_out(:,4) & LPC2_cat, :);
a_LPC2_int = a_interp(cat_out(:,4) & LPC2_cat, :);
ma_LPC2_int = m_interp(cat_out(:,4) & LPC2_cat, :);
kin_LPC2_int = {y_LPC2_int, v_LPC2_int, a_LPC2_int, ma_LPC2_int};

y_HUC1_up = y_interp(cat_out(:,2) & HPC1_cat, :);
v_HUC1_up = v_interp(cat_out(:,2) & HPC1_cat, :);
a_HUC1_up = a_interp(cat_out(:,2) & HPC1_cat, :);
ma_HUC1_up = m_interp(cat_out(:,2) & HPC1_cat, :);
kin_HUC1_up = {y_HUC1_up, v_HUC1_up, a_HUC1_up, ma_HUC1_up};

y_HUC1_down = y_interp(cat_out(:,3) & HPC1_cat, :);
v_HUC1_down = v_interp(cat_out(:,3) & HPC1_cat, :);
a_HUC1_down = a_interp(cat_out(:,3) & HPC1_cat, :);
ma_HUC1_down = m_interp(cat_out(:,3) & HPC1_cat, :);
kin_HUC1_down = {y_HUC1_down, v_HUC1_down, a_HUC1_down, ma_HUC1_down};

y_HUC1_int = y_interp(cat_out(:,4) & HPC1_cat, :);
v_HUC1_int = v_interp(cat_out(:,4) & HPC1_cat, :);
a_HUC1_int = a_interp(cat_out(:,4) & HPC1_cat, :);
ma_HUC1_int = m_interp(cat_out(:,4) & HPC1_cat, :);
kin_HUC1_int = {y_HUC1_int, v_HUC1_int, a_HUC1_int, ma_HUC1_int};

y_HUC2_up = y_interp(cat_out(:,2) & HUC2_cat, :);
v_HUC2_up = v_interp(cat_out(:,2) & HUC2_cat, :);
a_HUC2_up = a_interp(cat_out(:,2) & HUC2_cat, :);
ma_HUC2_up = m_interp(cat_out(:,2) & HUC2_cat, :);
kin_HUC2_up = {y_HUC2_up, v_HUC2_up, a_HUC2_up, ma_HUC2_up};

y_HUC2_down = y_interp(cat_out(:,3) & HUC2_cat, :);
v_HUC2_down = v_interp(cat_out(:,3) & HUC2_cat, :);
a_HUC2_down = a_interp(cat_out(:,3) & HUC2_cat, :);
ma_HUC2_down = m_interp(cat_out(:,3) & HUC2_cat, :);
kin_HUC2_down = {y_HUC2_down, v_HUC2_down, a_HUC2_down, ma_HUC2_down};

y_HUC2_int = y_interp(cat_out(:,4) & HUC2_cat, :);
v_HUC2_int = v_interp(cat_out(:,4) & HUC2_cat, :);
a_HUC2_int = a_interp(cat_out(:,4) & HUC2_cat, :);
ma_HUC2_int = m_interp(cat_out(:,4) & HUC2_cat, :);
kin_HUC2_int = {y_HUC2_int, v_HUC2_int, a_HUC2_int, ma_HUC2_int};

y_LUC1_up = y_interp(cat_out(:,2) & LUC1_cat, :);
v_LUC1_up = v_interp(cat_out(:,2) & LUC1_cat, :);
a_LUC1_up = a_interp(cat_out(:,2) & LUC1_cat, :);
ma_LUC1_up = m_interp(cat_out(:,2) & LUC1_cat, :);
kin_LUC1_up = {y_LUC1_up, v_LUC1_up, a_LUC1_up, ma_LUC1_up};

y_LUC1_down = y_interp(cat_out(:,3) & LUC1_cat, :);
v_LUC1_down = v_interp(cat_out(:,3) & LUC1_cat, :);
a_LUC1_down = a_interp(cat_out(:,3) & LUC1_cat, :);
ma_LUC1_down = m_interp(cat_out(:,3) & LUC1_cat, :);
kin_LUC1_down = {y_LUC1_down, v_LUC1_down, a_LUC1_down, ma_LUC1_down};

y_LUC1_int = y_interp(cat_out(:,4) & LUC1_cat, :);
v_LUC1_int = v_interp(cat_out(:,4) & LUC1_cat, :);
a_LUC1_int = a_interp(cat_out(:,4) & LUC1_cat, :);
ma_LUC1_int = m_interp(cat_out(:,4) & LUC1_cat, :);
kin_LUC1_int = {y_LUC1_int, v_LUC1_int, a_LUC1_int, ma_LUC1_int};

y_LUC2_up = y_interp(cat_out(:,2) & LUC2_cat, :);
v_LUC2_up = v_interp(cat_out(:,2) & LUC2_cat, :);
a_LUC2_up = a_interp(cat_out(:,2) & LUC2_cat, :);
ma_LUC2_up = m_interp(cat_out(:,2) & LUC2_cat, :);
kin_LUC2_up = {y_LUC2_up, v_LUC2_up, a_LUC2_up, ma_LUC2_up};

y_LUC2_down = y_interp(cat_out(:,3) & LUC2_cat, :);
v_LUC2_down = v_interp(cat_out(:,3) & LUC2_cat, :);
a_LUC2_down = a_interp(cat_out(:,3) & LUC2_cat, :);
ma_LUC2_down = m_interp(cat_out(:,3) & LUC2_cat, :);
kin_LUC2_down = {y_LUC2_down, v_LUC2_down, a_LUC2_down, ma_LUC2_down};

y_LUC2_int = y_interp(cat_out(:,4) & LUC2_cat, :);
v_LUC2_int = v_interp(cat_out(:,4) & LUC2_cat, :);
a_LUC2_int = a_interp(cat_out(:,4) & LUC2_cat, :);
ma_LUC2_int = m_interp(cat_out(:,4) & LUC2_cat, :);
kin_LUC2_int = {y_LUC2_int, v_LUC2_int, a_LUC2_int, ma_LUC2_int};



y_HPC_fin_up = y_interp(cat_out(:,2) & HPC_fin_cat, :);
v_HPC_fin_up = v_interp(cat_out(:,2) & HPC_fin_cat, :);
a_HPC_fin_up = a_interp(cat_out(:,2) & HPC_fin_cat, :);
ma_HPC_fin_up = m_interp(cat_out(:,2) & HPC_fin_cat, :);
kin_HPC_fin_up = {y_HPC_fin_up, v_HPC_fin_up, a_HPC_fin_up, ma_HPC_fin_up};

y_HPC_fin_down = y_interp(cat_out(:,3) & HPC_fin_cat, :);
v_HPC_fin_down = v_interp(cat_out(:,3) & HPC_fin_cat, :);
a_HPC_fin_down = a_interp(cat_out(:,3) & HPC_fin_cat, :);
ma_HPC_fin_down = m_interp(cat_out(:,3) & HPC_fin_cat, :);
kin_HPC_fin_down = {y_HPC_fin_down, v_HPC_fin_down, a_HPC_fin_down, ma_HPC_fin_down};

y_HPC_fin_int = y_interp(cat_out(:,4) & HPC_fin_cat, :);
v_HPC_fin_int = v_interp(cat_out(:,4) & HPC_fin_cat, :);
a_HPC_fin_int = a_interp(cat_out(:,4) & HPC_fin_cat, :);
ma_HPC_fin_int = m_interp(cat_out(:,4) & HPC_fin_cat, :);
kin_HPC_fin_int = {y_HPC_fin_int, v_HPC_fin_int, a_HPC_fin_int, ma_HPC_fin_int};

y_LPC_fin_up = y_interp(cat_out(:,2) & LPC_fin_cat, :);
v_LPC_fin_up = v_interp(cat_out(:,2) & LPC_fin_cat, :);
a_LPC_fin_up = a_interp(cat_out(:,2) & LPC_fin_cat, :);
ma_LPC_fin_up = m_interp(cat_out(:,2) & LPC_fin_cat, :);
kin_LPC_fin_up = {y_LPC_fin_up, v_LPC_fin_up, a_LPC_fin_up, ma_LPC_fin_up};

y_LPC_fin_down = y_interp(cat_out(:,3) & LPC_fin_cat, :);
v_LPC_fin_down = v_interp(cat_out(:,3) & LPC_fin_cat, :);
a_LPC_fin_down = a_interp(cat_out(:,3) & LPC_fin_cat, :);
ma_LPC_fin_down = m_interp(cat_out(:,3) & LPC_fin_cat, :);
kin_LPC_fin_down = {y_LPC_fin_down, v_LPC_fin_down, a_LPC_fin_down, ma_LPC_fin_down};

y_LPC_fin_int = y_interp(cat_out(:,4) & LPC_fin_cat, :);
v_LPC_fin_int = v_interp(cat_out(:,4) & LPC_fin_cat, :);
a_LPC_fin_int = a_interp(cat_out(:,4) & LPC_fin_cat, :);
ma_LPC_fin_int = m_interp(cat_out(:,4) & LPC_fin_cat, :);
kin_LPC_fin_int = {y_LPC_fin_int, v_LPC_fin_int, a_LPC_fin_int, ma_LPC_fin_int};

y_HUC_fin_up = y_interp(cat_out(:,2) & HUC_fin_cat, :);
v_HUC_fin_up = v_interp(cat_out(:,2) & HUC_fin_cat, :);
a_HUC_fin_up = a_interp(cat_out(:,2) & HUC_fin_cat, :);
ma_HUC_fin_up = m_interp(cat_out(:,2) & HUC_fin_cat, :);
kin_HUC_fin_up = {y_HUC_fin_up, v_HUC_fin_up, a_HUC_fin_up, ma_HUC_fin_up};

y_HUC_fin_down = y_interp(cat_out(:,3) & HUC_fin_cat, :);
v_HUC_fin_down = v_interp(cat_out(:,3) & HUC_fin_cat, :);
a_HUC_fin_down = a_interp(cat_out(:,3) & HUC_fin_cat, :);
ma_HUC_fin_down = m_interp(cat_out(:,3) & HUC_fin_cat, :);
kin_HUC_fin_down = {y_HUC_fin_down, v_HUC_fin_down, a_HUC_fin_down, ma_HUC_fin_down};

y_HUC_fin_int = y_interp(cat_out(:,4) & HUC_fin_cat, :);
v_HUC_fin_int = v_interp(cat_out(:,4) & HUC_fin_cat, :);
a_HUC_fin_int = a_interp(cat_out(:,4) & HUC_fin_cat, :);
ma_HUC_fin_int = m_interp(cat_out(:,4) & HUC_fin_cat, :);
kin_HUC_fin_int = {y_HUC_fin_int, v_HUC_fin_int, a_HUC_fin_int, ma_HUC_fin_int};

y_LUC_fin_up = y_interp(cat_out(:,2) & LUC_fin_cat, :);
v_LUC_fin_up = v_interp(cat_out(:,2) & LUC_fin_cat, :);
a_LUC_fin_up = a_interp(cat_out(:,2) & LUC_fin_cat, :);
ma_LUC_fin_up = m_interp(cat_out(:,2) & LUC_fin_cat, :);
kin_LUC_fin_up = {y_LUC_fin_up, v_LUC_fin_up, a_LUC_fin_up, ma_LUC_fin_up};

y_LUC_fin_down = y_interp(cat_out(:,3) & LUC_fin_cat, :);
v_LUC_fin_down = v_interp(cat_out(:,3) & LUC_fin_cat, :);
a_LUC_fin_down = a_interp(cat_out(:,3) & LUC_fin_cat, :);
ma_LUC_fin_down = m_interp(cat_out(:,3) & LUC_fin_cat, :);
kin_LUC_fin_down = {y_LUC_fin_down, v_LUC_fin_down, a_LUC_fin_down, ma_LUC_fin_down};

y_LUC_fin_int = y_interp(cat_out(:,4) & LUC_fin_cat, :);
v_LUC_fin_int = v_interp(cat_out(:,4) & LUC_fin_cat, :);
a_LUC_fin_int = a_interp(cat_out(:,4) & LUC_fin_cat, :);
ma_LUC_fin_int = m_interp(cat_out(:,4) & LUC_fin_cat, :);
kin_LUC_fin_int = {y_LUC_fin_int, v_LUC_fin_int, a_LUC_fin_int, ma_LUC_fin_int};





y_down = y_interp(cat_out(:,3), :);
v_down = v_interp(cat_out(:,3), :);
a_down = a_interp(cat_out(:,3), :);
ma_down = m_interp(cat_out(:,3), :);
kin_down = {y_down, v_down, a_down, ma_down};

y_int = y_interp(cat_out(:,4), :);
v_int = v_interp(cat_out(:,4), :);
a_int = a_interp(cat_out(:,4), :);
ma_int = m_interp(cat_out(:,4), :);
kin_int = {y_int, v_int, a_int, ma_int};

meas_limits = {[-250 250], [-.2, .2], [-0.00031, 0.00031], [-0.00031 0.00031]};

figure;
for i_meas = 1:4
    subplot(1,4,i_meas)
    plot(t_main, nanmean(kin_0{i_meas}));
    ylim(meas_limits{i_meas});
end
f = gcf; 
set(f, 'Position', [1 683 1679 272]);
if nargin > 1
    saveas(f, ['targ0_', varargin{1}, '.png'])
end

figure;
for i_meas = 1:4
    subplot(1,4,i_meas)
    plot(t_main, nanmean(kin_up{i_meas}));
    ylim(meas_limits{i_meas});
end
f = gcf; 
set(f, 'Position', [1 683 1679 272]);
if nargin > 1
    saveas(f, ['targ1_', varargin{1}, '.png'])
end

figure;
for i_meas = 1:4
    subplot(1,4,i_meas)
    plot(t_main, nanmean(kin_down{i_meas}));
    ylim(meas_limits{i_meas});
end
f = gcf; 
set(f, 'Position', [1 683 1679 272]);
if nargin > 1
    saveas(f, ['targ2_', varargin{1}, '.png'])
end

figure;
for i_meas = 1:4
    subplot(1,4,i_meas)
    plot(t_main, nanmean(kin_int{i_meas}));
    ylim(meas_limits{i_meas});
end
f = gcf; 
set(f, 'Position', [1 683 1679 272]);
if nargin > 1
    saveas(f, ['targ_intermediate_', varargin{1}, '.png'])
end

% plot and save all on same axes:
figure;
for i_meas = 1:4
    subplot(1,4,i_meas); hold on;
    
    plot(t_main, nanmean(kin_0{i_meas}), 'k');
%     plot(t_main, nanmean(kin_up{i_meas}), 'g');
%     plot(t_main, nanmean(kin_down{i_meas}), 'r');
%     plot(t_main, nanmean(kin_int{i_meas}), 'Color', [.5 .5 .5]);
    ylim(meas_limits{i_meas});
end
f = gcf; 
set(f, 'Position', [1 683 1679 272]);
if nargin > 1
    saveas(f, ['targ_all_', varargin{1}, '.png'])
end

% figure;
% subplot(411); hold on
% [t_, a_] = timeaverage(t_HPN', y_HPN');
% plot(t_, a_);
% xlim([-.5 .5])
% 
% subplot(412); hold on
% [t_, a_] = timeaverage(t_HPN', vy_HPN');
% plot(t_, a_);
% xlim([-.5 .5])
% 
% subplot(413); hold on
% [t_, a_] = timeaverage(t_HPN', ay_HPN');
% plot(t_, a_);
% xlim([-.5 .5])
% 
% subplot(414); hold on
% [t_, a_] = timeaverage(t_HPN', am_HPN');
% plot(t_, a_);
% xlim([-.5 .5])
% 
% figure;
% subplot(411); hold on
% [t_, a_] = timeaverage(t_LUN', y_LUN');
% plot(t_, a_);
% xlim([-.5 .5])
% 
% subplot(412); hold on
% [t_, a_] = timeaverage(t_LUN', vy_LUN');
% plot(t_, a_);
% xlim([-.5 .5])
% 
% subplot(413); hold on
% [t_, a_] = timeaverage(t_LUN', ay_LUN');
% plot(t_, a_);
% xlim([-.5 .5])
% 
% subplot(414); hold on
% [t_, a_] = timeaverage(t_LUN', am_LUN');
% plot(t_, a_);
% xlim([-.5 .5])


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

function plot_special(result_cell)


