function H_results = process_group
% process the data for all subjects, calculate features from data, and do
% statistics on those features.
%
% David Huberdeau, 04/10/2019
load('data_sets.mat')

% rand_set = randperm(20);
% data_sets = data_sets(rand_set(1:10));

b = cell(1, length(data_sets));
for i_sub = 1:length(data_sets)
    b{i_sub} = process_subject(data_sets{i_sub});
end
N_f = size(b{1}.h_means_1{1,1},1);
N_blocks = 5;
N_trials_block = 48;
N_FACTORS = 8;

H_l_nc_2 = nan(N_f, length(data_sets));
H_l_p_2 = nan(N_f, length(data_sets), N_blocks);
H_l_np_2 = nan(N_f, length(data_sets), N_blocks);
H_h_nc_2 = nan(N_f, length(data_sets));
H_h_p_2 = nan(N_f, length(data_sets), N_blocks);
H_h_np_2 = nan(N_f, length(data_sets), N_blocks);

H_metric_all = nan(length(data_sets)*(N_blocks+1)*N_trials_block, N_FACTORS+1);
Kinematics_all = nan(length(data_sets)*(N_blocks+1)*N_trials_block, 500, 8);
V_metric_all = nan(length(data_sets)*(N_blocks+1)*N_trials_block, N_FACTORS+1);
k_split_all = nan(length(data_sets)*(N_blocks+1)*N_trials_block, 1);
k_all = 1;

for i_sub = 1:length(data_sets)
    
    H_l_nc_2(:, i_sub) = b{i_sub}.h_means_2{2,3}(:, 1);
    for i_block = 1:(size(b{i_sub}.h_means_1{2,1}, 2)-2)
        H_l_p_2(:, i_sub, i_block) = b{i_sub}.h_means_2{2,1}(:, i_block + 1);
    end
    H_l_p_2(:, i_sub, end) = b{i_sub}.h_means_2{2,1}(:, end);
    for i_block = 1:(size(b{i_sub}.h_means_1{2,2}, 2)-2)
        H_l_np_2(:, i_sub, i_block) = b{i_sub}.h_means_2{2,2}(:, i_block + 1);
    end
    H_l_np_2(:, i_sub, end) = b{i_sub}.h_means_2{2,2}(:, end);
    
    
    H_h_nc_2(:, i_sub) = b{i_sub}.h_means_2{1,3}(:, 1);
    for i_block = 1:(size(b{i_sub}.h_means_1{1,1}, 2)-2)
        H_h_p_2(:, i_sub, i_block) = b{i_sub}.h_means_2{1,1}(:, i_block + 1);
    end
    H_h_p_2(:, i_sub, end) = b{i_sub}.h_means_2{1,1}(:, end);
    for i_block = 1:(size(b{i_sub}.h_means_1{1,2}, 2)-2)
        H_h_np_2(:, i_sub, i_block) = b{i_sub}.h_means_2{1,2}(:, i_block + 1);
    end
    H_h_np_2(:, i_sub, end) = b{i_sub}.h_means_2{1,2}(:, end);
    
    H_metric_all(k_all - 1 + (1:size(b{i_sub}.h_all,1)), 1:N_FACTORS) = b{i_sub}.h_all;
    H_metric_all(k_all - 1 + (1:size(b{i_sub}.h_all,1)), (N_FACTORS+1)) = i_sub;
    V_metric_all(k_all - 1 + (1:size(b{i_sub}.h_all,1)), 1:N_FACTORS) = b{i_sub}.v_err_all;
    V_metric_all(k_all - 1 + (1:size(b{i_sub}.h_all,1)), (N_FACTORS+1)) = i_sub;
    Kinematics_all(k_all - 1 + (1:size(b{i_sub}.h_all,1)), :, :) = b{i_sub}.kin_all;
    k_split_all(k_all - 1 + (1:size(b{i_sub}.h_all,1)), :, :) = b{i_sub}.k_split_all;
    
    k_all = k_all + size(b{i_sub}.h_all,1);
end

H_results.high_pt.nc = H_h_nc_2;
H_results.high_pt.np = H_h_np_2;
H_results.high_pt.p = H_h_p_2;
H_results.low_pt.nc = H_l_nc_2;
H_results.low_pt.np = H_l_np_2;
H_results.low_pt.p = H_l_p_2;
H_results.freq = b{1}.f;
H_results.metric_all = H_metric_all;
H_results.kin_all = Kinematics_all;
H_results.v_err_all = V_metric_all;
H_results.k_split_all = k_split_all;

%%
H_temp = H_results.metric_all;
csvwrite('H_metric_all_v1.txt', H_temp);

V_temp = H_results.v_err_all;
csvwrite('V_metric_all_v1.txt', V_temp);

% 
% figure; errorbar(b{1}.f, mean(H_hpt1,2), std(H_hpt1,[],2)./sqrt(length(data_sets))); 
% hold on; 
% errorbar(b{1}.f, mean(H_lpt1,2), std(H_lpt1,[],2)./sqrt(length(data_sets)));
% 
% figure; errorbar(b{1}.f, mean(H_hpt2,2), std(H_hpt2,[],2)./sqrt(length(data_sets))); 
% hold on; 
% errorbar(b{1}.f, mean(H_lpt2,2), std(H_lpt2,[],2)./sqrt(length(data_sets)));

%% plot by High or Low PT
% figure; 
% subplot(2,1,1); hold on
% plot(b{1}.f, nanmean(H_h_nc_2, 2), 'r-', 'Linewidth',2)
% plot(b{1}.f, nanmean(H_h_p_2(:,:,1), 2), 'b.-', 'Linewidth',1)
% plot(b{1}.f, nanmean(nanmean(H_h_p_2(:,:,2:(end-1)),3), 2), 'b--', 'Linewidth',2)
% plot(b{1}.f, nanmean(H_h_p_2(:,:,end), 2), 'b-', 'Linewidth',2)
% plot(b{1}.f, nanmean(H_h_np_2(:,:,1), 2), 'g.-', 'Linewidth',1)
% plot(b{1}.f, nanmean(nanmean(H_h_np_2(:,:,2:(end-1)),3), 2), 'g--', 'Linewidth',2)
% plot(b{1}.f, nanmean(H_h_np_2(:,:,end), 2), 'g-', 'Linewidth',2)
% title('High PT')
% axis([0 10 -110 -80])
% 
% subplot(2,1,2); hold on
% plot(b{1}.f, nanmean(H_l_nc_2, 2), 'r-', 'Linewidth',2)
% plot(b{1}.f, nanmean(H_l_p_2(:,:,1), 2), 'b.-', 'Linewidth',1)
% plot(b{1}.f, nanmean(nanmean(H_l_p_2(:,:,2:(end-1)),3), 2), 'b--', 'Linewidth',2)
% plot(b{1}.f, nanmean(H_l_p_2(:,:,end), 2), 'b-', 'Linewidth',2)
% plot(b{1}.f, nanmean(H_l_np_2(:,:,1), 2), 'g.-', 'Linewidth',1)
% plot(b{1}.f, nanmean(nanmean(H_l_np_2(:,:,2:(end-1)),3), 2), 'g--', 'Linewidth',2)
% plot(b{1}.f, nanmean(H_l_np_2(:,:,end), 2), 'g-', 'Linewidth',2)
% title('Low PT')
% axis([0 10 -110 -80])

%% plot by predictive or non-predictive 
% figure; 
% subplot(2,2,1); hold on
% plot(b{1}.f, nanmean(H_h_nc_2, 2), 'b-', 'Linewidth',2)
% plot(b{1}.f, nanmean(H_l_nc_2, 2), 'g-', 'Linewidth',2)
% title('No-cue')
% legend('high-pt / no-cue',...
%     'low-pt / no-cue');
% axis([0 10 -110 -80])
% 
% subplot(2,2,2); hold on
% % plot(b{1}.f, nanmean(H_h_nc_2, 2), 'r-', 'Linewidth',2)
% plot(b{1}.f, nanmean(H_h_p_2(:,:,1), 2), 'b.-', 'Linewidth',1)
% plot(b{1}.f, nanmean(nanmean(H_h_p_2(:,:,2:(end-1)),3), 2), 'b--', 'Linewidth',2)
% plot(b{1}.f, nanmean(H_h_p_2(:,:,end), 2), 'b-', 'Linewidth',2)
% plot(b{1}.f, nanmean(H_l_p_2(:,:,1), 2), 'g.-', 'Linewidth',1)
% plot(b{1}.f, nanmean(nanmean(H_l_p_2(:,:,2:(end-1)),3), 2), 'g--', 'Linewidth',2)
% plot(b{1}.f, nanmean(H_l_p_2(:,:,end), 2), 'g-', 'Linewidth',2)
% title('Predictive')
% legend('high-pt / block1', 'high-pt / mid-block', 'high-pt / known',...
%     'low-pt / block1', 'low-pt / mid-block', 'low-pt / known')
% axis([0 10 -110 -80])
% 
% 
% subplot(2,2,4); hold on
% plot(b{1}.f, nanmean(H_h_np_2(:,:,1), 2), 'b.-', 'Linewidth',1)
% plot(b{1}.f, nanmean(nanmean(H_h_np_2(:,:,2:(end-1)),3), 2), 'b--', 'Linewidth',2)
% plot(b{1}.f, nanmean(H_h_np_2(:,:,end), 2), 'b-', 'Linewidth',2)
% plot(b{1}.f, nanmean(H_l_np_2(:,:,1), 2), 'g.-', 'Linewidth',1)
% plot(b{1}.f, nanmean(nanmean(H_l_np_2(:,:,2:(end-1)),3), 2), 'g--', 'Linewidth',2)
% plot(b{1}.f, nanmean(H_l_np_2(:,:,end), 2), 'g-', 'Linewidth',2)
% title('Non-predictive')
% legend('high-pt / block1', 'high-pt / mid-block', 'high-pt / known',...
%     'low-pt / block1', 'low-pt / mid-block', 'low-pt / known')
% axis([0 10 -110 -80])
%% extract features and do stats:
% 
% f_feat = [2 5]; %ffrequency of interest 
% h_feat = nan(2, length(data_sets));
% h_feat(1, :) = mean(H_hpt1(b{1}.f >= f_feat(1) & b{1}.f <= f_feat(2), :), 1);
% h_feat(2, :) = mean(H_lpt1(b{1}.f >= f_feat(1) & b{1}.f <= f_feat(2), :), 1);
% 
% [a_1,b_1,c_1,d_1] = ttest(h_feat(2,:) - h_feat(1,:));
% 
% 
% f_feat = [2 3 4 5]; %ffrequency of interest 
% h_feat = nan(2, length(data_sets));
% h_feat(1, :) = mean(H_hpt2(b{1}.f >= f_feat(1) & b{1}.f <= f_feat(2), :), 1);
% h_feat(2, :) = mean(H_lpt2(b{1}.f >= f_feat(1) & b{1}.f <= f_feat(2), :), 1);
% 
% [a_2,b_2,c_2,d_2] = ttest(h_feat(2,:) - h_feat(1,:));


%% plot only last block:

% 
% figure; 
% subplot(2,1,1);
% hold on
% plot(b{1}.f, nanmean(H_h_nc_2(:,:,end), 2), 'k-', 'Linewidth',2)
% plot(b{1}.f, nanmean(H_l_nc_2(:,:,end), 2), 'r-', 'Linewidth',2)
% title('No-cue');
% legend('High-pt', 'Low-pt');
% subplot(2,1,2); hold on;
% plot(b{1}.f, nanmean(H_h_np_2(:,:,end) - H_l_np_2(:,:,end), 2),...
%     'r-', 'Linewidth',2)
% plot([0 30], [0 0], '-', 'Color', [.5 .5 .5]);
% title('No-cue');
% legend('Diff: High - Low');
% axis([0 10 -5 4])
% 
% 
% figure; 
% subplot(2,2,1); hold on
% plot(b{1}.f, nanmean(H_l_nc_2(:,:,end), 2), 'b-', 'Linewidth',2)
% plot(b{1}.f, nanmean(H_h_p_2(:,:,1), 2), 'k.-', 'Linewidth',1)
% plot(b{1}.f, nanmean(H_h_p_2(:,:,end), 2), 'k-', 'Linewidth',2)
% plot(b{1}.f, nanmean(H_l_p_2(:,:,1), 2), 'r.-', 'Linewidth',1)
% plot(b{1}.f, nanmean(H_l_p_2(:,:,end), 2), 'r-', 'Linewidth',2)
% title('Predictive')
% legend('low-pt / no-cue', 'high-pt / block1', 'high-pt / known',...
%     'low-pt / block1', 'low-pt / known')
% axis([0 10 -120 -80])
% 
% subplot(2,2,2); hold on
% plot(b{1}.f, nanmean(H_l_nc_2(:,:,end), 2), 'b-', 'Linewidth',2)
% plot(b{1}.f, nanmean(H_h_np_2(:,:,1), 2), 'k.-', 'Linewidth',1)
% plot(b{1}.f, nanmean(H_h_np_2(:,:,end), 2), 'k-', 'Linewidth',2)
% plot(b{1}.f, nanmean(H_l_np_2(:,:,1), 2), 'r.-', 'Linewidth',1)
% plot(b{1}.f, nanmean(H_l_np_2(:,:,end), 2), 'r-', 'Linewidth',2)
% title('Non-predictive')
% legend('low-pt / no-cue', 'high-pt / block1', 'high-pt / known',...
%     'low-pt / block1', 'low-pt / known')
% axis([0 10 -120 -80])
% 
% subplot(2,2,3); hold on;
% plot(b{1}.f, nanmean(H_h_p_2(:,:,end) - H_l_nc_2(:,:,end), 2), 'b-', 'Linewidth',1)
% plot(b{1}.f, nanmean(H_h_p_2(:,:,end) - H_h_p_2(:,:,1), 2), 'k.-', 'Linewidth',1)
% plot(b{1}.f, nanmean(H_h_p_2(:,:,end) - H_l_p_2(:,:,1), 2), 'r.-', 'Linewidth',1)
% plot(b{1}.f, nanmean(H_h_p_2(:,:,end) - H_l_p_2(:,:,end), 2), 'r-', 'Linewidth',2)
% title('Diff: predictive')
% plot([0 30], [0 0], '-', 'Color', [.5 .5 .5]);
% axis([0 10 -5 4])
% 
% subplot(2,2,4); hold on;
% plot(b{1}.f, nanmean(H_h_np_2(:,:,end) - H_l_nc_2(:,:,end), 2), 'b-', 'Linewidth',1)
% plot(b{1}.f, nanmean(H_h_np_2(:,:,end) - H_h_np_2(:,:,1), 2), 'k.-', 'Linewidth',1)
% 
% plot(b{1}.f, nanmean(H_h_np_2(:,:,end) - H_l_np_2(:,:,1), 2), 'r.-', 'Linewidth',1)
% plot(b{1}.f, nanmean(H_h_np_2(:,:, end) - H_l_np_2(:,:,end), 2), 'r-', 'Linewidth',2)
% title('Diff: Non-predictive')
% plot([0 30], [0 0], '-', 'Color', [.5 .5 .5]);
% axis([0 10 -5 4])