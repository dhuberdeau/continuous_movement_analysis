function H_results = process_group_v4
% process the data for all subjects, calculate features from data, and do
% statistics on those features.
%
% David Huberdeau, 05/01/2019
load('data_sets_v4.mat')

% data_sets = data_sets([1:6, 13]); %symbol block 1st
% data_sets = data_sets([7:12, 14]); %no-cue block 1st

b = cell(1, length(data_sets));
for i_sub = 1:length(data_sets)
    b{i_sub} = process_subject(data_sets{i_sub});
end
N_f = size(b{1}.h_means_1{1,1},1);

N_blocks_learn = 3;
blocks_learn = [2 3 4];
N_blocks_practice = 1;
blocks_practice = [1];
N_blocks_known = 1;
blocks_known = [5];
N_blocks = N_blocks_known + N_blocks_learn;
N_trials_block = 48;
N_FACTORS = 8;

H_l_nc_2 = nan(N_f, length(data_sets), N_blocks_practice + N_blocks_known);
H_l_p_2 = nan(N_f, length(data_sets), N_blocks_learn + N_blocks_known);
H_l_np_2 = nan(N_f, length(data_sets), N_blocks_learn + N_blocks_known);
H_h_nc_2 = nan(N_f, length(data_sets), N_blocks_practice + N_blocks_known);
H_h_p_2 = nan(N_f, length(data_sets), N_blocks_learn + N_blocks_known);
H_h_np_2 = nan(N_f, length(data_sets), N_blocks_learn + N_blocks_known);

H_metric_all = nan(length(data_sets)*N_blocks*N_trials_block, 7);
V_metric_all = nan(length(data_sets)*(N_blocks+1)*N_trials_block, 7);
Kinematics_all = nan(length(data_sets)*(N_blocks+1)*N_trials_block, 500, 8);
k_split_all = nan(length(data_sets)*(N_blocks+1)*N_trials_block, 1);
k_all = 1;
for i_sub = 1:length(data_sets)
    
    blks = union(blocks_practice, blocks_known);
    for i_block = 1:length(blks)
        try
            H_l_nc_2(:, i_sub, i_block) = b{i_sub}.h_means_2{2,3}(:, blks(i_block));
        catch
        end
    end
    blks = union(blocks_learn, blocks_known);
    for i_block = 1:length(blks)
        try
            H_l_p_2(:, i_sub, i_block) = b{i_sub}.h_means_2{2,1}(:, blks(i_block));
        catch
        end
    end
    blks = union(blocks_learn, blocks_known);
    for i_block = 1:length(blks)
        try
            H_l_np_2(:, i_sub, i_block) = b{i_sub}.h_means_2{2,2}(:, blks(i_block));
        catch
        end
    end
    
    blks = union(blocks_practice, blocks_known);
    for i_block = 1:length(blks)
        try
            H_h_nc_2(:, i_sub, i_block) = b{i_sub}.h_means_2{1,3}(:, blks(i_block));
        catch
        end
    end
    blks = union(blocks_learn, blocks_known);
    for i_block = 1:length(blks)
        try
            H_h_p_2(:, i_sub, i_block) = b{i_sub}.h_means_2{1,1}(:, blks(i_block));
        catch
        end
    end
    blks = union(blocks_learn, blocks_known);
    for i_block = 1:length(blks)
        try
            H_h_np_2(:, i_sub, i_block) = b{i_sub}.h_means_2{1,2}(:, blks(i_block));
        catch
        end
    end
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
csvwrite('H_metric_all_v4.txt', H_temp);

V_temp = H_results.v_err_all;
csvwrite('V_metric_all_v4.txt', V_temp);

%% plot by High or Low PT
% figure; 
% subplot(2,1,1); hold on
% plot(b{1}.f, nanmean(H_h_nc_2(:,:,1), 2), 'r.-', 'Linewidth',1)
% plot(b{1}.f, nanmean(H_h_nc_2(:,:,2), 2), 'r-', 'Linewidth',2)
% 
% plot(b{1}.f, nanmean(H_h_p_2(:,:,1), 2), 'b.-', 'Linewidth',1)
% plot(b{1}.f, nanmean(nanmean(H_h_p_2(:,:,2:3),3), 2), 'b--', 'Linewidth',2)
% plot(b{1}.f, nanmean(H_h_p_2(:,:,end), 2), 'b-', 'Linewidth',2)
% 
% plot(b{1}.f, nanmean(H_h_np_2(:,:,1), 2), 'g.-', 'Linewidth',1)
% plot(b{1}.f, nanmean(nanmean(H_h_np_2(:,:,2:3),3), 2), 'g--', 'Linewidth',2)
% plot(b{1}.f, nanmean(H_h_np_2(:,:,end), 2), 'g-', 'Linewidth',2)
% title('High PT')
% axis([0 10 -110 -80])
% 
% subplot(2,1,2); hold on
% plot(b{1}.f, nanmean(H_l_nc_2(:,:,1), 2), 'r.-', 'Linewidth',1)
% plot(b{1}.f, nanmean(H_l_nc_2(:,:,2), 2), 'r-', 'Linewidth',2)
% 
% plot(b{1}.f, nanmean(H_l_p_2(:,:,1), 2), 'b.-', 'Linewidth',1)
% plot(b{1}.f, nanmean(nanmean(H_l_p_2(:,:,2:3),3), 2), 'b--', 'Linewidth',2)
% plot(b{1}.f, nanmean(H_l_p_2(:,:,end), 2), 'b-', 'Linewidth',2)
% 
% plot(b{1}.f, nanmean(H_l_np_2(:,:,1), 2), 'g.-', 'Linewidth',1)
% plot(b{1}.f, nanmean(nanmean(H_l_np_2(:,:,2:3),3), 2), 'g--', 'Linewidth',2)
% plot(b{1}.f, nanmean(H_l_np_2(:,:,end), 2), 'g-', 'Linewidth',2)
% title('Low PT')
% axis([0 10 -110 -80])

%% plot by predictive or non-predictive 
% figure; 
% subplot(2,2,1); hold on
% plot(b{1}.f, nanmean(H_h_nc_2(:,:,2), 2), 'b-', 'Linewidth',2)
% plot(b{1}.f, nanmean(H_l_nc_2(:,:,2), 2), 'g-', 'Linewidth',2)
% title('No-cue')
% legend('high-pt / no-cue',...
%     'low-pt / no-cue');
% axis([0 10 -110 -80])
% 
% subplot(2,2,2); hold on
% % plot(b{1}.f, nanmean(H_h_nc_2, 2), 'r-', 'Linewidth',2)
% plot(b{1}.f, nanmean(H_h_p_2(:,:,1), 2), 'b.-', 'Linewidth',1)
% % plot(b{1}.f, nanmean(nanmean(H_h_p_2(:,:,2:(end-1)),3), 2), 'b--', 'Linewidth',2)
% plot(b{1}.f, nanmean(H_h_p_2(:,:,end), 2), 'b-', 'Linewidth',2)
% plot(b{1}.f, nanmean(H_l_p_2(:,:,1), 2), 'g.-', 'Linewidth',1)
% % plot(b{1}.f, nanmean(nanmean(H_l_p_2(:,:,2:(end-1)),3), 2), 'g--', 'Linewidth',2)
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
% % plot(b{1}.f, nanmean(H_h_nc_2, 2), 'k-', 'Linewidth',2)
% plot(b{1}.f, nanmean(H_h_nc_2(:,:,end), 2), 'k-', 'Linewidth',2)
% % plot(b{1}.f, nanmean(H_l_nc_2, 2), 'r-', 'Linewidth',1)
% plot(b{1}.f, nanmean(H_l_nc_2(:,:,end), 2), 'r-', 'Linewidth',2)
% title('No-cue');
% legend('High-pt', 'Low-pt');
% subplot(2,1,2); hold on;
% plot(b{1}.f, nanmean(H_h_np_2(:,:,end) - H_l_np_2(:,:,end), 2),...
%     'r-', 'Linewidth',2)
% plot([0 30], [0 0], '-', 'Color', [.5 .5 .5]);
% title('No-cue');
% legend('Diff: High - Low');
% axis([0 20 -5 4])
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
% legend('no-cue', 'high-pt / block1', 'high-pt / known',...
%     'low-pt / block1', 'low-pt / known')
% axis([0 20 -120 -80])
% 
% subplot(2,2,2); hold on
% plot(b{1}.f, nanmean(H_l_nc_2(:,:,end), 2), 'b-', 'Linewidth',2)
% plot(b{1}.f, nanmean(H_h_np_2(:,:,1), 2), 'k.-', 'Linewidth',1)
% plot(b{1}.f, nanmean(H_h_np_2(:,:,end), 2), 'k-', 'Linewidth',2)
% plot(b{1}.f, nanmean(H_l_np_2(:,:,1), 2), 'r.-', 'Linewidth',1)
% plot(b{1}.f, nanmean(H_l_np_2(:,:,end), 2), 'r-', 'Linewidth',2)
% title('Non-predictive')
% legend('no-cue', 'high-pt / block1', 'high-pt / known',...
%     'low-pt / block1', 'low-pt / known')
% axis([0 20 -120 -80])
% 
% subplot(2,2,3); hold on;
% plot(b{1}.f, nanmean(H_h_p_2(:,:,end) - H_l_nc_2(:,:,end), 2), 'b-', 'Linewidth',1)
% plot(b{1}.f, nanmean(H_h_p_2(:,:,end) - H_h_p_2(:,:,1), 2), 'k.-', 'Linewidth',1)
% plot(b{1}.f, nanmean(H_h_p_2(:,:,end) - H_l_p_2(:,:,1), 2), 'r.-', 'Linewidth',1)
% plot(b{1}.f, nanmean(H_h_p_2(:,:,end) - H_l_p_2(:,:,end), 2), 'r-', 'Linewidth',2)
% title('Diff: predictive')
% plot([0 30], [0 0], '-', 'Color', [.5 .5 .5]);
% axis([0 20 -5 4])
% 
% subplot(2,2,4); hold on;
% plot(b{1}.f, nanmean(H_h_np_2(:,:,end) - H_l_nc_2(:,:,end), 2), 'b-', 'Linewidth',1)
% plot(b{1}.f, nanmean(H_h_np_2(:,:,end) - H_h_np_2(:,:,1), 2), 'k.-', 'Linewidth',1)
% 
% plot(b{1}.f, nanmean(H_h_np_2(:,:,end) - H_l_np_2(:,:,1), 2), 'r.-', 'Linewidth',1)
% plot(b{1}.f, nanmean(H_h_np_2(:,:, end) - H_l_np_2(:,:,end), 2), 'r-', 'Linewidth',2)
% title('Diff: Non-predictive')
% plot([0 30], [0 0], '-', 'Color', [.5 .5 .5]);
% axis([0 20 -5 4])