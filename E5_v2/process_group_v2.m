% process the data for all subjects, calculate features from data, and do
% statistics on those features.
%
% David Huberdeau, 04/10/2019
load('data_sets_v2.mat')

b = cell(1, length(data_sets));
for i_sub = 1:length(data_sets)
    b{i_sub} = process_subject(data_sets{i_sub});
end
N_f = size(b{1}.h_means_1{1,1},1);
N_blocks = 5;

H_l_nc_2 = nan(N_f, length(data_sets));
H_l_p_2 = nan(N_f, length(data_sets), N_blocks);
H_l_np_2 = nan(N_f, length(data_sets), N_blocks);
H_h_nc_2 = nan(N_f, length(data_sets));
H_h_p_2 = nan(N_f, length(data_sets), N_blocks);
H_h_np_2 = nan(N_f, length(data_sets), N_blocks);
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
end
% 
% figure; errorbar(b{1}.f, mean(H_hpt1,2), std(H_hpt1,[],2)./sqrt(length(data_sets))); 
% hold on; 
% errorbar(b{1}.f, mean(H_lpt1,2), std(H_lpt1,[],2)./sqrt(length(data_sets)));
% 
% figure; errorbar(b{1}.f, mean(H_hpt2,2), std(H_hpt2,[],2)./sqrt(length(data_sets))); 
% hold on; 
% errorbar(b{1}.f, mean(H_lpt2,2), std(H_lpt2,[],2)./sqrt(length(data_sets)));

%% plot by High or Low PT
figure; 
subplot(2,1,1); hold on
plot(b{1}.f, nanmean(H_h_nc_2, 2), 'r-', 'Linewidth',2)
plot(b{1}.f, nanmean(H_h_p_2(:,:,1), 2), 'b.-', 'Linewidth',1)
plot(b{1}.f, nanmean(nanmean(H_h_p_2(:,:,2:(end-1)),3), 2), 'b--', 'Linewidth',2)
plot(b{1}.f, nanmean(H_h_p_2(:,:,end), 2), 'b-', 'Linewidth',2)
plot(b{1}.f, nanmean(H_h_np_2(:,:,1), 2), 'g.-', 'Linewidth',1)
plot(b{1}.f, nanmean(nanmean(H_h_np_2(:,:,2:(end-1)),3), 2), 'g--', 'Linewidth',2)
plot(b{1}.f, nanmean(H_h_np_2(:,:,end), 2), 'g-', 'Linewidth',2)
title('High PT')
axis([0 10 -110 -80])

subplot(2,1,2); hold on
plot(b{1}.f, nanmean(H_l_nc_2, 2), 'r-', 'Linewidth',2)
plot(b{1}.f, nanmean(H_l_p_2(:,:,1), 2), 'b.-', 'Linewidth',1)
plot(b{1}.f, nanmean(nanmean(H_l_p_2(:,:,2:(end-1)),3), 2), 'b--', 'Linewidth',2)
plot(b{1}.f, nanmean(H_l_p_2(:,:,end), 2), 'b-', 'Linewidth',2)
plot(b{1}.f, nanmean(H_l_np_2(:,:,1), 2), 'g.-', 'Linewidth',1)
plot(b{1}.f, nanmean(nanmean(H_l_np_2(:,:,2:(end-1)),3), 2), 'g--', 'Linewidth',2)
plot(b{1}.f, nanmean(H_l_np_2(:,:,end), 2), 'g-', 'Linewidth',2)
title('Low PT')
axis([0 10 -110 -80])

%% plot by predictive or non-predictive 
figure; 
subplot(2,2,1); hold on
plot(b{1}.f, nanmean(H_h_nc_2, 2), 'b-', 'Linewidth',2)
plot(b{1}.f, nanmean(H_l_nc_2, 2), 'g-', 'Linewidth',2)
title('No-cue')
legend('high-pt / no-cue',...
    'low-pt / no-cue');
axis([0 10 -110 -80])

subplot(2,2,2); hold on
% plot(b{1}.f, nanmean(H_h_nc_2, 2), 'r-', 'Linewidth',2)
plot(b{1}.f, nanmean(H_h_p_2(:,:,1), 2), 'b.-', 'Linewidth',1)
plot(b{1}.f, nanmean(nanmean(H_h_p_2(:,:,2:(end-1)),3), 2), 'b--', 'Linewidth',2)
plot(b{1}.f, nanmean(H_h_p_2(:,:,end), 2), 'b-', 'Linewidth',2)
plot(b{1}.f, nanmean(H_l_p_2(:,:,1), 2), 'g.-', 'Linewidth',1)
plot(b{1}.f, nanmean(nanmean(H_l_p_2(:,:,2:(end-1)),3), 2), 'g--', 'Linewidth',2)
plot(b{1}.f, nanmean(H_l_p_2(:,:,end), 2), 'g-', 'Linewidth',2)
title('Predictive')
legend('high-pt / block1', 'high-pt / mid-block', 'high-pt / known',...
    'low-pt / block1', 'low-pt / mid-block', 'low-pt / known')
axis([0 10 -110 -80])


subplot(2,2,4); hold on
plot(b{1}.f, nanmean(H_h_np_2(:,:,1), 2), 'b.-', 'Linewidth',1)
plot(b{1}.f, nanmean(nanmean(H_h_np_2(:,:,2:(end-1)),3), 2), 'b--', 'Linewidth',2)
plot(b{1}.f, nanmean(H_h_np_2(:,:,end), 2), 'b-', 'Linewidth',2)
plot(b{1}.f, nanmean(H_l_np_2(:,:,1), 2), 'g.-', 'Linewidth',1)
plot(b{1}.f, nanmean(nanmean(H_l_np_2(:,:,2:(end-1)),3), 2), 'g--', 'Linewidth',2)
plot(b{1}.f, nanmean(H_l_np_2(:,:,end), 2), 'g-', 'Linewidth',2)
title('Non-predictive')
legend('high-pt / block1', 'high-pt / mid-block', 'high-pt / known',...
    'low-pt / block1', 'low-pt / mid-block', 'low-pt / known')
axis([0 10 -110 -80])
%% extract features and do stats:
% 
% f_feat = [2 3 4 5]; %ffrequency of interest 
% h_feat = nan(2, length(data_sets));
% h_feat(1, :) = mean(H_hpt1(f_feat, :), 1);
% h_feat(2, :) = mean(H_lpt1(f_feat, :), 1);
% 
% [a,b,c,d] = ttest(h_feat(2,:) - h_feat(1,:));
% 
% 
% f_feat = [2 3 4 5]; %ffrequency of interest 
% h_feat = nan(2, length(data_sets));
% h_feat(1, :) = mean(H_hpt2(f_feat, :), 1);
% h_feat(2, :) = mean(H_lpt2(f_feat, :), 1);
% 
% [a_,b_,c_,d_] = ttest(h_feat(2,:) - h_feat(1,:));
