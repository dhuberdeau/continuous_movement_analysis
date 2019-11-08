function stats_results = plot_results_learning(H_results, varargin)
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
% metric_all has columns: 1) metric-value; 2) pred; 3) hpt; 4) nocue; 5) total_trial; 6) block
hpt = H_results.metric_all(:, 3);
pred = H_results.metric_all(:, 2);
nocue = H_results.metric_all(:, 4);
block_ind = H_results.metric_all(:, 6);

hpt(isnan(hpt)) = 0;
pred(isnan(pred)) = 0;
nocue(isnan(nocue)) = 0;
hpt = boolean(hpt);
pred = boolean(pred);
nocue = boolean(nocue); % these are a combination that will never be analyzed (because it didn't happen in the exp.)
% block = boolean(block);

basis_cat = block_ind > 4 & hpt == 1 & pred == 1 & nocue == 0;
% t_basis = temp(basis_cat, :,1);
% y_basis = temp(basis_cat, :,2);
% v_basis = temp(basis_cat, :,3);
% a_basis = temp(basis_cat, :,4);
% m_basis = temp(basis_cat, :,5);
for i_block = 1:6
    for pt = 1:2
        for prd = 1:2
            nc = 1;
            sub_cat = block_ind == i_block & ...
                hpt == (pt - 1) & ...
                pred == (prd - 1) & ...
                nocue == (nc - 1);


            if nargin > 1
                plot_special(temp, sub_cat, basis_cat, varargin{1});
            else
                plot_special(temp, sub_cat, basis_cat, varargin{1});
            end
%             if sum(sub_cat > 0)
%                 t_temp = temp(sub_cat, :,1);
%                 y_temp = temp(sub_cat, :,2);
%                 v_temp = temp(sub_cat, :,3);
%                 a_temp = temp(sub_cat, :,4);
%                 m_temp = temp(sub_cat, :,5);
% 
%                 t_ = [t_basis; t_temp];
%                 y_ = [y_basis; y_temp];
%                 v_ = [v_basis; v_temp];
%                 a_ = [a_basis; a_temp];
%                 m_ = [m_basis; m_temp];
% 
%                 cat_ = [zeros(size(t_basis,1),1); ones(size(t_temp,1),1)];
% 
%                 [t_main, kin_0, kin_up, kin_down, kin_int] = ...
%                     fit_transform_movement_direction(t_, y_, v_, a_, m_, cat_);
% 
%                 % plot and save all on same axes:
%                 figure;
%                 for i_meas = 1:4
%                     subplot(1,4,i_meas); hold on;
% 
%                     if i_meas == 2
%                         title(['block:', num2str(i_block),...
%                             ' HPT:', num2str(pt), ...
%                             ' Pred:', num2str(prd), ...
%                             ' Cue:', num2str(nc)...
%                             ]);
%                     end
% 
%                     plot(t_main, nanmean(kin_0{i_meas}), 'k');
%                     plot(t_main, nanmean(kin_up{i_meas}), 'g');
%                     plot(t_main, nanmean(kin_down{i_meas}), 'r');
%                     plot(t_main, nanmean(kin_int{i_meas}), 'Color', [.5 .5 .5]);
%                     ylim(meas_limits{i_meas});
% 
%                 end
%                 f = gcf; 
%                 set(f, 'Position', [1 683 1679 272]);
%                 if nargin > 1
%                     saveas(f, ['targ_all_', varargin{1}, '.png'])
%                 end
%             end
            title(['block:', num2str(i_block),...
                    ' HPT:', num2str(pt-1), ...
                    ' Pred:', num2str(prd-1), ...
                    ' Cue:', num2str(mod(nc,2))...
                    ]);
        end
    end
end

for block_ind = 1:6
    for pt = 1:2
        nc = 2;
        sub_cat = block_ind == i_block & ...
                    hpt == (pt - 1) & ...
                    nocue == (nc - 1);

        if nargin > 1
            plot_special(temp, sub_cat, basis_cat, varargin{1});
        else
            plot_special(temp, sub_cat, basis_cat, varargin{1});
        end
        
        title(['block:', num2str(i_block),...
            ' HPT:', num2str(pt-1), ...
            ' Cue:', num2str(mod(nc,2))...
            ]);

    end
end


%%
stats_results = [];


function plot_special(temp, sub_cat, basis_cat, varargin)

meas_limits = {[-250 250], [-.2, .2], [-0.00031, 0.00031], [-0.00031 0.00031]};

t_basis = temp(basis_cat, :,1);
y_basis = temp(basis_cat, :,2);
v_basis = temp(basis_cat, :,3);
a_basis = temp(basis_cat, :,4);
m_basis = temp(basis_cat, :,5);

if sum(sub_cat > 0)
    t_temp = temp(sub_cat, :,1);
    y_temp = temp(sub_cat, :,2);
    v_temp = temp(sub_cat, :,3);
    a_temp = temp(sub_cat, :,4);
    m_temp = temp(sub_cat, :,5);

    t_ = [t_basis; t_temp];
    y_ = [y_basis; y_temp];
    v_ = [v_basis; v_temp];
    a_ = [a_basis; a_temp];
    m_ = [m_basis; m_temp];

    cat_ = [zeros(size(t_basis,1),1); ones(size(t_temp,1),1)];

    [t_main, kin_0, kin_up, kin_down, kin_int] = ...
        fit_transform_movement_direction(t_, y_, v_, a_, m_, cat_);

    % plot and save all on same axes:
    figure;
    for i_meas = 1:4
        subplot(1,4,i_meas); hold on;

        plot(t_main, nanmean(kin_0{i_meas}), 'k');
        plot(t_main, nanmean(kin_up{i_meas}), 'g');
        plot(t_main, nanmean(kin_down{i_meas}), 'r');
        plot(t_main, nanmean(kin_int{i_meas}), 'Color', [.5 .5 .5]);
        ylim(meas_limits{i_meas});

    end
    f = gcf; 
    set(f, 'Position', [1 683 1679 272]);
    if nargin > 1
        saveas(f, ['targ_all_', varargin{1}, '.png'])
    end
end
