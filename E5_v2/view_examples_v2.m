function varargout = view_examples_v2(data_set, varargin)
% view_examples
%
% optional inputs:
%  (1) block
% 
% David Huberdeau, 04/15/2019

pred_combo = [...
    5 5 6 6 7 7 8 8; ...
    3 1 4 2 1 3 2 4];

unpred_combo = [...
    5 6 7 8;...
    4 1 2 3];
    
nocue_blocks = 1;

% plot_matrix = [1 5 9 13 17 21 25;...
%     2 6 10 14 18 22 26;...
%     3 7 11 15 19 23 27;...
%     4 8 12 16 20 24 28];

plot_matrix = [...
    1 2 3 4 5 6 7;...
    8 9 10 11 12 13 14;...
    ];

colors_matrix = {'k', 'g'; 'b', 'c'}; 

Fs = 60;
dist_th = 70;
f_all = .5:.5:(Fs/2);
pt_th = .325;

% compute features:
%   - error at 70 pixels from center
%   - psd of second half
%   - psd of first half
h_means_1 = cell(2, 3); %rows: L- vs. H-PT, cols: pred, unpred, or no-cue
for i_pt = 1:2
    for i_cue = 1:3
        h_means_1{i_pt, i_cue} = nan(length(f_all), length(data_set));
    end
end

h_means_2 = cell(2, 2);
for i_pt = 1:2
    for i_cue = 1:3
        h_means_2{i_pt, i_cue} = nan(length(f_all), length(data_set));
    end
end

figure;
for i_block = 1:length(data_set)
    try 
        % process block:
        b_data = process_block(data_set{i_block});

        % compute actual PT:
        pt = (1/Fs)*(b_data.k_split - b_data.k_targ);
        
        lpt_tr = pt < pt_th;
        hpt_tr = pt >= pt_th;
        
        combo = [data_set{i_block}.params.trial_home_numbers; data_set{i_block}.params.trial_stimA_numbers];
        
        pred_tr = zeros(1, length(b_data.k_split));
        unpred_tr = zeros(1, length(b_data.k_split));
        nocue_tr = zeros(1, length(b_data.k_split));
        
        if ismember(nocue_blocks, i_block)
            nocue_tr = ones(1, length(b_data.k_split));
        else
            for i_tr = 1:length(pred_tr)
                pred_check = nan(1, size(pred_combo,2));
                for i_pred = 1:size(pred_combo,2)
                    pred_check(i_pred) = isequal(pred_combo(:, i_pred), combo(:, i_tr));
                end
                unpred_check = nan(1, size(unpred_combo,2));
                for i_unpred = 1:size(unpred_combo,2)
                    unpred_check(i_unpred) = isequal(unpred_combo(:, i_unpred), combo(:, i_tr));
                end
                if sum(pred_check) > 0
                    pred_tr(i_tr) = 1;
                end
                if sum(unpred_check) > 0
                    unpred_tr(i_tr) = 1;
                end
                nocue_tr(i_tr) = ~pred_tr(i_tr) & ~unpred_tr(i_tr);
            end
        end
        pred_tr = logical(pred_tr);
        unpred_tr = logical(unpred_tr);
        nocue_tr = logical(nocue_tr);
        
        tr_inds = 1:length(b_data.k_split);
        
        %% plot each block and each type {hpt, lpt, pred, unpred}
        if i_block == 1
            subplot(2, 7, plot_matrix(1, i_block)); hold on;
            x_temp =  b_data.x; y_temp = b_data.y;
            plot(x_temp(:, hpt_tr & nocue_tr), y_temp(:, hpt_tr & nocue_tr), colors_matrix{1, 1});

            subplot(2, 7, plot_matrix(2, i_block)); hold on;
            x_temp =  b_data.x; y_temp = b_data.y;
            plot(x_temp(:, lpt_tr & nocue_tr), y_temp(:, lpt_tr & nocue_tr), colors_matrix{2, 1});
        
        elseif i_block < 6

            subplot(2, 7, plot_matrix(1, i_block)); hold on;
            x_temp =  b_data.x; y_temp = b_data.y;
            plot(x_temp(:, hpt_tr & pred_tr), y_temp(:, hpt_tr & pred_tr), colors_matrix{1, 1});

            subplot(2, 7, plot_matrix(1, i_block)); hold on;
            x_temp =  b_data.x; y_temp = b_data.y;
            plot(x_temp(:, hpt_tr & unpred_tr), y_temp(:, hpt_tr & unpred_tr), colors_matrix{1, 2});

            subplot(2, 7, plot_matrix(2, i_block)); hold on;
            x_temp =  b_data.x; y_temp = b_data.y;
            plot(x_temp(:, lpt_tr & pred_tr), y_temp(:, lpt_tr & pred_tr), colors_matrix{2, 1});

            subplot(2, 7, plot_matrix(2, i_block)); hold on;
            x_temp =  b_data.x; y_temp = b_data.y;
            plot(x_temp(:, lpt_tr & unpred_tr), y_temp(:, lpt_tr & unpred_tr), colors_matrix{2, 2});
        
            
        else
            subplot(2, 7, plot_matrix(1, i_block)); hold on;
            x_temp =  b_data.x; y_temp = b_data.y;
            plot(x_temp(:, hpt_tr & pred_tr), y_temp(:, hpt_tr & pred_tr), colors_matrix{1, 1});

            subplot(2, 7, plot_matrix(1, i_block)); hold on;
            x_temp =  b_data.x; y_temp = b_data.y;
            plot(x_temp(:, hpt_tr & unpred_tr), y_temp(:, hpt_tr & unpred_tr), colors_matrix{1, 2});

            subplot(2, 7, plot_matrix(2, i_block)); hold on;
            x_temp =  b_data.x; y_temp = b_data.y;
            plot(x_temp(:, lpt_tr & pred_tr), y_temp(:, lpt_tr & pred_tr), colors_matrix{2, 1});

            subplot(2, 7, plot_matrix(2, i_block)); hold on;
            x_temp =  b_data.x; y_temp = b_data.y;
            plot(x_temp(:, lpt_tr & unpred_tr), y_temp(:, lpt_tr & unpred_tr), colors_matrix{2, 2});
            
            
            subplot(2, 7, plot_matrix(1, i_block+1)); hold on;
            x_temp =  b_data.x; y_temp = b_data.y;
            plot(x_temp(:, hpt_tr & nocue_tr), y_temp(:, hpt_tr & nocue_tr), colors_matrix{1, 1});

            subplot(2, 7, plot_matrix(2, i_block+1)); hold on;
            x_temp =  b_data.x; y_temp = b_data.y;
            plot(x_temp(:, lpt_tr & nocue_tr), y_temp(:, lpt_tr & nocue_tr), colors_matrix{2, 1});
            
        end
        
    catch
        
    end
end
varargout = {};