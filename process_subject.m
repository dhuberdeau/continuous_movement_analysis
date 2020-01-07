function s_data = process_subject(data_set)

% function s_data = process_subject(data_set)
%
% INPUT:
%   - data_set - a cell array of Data structures, one from each block
%
% OUTPUT:
%   - s_data - all raw data from this subject, with fields t, x,  y, vx,
%               vy, ax, ay, am, k_split
%
% David Huberdeau, 04/10/2019

% load('data_sets.mat'); %this command must be called by this function's
% caller.

% Define the combinations of start positions and symbols that together make
%   a trial predictive (pred_combo) or un-predictive (unpred_combo):
pred_combo = [...
    5 5 6 6 7 7 8 8; ...
    3 1 4 2 1 3 2 4];

unpred_combo = [...
    5 6 7 8;...
    4 1 2 3];

flip_combos = [...
    5 6 7 8;...
    3 4 1 2];

Fs = 60;
dist_th = 70;
TARG_DIR = 45;
f_all = .5:.5:(Fs/2);
pt_th = .325;
max_signal_len = 500;

N_FACTORS = 8;
% factors: [metric, pred vs unpred, numeric PT, Binary HPT, Cue, Trial,
% Block, Rewarded, trial needs flipped (i.e. target)]

% compute features:
%   - error at 70 pixels from center
%   - psd of second half
%   - psd of first half
h_means_1 = cell(2, 3); %rows: L- vs. H-PT, cols: pred, unpred, or no-cue
h_metric_trials = nan(length(data_set)*48, N_FACTORS); %columns: metric, pred vs unpred, numeric PT, Binary Low vs high pt, cue vs. no cue, trial, block, trial-needing-trajectory-flipped
kin_metric_trials = nan(length(data_set)*48, max_signal_len, 5); %trials match with factors in h_metric_trials 
v_metric_trials = nan(length(data_set)*48, N_FACTORS);
k_split_trials = nan(length(data_set)*48, 1);
            % (trials, signal-length, [time, y, vy, ay, |a|)
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

v_err = nan(4, length(data_set)); %4 for cue high PT, cue low PT, no cue, and known cue
total_trial_k = 1;
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
        flip_tr = zeros(1, length(b_data.k_split));
        reward_tr = zeros(1, length(b_data.k_split));
        
        if isfield(data_set{i_block}, 'Reward_assigned')
            reward_tr = data_set{i_block}.Reward_assigned;
        end
    
        % label each trial as being of type P (predictive), NP
        % (unpredictive) or NC (no cue):
        for i_tr = 1:length(pred_tr)
            pred_check = nan(1, size(pred_combo,2));
            for i_pred = 1:size(pred_combo,2)
                pred_check(i_pred) = isequal(pred_combo(:, i_pred), combo(:, i_tr));
            end
            unpred_check = nan(1, size(unpred_combo,2));
            for i_unpred = 1:size(unpred_combo,2)
                unpred_check(i_unpred) = isequal(unpred_combo(:, i_unpred), combo(:, i_tr));
            end
            flip_check = nan(1, size(flip_combos,2));
            for i_flip = 1:size(flip_combos, 2)
                flip_check(i_flip) = isequal(flip_combos(:, i_flip), combo(:, i_tr));
            end
            if sum(pred_check) > 0
                pred_tr(i_tr) = 1;
            end
            if sum(unpred_check) > 0
                unpred_tr(i_tr) = 1;
            end
            if sum(flip_check) > 0
                flip_tr(i_tr) = 1;
            end
            nocue_tr(i_tr) = ~pred_tr(i_tr) & ~unpred_tr(i_tr);
        end
        
        % supercede above designation if it is known that trial type is 1
        % (no cue type):
        nocue_tr(data_set{i_block}.params.trial_type == 1) = 1;
        pred_tr(data_set{i_block}.params.trial_type == 1) = 0;
        unpred_tr(data_set{i_block}.params.trial_type == 1) = 0;

        pred_tr = logical(pred_tr);
        unpred_tr = logical(unpred_tr);
        nocue_tr = logical(nocue_tr);
        flip_tr = logical(flip_tr);

        tr_inds = 1:length(b_data.k_split);

        %% extract PSD: from high-pt and first half of movement

        % for pred:
        [f1_hpt,h1_hpt] = compute_psd(b_data, Fs, hpt_tr & pred_tr, 1); %first half
        h1_ = nan(length(f_all), size(f1_hpt,2));
        for i_tr = 1:length(tr_inds(pred_tr & hpt_tr))
            try
            h1_(:, i_tr) = interp1(f1_hpt(~isnan(f1_hpt(:,i_tr)), i_tr), ...
                h1_hpt(~isnan(f1_hpt(:,i_tr)), i_tr), f_all);
            catch er__
                warning('a');
            end
        end
        h_means_1{1, 1}(:, i_block) = nanmean(h1_,2);


        % for unpred:
        [f1_hpt,h1_hpt] = compute_psd(b_data, Fs, hpt_tr & unpred_tr, 1); %first half
        h1_ = nan(length(f_all), size(f1_hpt,2));
        for i_tr = 1:length(tr_inds(unpred_tr & hpt_tr))
            try
            h1_(:, i_tr) = interp1(f1_hpt(~isnan(f1_hpt(:,i_tr)), i_tr), ...
                h1_hpt(~isnan(f1_hpt(:,i_tr)), i_tr), f_all);
            catch er__
                warning('a');
            end
        end
        h_means_1{1, 2}(:, i_block) = nanmean(h1_,2);

        % for nocue:
        [f1_hpt,h1_hpt] = compute_psd(b_data, Fs, hpt_tr & nocue_tr, 1); %first half
        h1_ = nan(length(f_all), size(f1_hpt,2));
        for i_tr = 1:length(tr_inds(nocue_tr & hpt_tr))
            try
                h1_(:, i_tr) = interp1(f1_hpt(~isnan(f1_hpt(:,i_tr)), i_tr), ...
                    h1_hpt(~isnan(f1_hpt(:,i_tr)), i_tr), f_all);
            catch er__
                warning('a');
            end
        end
        h_means_1{1, 3}(:, i_block) = nanmean(h1_,2);

        %% extract PSD: from low-pt and first half of movement

        % for pred:
        [f1_hpt,h1_hpt] = compute_psd(b_data, Fs, pred_tr & lpt_tr, 1); %first half
        h1_ = nan(length(f_all), size(f1_hpt,2));
        for i_tr = 1:length(tr_inds(pred_tr & lpt_tr))
            try
            h1_(:, i_tr) = interp1(f1_hpt(~isnan(f1_hpt(:,i_tr)), i_tr), ...
                h1_hpt(~isnan(f1_hpt(:,i_tr)), i_tr), f_all);
            catch er__
                warning('a');
            end
        end
        h_means_1{2, 1}(:, i_block) = nanmean(h1_,2);

        % for unpred:
        [f1_hpt,h1_hpt] = compute_psd(b_data, Fs, unpred_tr & lpt_tr, 1); %first half
        h1_ = nan(length(f_all), size(f1_hpt,2));
        for i_tr = 1:length(tr_inds(unpred_tr & lpt_tr))
            try
            h1_(:, i_tr) = interp1(f1_hpt(~isnan(f1_hpt(:,i_tr)), i_tr), ...
                h1_hpt(~isnan(f1_hpt(:,i_tr)), i_tr), f_all);
            catch er__
                warning('a');
            end
        end
        h_means_1{2, 2}(:, i_block) = nanmean(h1_,2);

        % for no-cue:
        [f1_hpt,h1_hpt] = compute_psd(b_data, Fs, nocue_tr & lpt_tr, 1); %first half
        h1_ = nan(length(f_all), size(f1_hpt,2));
        for i_tr = 1:length(tr_inds(nocue_tr & lpt_tr))
            try
            h1_(:, i_tr) = interp1(f1_hpt(~isnan(f1_hpt(:,i_tr)), i_tr), ...
                h1_hpt(~isnan(f1_hpt(:,i_tr)), i_tr), f_all);
            catch er__
                warning('a');
            end
        end
        h_means_1{2, 3}(:, i_block) = nanmean(h1_,2);

        %% extract PSD: from high-pt and all of movement

        % for pred:
        [f1_hpt,h1_hpt] = compute_psd(b_data, Fs, pred_tr & hpt_tr, 0); %all of movement
        h1_ = nan(length(f_all), size(f1_hpt,2));
        for i_tr = 1:length(tr_inds(pred_tr & hpt_tr))
            try
            h1_(:, i_tr) = interp1(f1_hpt(~isnan(f1_hpt(:,i_tr)), i_tr), ...
                h1_hpt(~isnan(f1_hpt(:,i_tr)), i_tr), f_all);
            catch er__
                warning('a');
            end
        end
        h_means_2{1, 1}(:, i_block) = nanmean(h1_,2);

        % for unpred:
        [f1_hpt,h1_hpt] = compute_psd(b_data, Fs, unpred_tr & hpt_tr, 0); %all of movement
        h1_ = nan(length(f_all), size(f1_hpt,2));
        for i_tr = 1:length(tr_inds(unpred_tr & hpt_tr))
            try
            h1_(:, i_tr) = interp1(f1_hpt(~isnan(f1_hpt(:,i_tr)), i_tr), ...
                h1_hpt(~isnan(f1_hpt(:,i_tr)), i_tr), f_all);
            catch er__
                warning('a');
            end
        end
        h_means_2{1, 2}(:, i_block) = nanmean(h1_,2);

        % for no-cue:
        [f1_hpt,h1_hpt] = compute_psd(b_data, Fs, nocue_tr & hpt_tr, 0); %all of movement
        h1_ = nan(length(f_all), size(f1_hpt,2));
        for i_tr = 1:length(tr_inds(nocue_tr & hpt_tr))
            try
            h1_(:, i_tr) = interp1(f1_hpt(~isnan(f1_hpt(:,i_tr)), i_tr), ...
                h1_hpt(~isnan(f1_hpt(:,i_tr)), i_tr), f_all);
            catch er__
                warning('a');
            end
        end
        h_means_2{1, 3}(:, i_block) = nanmean(h1_,2);

        % for each trial without discrimination:
        [f1_all, h1_all] = compute_psd(b_data, Fs, boolean(ones(size(tr_inds))), 0); %all movement
        temp_f = f_all >= 3 & f_all <= 7;
        for i_tr = 1:length(tr_inds)
            try
                temp_h = interp1(f1_all(~isnan(f1_all(:,i_tr)), i_tr), ...
                    h1_all(~isnan(f1_all(:,i_tr)), i_tr), f_all);
            catch er__
                warning('interpolation error')
            end
            try
                h_metric_trials(total_trial_k, 1) = nanmean(temp_h(temp_f));
                h_metric_trials(total_trial_k, 2) = pred_tr(i_tr);
                h_metric_trials(total_trial_k, 3) = pt(i_tr);
                h_metric_trials(total_trial_k, 4) = hpt_tr(i_tr);
                h_metric_trials(total_trial_k, 5) = nocue_tr(i_tr);
                h_metric_trials(total_trial_k, 6) = total_trial_k;
                h_metric_trials(total_trial_k, 7) = i_block;
                h_metric_trials(total_trial_k, 8) = reward_tr(i_tr);
                
                t_offset_temp = b_data.t(b_data.k_split(i_tr), i_tr);
                t_shifted_temp = b_data.t(:, i_tr) - t_offset_temp;
                kin_metric_trials(total_trial_k, 1:length(t_shifted_temp), 1) = t_shifted_temp;
                if flip_tr(i_tr)
                    kin_metric_trials(total_trial_k, 1:length(t_shifted_temp), 2) =...
                        -b_data.y(:, i_tr);
                    kin_metric_trials(total_trial_k, 1:length(t_shifted_temp), 3) =...
                        -b_data.vy(:, i_tr);
                    kin_metric_trials(total_trial_k, 1:length(t_shifted_temp), 4) =...
                        -b_data.ay(:, i_tr);
                else
                    kin_metric_trials(total_trial_k, 1:length(t_shifted_temp), 2) =...
                        b_data.y(:, i_tr);
                    kin_metric_trials(total_trial_k, 1:length(t_shifted_temp), 3) =...
                        b_data.vy(:, i_tr);
                    kin_metric_trials(total_trial_k, 1:length(t_shifted_temp), 4) =...
                        b_data.ay(:, i_tr);
                end
                kin_metric_trials(total_trial_k, 1:length(t_shifted_temp), 5) = b_data.ma(:, i_tr);
                
                kin_metric_trials(total_trial_k, 1:length(t_shifted_temp), 6) =...
                    b_data.x(:, i_tr);
                kin_metric_trials(total_trial_k, 1:length(t_shifted_temp), 7) =...
                    b_data.vx(:, i_tr);
                kin_metric_trials(total_trial_k, 1:length(t_shifted_temp), 8) =...
                    b_data.ax(:, i_tr);
                
                % measure directional error of velocity from target:
                ksplit_ = b_data.k_split(i_tr);
                x__ = b_data.x(ksplit_:end, i_tr);
                vx_ = b_data.vx(ksplit_:end, i_tr);
                if flip_tr(i_tr)
                    y__ = -b_data.y(ksplit_:end, i_tr);
                    vy_ = -b_data.vy(ksplit_:end, i_tr);
                else
                    y__ = b_data.y(ksplit_:end, i_tr);
                    vy_ = b_data.vy(ksplit_:end, i_tr);
                end
                d__ = sqrt(x__.^2 + y__.^2);
                [~, k__th] = min(abs(d__ - dist_th));
                v__dir = atan2d(vy_(k__th), vx_(k__th));
                v_metric_trials(total_trial_k, 1) = v__dir - TARG_DIR;
                v_metric_trials(total_trial_k, 2) = pred_tr(i_tr);
                v_metric_trials(total_trial_k, 3) = pt(i_tr);
                v_metric_trials(total_trial_k, 4) = hpt_tr(i_tr);
                v_metric_trials(total_trial_k, 5) = nocue_tr(i_tr);
                v_metric_trials(total_trial_k, 6) = total_trial_k;
                v_metric_trials(total_trial_k, 7) = i_block;
                v_metric_trials(total_trial_k, 8) = reward_tr(i_tr);
                
                k_split_trials(total_trial_k, 1) = ksplit_;
                
                total_trial_k = total_trial_k + 1;
                temp_h = nan(size(f_all));
                t_shifted_temp = nan;
            catch er__
                warning('data assignment error')
                total_trial_k = total_trial_k + 1;
                temp_h = nan(size(f_all));
                t_shifted_temp = nan;
            end
        end

        %% extract PSD: from low-pt and all of movement

        % for pred:
        [f1_hpt,h1_hpt] = compute_psd(b_data, Fs, pred_tr & lpt_tr, 0); %all of movement
        h1_ = nan(length(f_all), size(f1_hpt,2));
        for i_tr = 1:length(tr_inds(pred_tr & lpt_tr))
            try
            h1_(:, i_tr) = interp1(f1_hpt(~isnan(f1_hpt(:,i_tr)), i_tr), ...
                h1_hpt(~isnan(f1_hpt(:,i_tr)), i_tr), f_all);
            catch er__
                warning('a');
            end
        end
        h_means_2{2, 1}(:, i_block) = nanmean(h1_,2);

        % for unpred:
        [f1_hpt,h1_hpt] = compute_psd(b_data, Fs, unpred_tr & lpt_tr, 0); %all of movement
        h1_ = nan(length(f_all), size(f1_hpt,2));
        for i_tr = 1:length(tr_inds(unpred_tr & lpt_tr))
            try
            h1_(:, i_tr) = interp1(f1_hpt(~isnan(f1_hpt(:,i_tr)), i_tr), ...
                h1_hpt(~isnan(f1_hpt(:,i_tr)), i_tr), f_all);
            catch er__
                warning('a');
            end
        end
        h_means_2{2, 2}(:, i_block) = nanmean(h1_,2);

        % for no-cue:
        [f1_hpt,h1_hpt] = compute_psd(b_data, Fs, nocue_tr & lpt_tr, 0); %all of movement
        h1_ = nan(length(f_all), size(f1_hpt,2));
        for i_tr = 1:length(tr_inds(nocue_tr & lpt_tr))
            try
            h1_(:, i_tr) = interp1(f1_hpt(~isnan(f1_hpt(:,i_tr)), i_tr), ...
                h1_hpt(~isnan(f1_hpt(:,i_tr)), i_tr), f_all);
            catch er__
                warning('a');
            end
        end
        h_means_2{2, 3}(:, i_block) = nanmean(h1_,2);


        %% compute movement error:
%         v_mag = nan(1, size(b_data.x,2));
%         for i_tr = 1:size(b_data.x,2)
%             dist_sig = sqrt((b_data.x(:, i_tr)).^2 + (b_data.y(:, i_tr)).^2);
%             [~, k_min] = min((dist_sig - dist_th).^2);
% %             v_mag(i_tr) = sqrt((b_data.vx(k_min, i_tr)).^2 + (b_data.vy(k_min, i_tr)).^2);
%             % get the direction of the velocity.
%             
%             % compute the error between direction of the velocity and
%             % direction of target:
%             
%         end
%         v_err(i_block) = nanmean(v_mag);

    catch er_
        warning('e');
    end
end

s_data.f = f_all;
s_data.h_means_1 = h_means_1;
s_data.h_means_2 = h_means_2;
% s_data.v_err = v_err;
s_data.h_all = h_metric_trials;
s_data.kin_all = kin_metric_trials;
s_data.v_err_all = v_metric_trials;
s_data.k_split_all = k_split_trials;

function [f,h] = compute_psd(sd, F, valid_trs, which_half)
% compute the first half of signal psd
% f = nan(size(sd.x,1), size(sd.ma,2));
% h = nan(size(sd.x,1), size(sd.ma,2));

f = nan(size(sd.x,1), sum(valid_trs));
h = nan(size(sd.x,1), sum(valid_trs));

all_trs = 1:length(valid_trs);
trs_to_analyze = all_trs(valid_trs);

for i_ = 1:length(trs_to_analyze)
    if ~isnan(sd.k_split(trs_to_analyze(i_)))
        % select out the magnitude of acceleration:
        ss_temp = sd.ma(:, trs_to_analyze(i_));
        ss = ss_temp(~isnan(ss_temp));
        switch which_half
            case 1
                % first half
                [f_,h_] = simple_psd(ss(1:sd.k_split(trs_to_analyze(i_))), F);
            case 2
                % second half
                [f_,h_] = simple_psd(ss((sd.k_split(trs_to_analyze(i_))+1):end), F);
            case 0
                % entire trajectory
                [f_,h_] = simple_psd(ss, F);
            otherwise
                error('The desired half of the signal to compute PSD of is not specified.');
        end
        f(1:length(f_), i_) = f_;
        h(1:length(h_), i_) = h_;
    end
end

function dir_err = compute_velocity_error(sd, target)
% compute the directional error: difference between direction of velocity
% and target direction.



