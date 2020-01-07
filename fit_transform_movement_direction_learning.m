function [t_main, y_out, v_out, a_out, ma_out, cat_out] = fit_transform_movement_direction_learning(t, y, v, a, ma, cat)

% function f = fit_transform_movement_direction(m)
%
% Classify the direction of movement of each continuous movement and flip it
% if necessary. All movements should be in the "up" direction. Return the
% kinematic variables with those trials needing to be flipped flipped.
%
% Outputs have same order as inputs; category columns are:
% 1) boolean: original category list
% 2) boolean: classified with original category list (up target; not
%    original category)
% 3) boolean: classified as opposite original category list (down target)
% 4) boolean: classified as intermediate 
%
% David Huberdeau
% October 11, 2019

N_COMPONENTS = 4;
N_CLUSTERS = 3;

Fs = 1/60;

t_main = -.5:Fs:.5;

y_interp = nan(size(y,1), length(t_main));
v_interp = nan(size(v,1), length(t_main));
a_interp = nan(size(a,1), length(t_main));
m_interp = nan(size(ma,1), length(t_main));

for i_tr = 1:size(y,1)
    try 
        t_temp = t(i_tr, :)';
        y_temp = y(i_tr, :)';
        assert(sum(~isnan(t_temp)) == sum(~isnan(y_temp)), 'Must be equal')
        t_temp = t_temp(~isnan(t_temp));
        y_temp = y_temp(~isnan(y_temp));        
        y_interp(i_tr, :) = interp1(t_temp, y_temp, t_main');
    catch err__
        warning('interpolation error');
    end
end


for i_tr = 1:size(v,1)
    try 
        t_temp = t(i_tr, :)';
        y_temp = v(i_tr, :)';
        assert(sum(~isnan(t_temp)) == sum(~isnan(y_temp)), 'Must be equal')
        t_temp = t_temp(~isnan(t_temp));
        y_temp = y_temp(~isnan(y_temp));        
        v_interp(i_tr, :) = interp1(t_temp, y_temp, t_main');
    catch err__
        warning('interpolation error');
    end
end

for i_tr = 1:size(a,1)
    try 
        t_temp = t(i_tr, :)';
        y_temp = a(i_tr, :)';
        assert(sum(~isnan(t_temp)) == sum(~isnan(y_temp)), 'Must be equal')
        t_temp = t_temp(~isnan(t_temp));
        y_temp = y_temp(~isnan(y_temp));        
        a_interp(i_tr, :) = interp1(t_temp, y_temp, t_main');
    catch err__
        warning('interpolation error');
    end
end

for i_tr = 1:size(ma,1)
    try 
        t_temp = t(i_tr, :)';
        y_temp = ma(i_tr, :)';
        assert(sum(~isnan(t_temp)) == sum(~isnan(y_temp)), 'Must be equal')
        t_temp = t_temp(~isnan(t_temp));
        y_temp = y_temp(~isnan(y_temp));        
        m_interp(i_tr, :) = interp1(t_temp, y_temp, t_main');
    catch err__
        warning('interpolation error');
    end
end

% decompose all trajectories into lower dimensional space:
[~, y_score] = pca(y_interp);

% determine k-means clusters of trajectories:
idk = kmeans(y_score(:, 1:N_COMPONENTS), N_CLUSTERS);
uid_ = unique(idk);
uid = uid_(~isnan(uid_)); % a short list of the cluster category labels

% find the cluster that corresponds most closesly to category 1 trials
% (which are the "best" trials - HPT with predictive cues on last block).
cat0_id = mode(idk(cat == 1));

% find the cluster furthest from the cat0 cluster. This will be called the
% 'opposite trajectory' that needs to be inverted. The remaining cluster
% will be considered an intermediate set of trajectories (perhaps trials that went down the
% middle initially) 

uid_remain = setdiff(uid, cat0_id); %the cluster categories that aren't cat0's
uid_dists = nan(1, length(uid_remain));
for k_id = 1:length(uid_remain)
    uid_dists(k_id) = norm(nanmean(y_score(idk == cat0_id, 1:N_COMPONENTS), 1) - ...
    nanmean(y_score(idk == uid_remain(k_id), 1:N_COMPONENTS), 1));
end
[~, ind_uid_max_dist] = max(uid_dists);
down_id = uid_remain(ind_uid_max_dist); % the cluster id furthest from cat0 - flip these trajectories.
int_id = setdiff(uid, [cat0_id, down_id]); % the remaining cluster.

y_out = y_interp;
v_out = v_interp;
a_out = a_interp;
ma_out = m_interp;
cat_out = boolean([cat == 1,...
    idk == cat0_id & ~(cat == 1),...
    idk == down_id,...
    idk == int_id]);
% y_0 = y_interp(cat == 0, :);
% v_0 = v_interp(cat == 0, :);
% a_0 = a_interp(cat == 0, :);
% ma_0 = m_interp(cat == 0, :);
% kin_0 = {y_0, v_0, a_0, ma_0};
% 
% y_up = y_interp(idk == cat0_id & ~(cat == 0), :);
% v_up = v_interp(idk == cat0_id & ~(cat == 0), :);
% a_up = a_interp(idk == cat0_id & ~(cat == 0), :);
% ma_up = m_interp(idk == cat0_id & ~(cat == 0), :);
% kin_up = {y_up, v_up, a_up, ma_up};
% 
% y_down = y_interp(idk == down_id, :);
% v_down = v_interp(idk == down_id, :);
% a_down = a_interp(idk == down_id, :);
% ma_down = m_interp(idk == down_id, :);
% kin_down = {y_down, v_down, a_down, ma_down};
% 
% y_int = y_interp(idk == int_id, :);
% v_int = v_interp(idk == int_id, :);
% a_int = a_interp(idk == int_id, :);
% ma_int = m_interp(idk == int_id, :);
% kin_int = {y_int, v_int, a_int, ma_int};

%% optional plot the trajectories of different groups:
figure; hold on
plot(t_main, y_interp(cat == 1, :)', 'k-');
plot(t_main, y_interp(idk == cat0_id & ~(cat == 1), :)', 'g-');
plot(t_main, y_interp(idk == int_id, :)', 'Color', [.5 .5 .5]);
plot(t_main, y_interp(idk == down_id, :)', 'r-');
