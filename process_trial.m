function [t_data, varargout] = process_trial(Data, i_tr)

% function p_data = process_trial(Data, i_tr)
%
% INPUTS:
%   - Data - data structure saved from experiment
%   - i_tr - the trial number to be processed
%
% OUTPUTS:
%   - t_data - trial data structure with fields t,x,y,vx,vy,ax,ay,am,k_split
%   - errors - cell array with any error messages
%
% Process all trials individually, including:
%   - Align to the template (rotate trajectory as if all start positions
%   were left horizonal).
%   - isolate trajectory from the time that the main movement started
%   - Compute velocity and acceleration; compute norm of those.
%   - Find the time between movement onset and when the target appeared.
%   - Find the time between stimulus onset and when movement began.
%   - Find the time between stimulus onset and target onset.
%   - Find the time between the movement reaching the center and target
%   appearance.
%   - Find the directional error of the second sub-movement relative the
%   the target position.
%   - Compute the PSD of the acceleration amplitude


mm_per_pix = 1;
screen_dims = [1920 1080]*mm_per_pix;
screen_center = screen_dims/2;
rot_start_point = [-350, 0]*mm_per_pix;
start_buffer = 5;
targ_diam = 15;
center_th = 100;
targ_length = 350;

targ_set = [5 6 7 8];
rotation_angle_targ = [135 45 -45 -135];


filt_params = [3 15];

err_set = {};
try
    x_temp = Data.Kinematics{i_tr}(:,2) - screen_center(1);
    y_temp = Data.Kinematics{i_tr}(:,3) - screen_center(2);
    t = Data.Kinematics{i_tr}(:,1) - Data.Kinematics{i_tr}(1,1);
    T = 1/mean(diff(t));
    H = 1/T;

    ra_temp = rotation_angle_targ(targ_set == Data.params.trial_home_numbers(i_tr));
    R = [cosd(ra_temp) -sind(ra_temp);...
        sind(ra_temp) cosd(ra_temp)];

    p_rot = R*[x_temp'; y_temp'];

    k = 1:length(t);
    
    % find relative time of target appearance:
    t_targ = Data.time_targ_disp(i_tr) - Data.Kinematics{i_tr}(1,1);
    [~, k_targ_] = min((t_targ - t).^2);

    % find movement start:
    near_start = sqrt((p_rot(1,:) - rot_start_point(1)).^2 + ...
        (p_rot(2,:) - rot_start_point(2)).^2);
    near_start_trig_ = near_start < targ_diam;
    near_start_trig = [0, diff(near_start_trig_)];
    k_start_ = k(near_start_trig < 0);
    k_start = k_start_(1) - start_buffer; 
    % find first time the subject leaves the start position (and take a bit before then);

    t__ = t(k_start:end); t_start_ = t__;
    x__ = p_rot(1, k_start:end); x_start_ = x__;
    y__ = p_rot(2, k_start:end); y_start_ = y__;
    k__ = 1:length(t__); k_start_ = k__;

    % find "completion" of sub-movement 1:
    near_center = sqrt(x__.^2 + y__.^2);
    near_center_trig_ = near_center < center_th;
    near_center_trig = [0, diff(near_center_trig_)];
    k_center_1_ = k__(near_center_trig > 0);
    k_center_2_ = k__(near_center_trig < 0);
    k_center_1 = k_center_1_(1);
    k_center_2 = k_center_2_(1);% find the first time (after movement onset) that subject enters
    % and leaves center position)
    [~, k_center_] = min(near_center(k_center_1:k_center_2));
    k_center = k_center_ + k_center_1;

    t__ = t__(k_center:end); t_center_ = t__;
    x__ = x__(k_center:end); x_center_ = x__;
    y__ = y__(k_center:end); y_center_ = y__;
    k__ = 1:length(t__);

    % find movement end (finish)
    position_magnitude = sqrt(x__.^2 + y__.^2);
    mag_trig_ = position_magnitude > targ_length;
    mag_trig = [0 diff(mag_trig_)];
    k_finish_ = k__(mag_trig > 0);
    if ~isempty(k_finish_)
        k_finish = k_finish_(1);
    else
        k_finish = k_center;
    end

    t_finish_ = t__(1:k_finish);
    x_finish_ = x__(1:k_finish);
    y_finish_ = y__(1:k_finish);

    k_targ = k_targ_ - (k_start - 1);
    
    inds_mov = k_start:(k_start + k_center + k_finish); % note: k_center is
%     relative to k_start, and k_finish is relative to k_center.
    t_ = t(inds_mov);
    x_ = p_rot(1, inds_mov);
    y_ = p_rot(2, inds_mov);

    vx_ = [0 diff(x_)]./T;
    vy_ = [0 diff(y_)]./T;

    vx = sgolayfilt(vx_, filt_params(1), filt_params(2));
    vy = sgolayfilt(vy_, filt_params(1), filt_params(2));

    ax_ = [0 diff(vx)]./T;
    ay_ = [0 diff(vy)]./T;

    ax = sgolayfilt(ax_, filt_params(1), filt_params(2));
    ay = sgolayfilt(ay_, filt_params(1), filt_params(2));

    am = sqrt(ax.^2 + ay.^2);        

catch err_tr
    warning(['Trial ', num2str(i_tr), ' failed.']);
    err_set{length(err_set) + 1} = err_tr;
end

t_data.t = t_;
t_data.x = x_;
t_data.y = y_;
t_data.vx = vx;
t_data.vy = vy;
t_data.ax = ax;
t_data.ay = ay;
t_data.am = am;
t_data.k_split = k_center;
t_data.k_targ = k_targ;

varargout{1} = err_set;