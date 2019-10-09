function [b_data, varargout] = process_block(Data)
% function p_data = process_block(Data)
%
% Process all trials in a block, including:
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

t_all = nan(300, length(Data.pPT));
x_all = nan(300, length(Data.pPT));
y_all = nan(300, length(Data.pPT));
vx_all = nan(300, length(Data.pPT));
vy_all = nan(300, length(Data.pPT));
ax_all = nan(300, length(Data.pPT));
ay_all = nan(300, length(Data.pPT));
ma_all = nan(300, length(Data.pPT));
k_split = nan(1, length(Data.pPT));
k_targ = nan(1, length(Data.pPT));

b_errs = {};
t_errs = cell(1, length(Data.pPT));
for i_tr = 1:length(Data.pPT)
    try        
        [t_data, t_err] = ...
            process_trial(Data, i_tr);
        
        t_all(1:length(t_data.t), i_tr) = t_data.t;
        x_all(1:length(t_data.t), i_tr) = t_data.x;
        y_all(1:length(t_data.t), i_tr) = t_data.y;
        vx_all(1:length(t_data.t), i_tr) = t_data.vx;
        vy_all(1:length(t_data.t), i_tr) = t_data.vy;
        ax_all(1:length(t_data.t), i_tr) = t_data.ax;
        ay_all(1:length(t_data.t), i_tr) = t_data.ay;
        ma_all(1:length(t_data.t), i_tr) = t_data.am;
        k_split(i_tr) = t_data.k_split;
        k_targ(i_tr) = t_data.k_targ;
        
        t_errs{i_tr} = t_err;
    catch err_tr
        warning(['Trial ', num2str(i_tr), ' failed.']);
        b_errs{length(b_errs) + 1} = err_tr;
    end
    
end

b_data.t = t_all;
b_data.x = x_all;
b_data.y = y_all;
b_data.vx = vx_all;
b_data.vy = vy_all;
b_data.ax = ax_all;
b_data.ay = ay_all;
b_data.ma = ma_all;
b_data.k_split = k_split;
b_data.k_targ = k_targ;

varargout{1} = {t_errs, b_errs};