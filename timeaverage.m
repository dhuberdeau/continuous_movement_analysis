function [time_, avg_] = timeaverage(time_mat, data_mat, varargin)

% function [t, ta] = timeaverage(time_mat, data_mat)
%
% Compute the average signal at each time from the global min of time_mat
% to the global max of time_mat. Specify a third argument to indicate the
% minimum number of signals going into the average (to avoid having only 1
% or a few signals being called the average at the bookends of the time
% window).
%
% time_mat and data_mat should be matricies of Nxns, where N is the length
% of the signal and ns is the number of signals.

if nargin > 2
    min_samps = varargin{1};
else
    min_samps = 1;
end

t1 = min(time_mat(:));
tf = max(time_mat(:));

% find the resolution of the time matrix:
iti_set = nan(1, size(time_mat,2));
for i_sample = 1:size(time_mat,2)
    iti_set(i_sample) = nanmean(diff(time_mat(:, i_sample)));
end

mean_iti = nanmean(iti_set);
time_ = t1:mean_iti:tf;
avg_ = nan(size(time_));

for i_t = 1:length(time_)
    samps = nan(size(data_mat,2), 1);
    for i_sample = 1:size(data_mat,2)
        try
            this_signal__ = data_mat(:,i_sample);
            this_time__ = time_mat(:, i_sample);
            this_time_ = this_time__(~isnan(this_time__));
            this_signal_ = this_signal__(~isnan(this_time__));

            this_time = time_(find(time_ >= this_time_(1) & time_ <= this_time_(end)));
            this_signal = interp1(this_time_(:), this_signal_(:), this_time);

    %         this_iti = nanmean(diff(this_time));
            this_samp = this_signal(find(this_time >= (time_(i_t) - mean_iti/2) & ...
                this_time < (time_(i_t) + mean_iti/2)));
            if ~isempty(this_samp)
                samps(i_sample) = this_samp(1);
            end
        catch
%             warning('sample failed');
        end
    end
    
    avg_(i_t) = nanmean(samps);
end
