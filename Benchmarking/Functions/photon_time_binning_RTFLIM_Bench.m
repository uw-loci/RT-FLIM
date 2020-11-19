function [ photon_data, time_bin_size ] = ...
    photon_time_binning_RTFLIM_Bench( ...
    photon_data, time_bin_size )
%% Real Time FLIM Photon Time Binner
%   By: Niklas Gahm
%   2020/11/16
%
%   This code takes the benchmark metrics and visualizes them in a readily
%   human interpretable manner. It is designed to be able to handle an
%   arbitrary number of benchmarked algorithms, as long as they all conform
%   to the benchmark output standards. This version is specifically for the
%   Real Time FLIM project.
%
%   2020/11/16 - Started
%   2020/11/18 - Finished




%% Check if Bin Size is 1 or Less and if it is an Integer
if time_bin_size == 1
    return;
elseif time_bin_size < 1
    warning(['Interpolating into time bins is not supported. ', ...
        'Time Bin Size used is 1.']);
    time_bin_size = 1;
elseif time_bin_size ~= round(time_bin_size)
    time_bin_size = round(time_bin_size);
    warning(['Non-integer time bin sizes are not supported. ', ...
        'New time bin size is: ' num2str(time_bin_size)]);
end



%% Calculate Variables
num_bins = size(photon_data(1).counts, 3) / time_bin_size;
% We assume that the number of raw bins does not change between consequtive
% data sets.

% If there is a remainder send a warning that the final time bin will be
% removed.
if num_bins ~= round(num_bins)
    warning(['Data Set does not fit perfectly into binning. ', ...
        'Therefore to prevent artifacts the final partially filled ', ...
        'time bin is removed.']);
    num_bins = floor(num_bins);
end

% Initialize a temporary placeholder that is iteratively filled
temp = zeros(size(photon_data(1).counts,1), ...
    size(photon_data(1).counts,2), num_bins );


%% Iterate through all data sets and bin
num_datasets = numel(photon_data);
loader_bar = waitbar((1/num_datasets), 'Time Binning Photon Data.');

for i = 1:num_datasets
    waitbar((i/num_datasets), loader_bar);
    
    temp = temp.*0; % Clears the temp
    
    % Iterate on Time Bins
    for j = 1:num_bins
        ind_offset = (j-1) * time_bin_size;
        temp(:,:,j) = sum(photon_data(i).counts( :, :, ...
            (ind_offset + 1):(ind_offset + time_bin_size) ), 3);
    end
    
    % Overwrite photon_data.counts with time_binned version
    photon_data(i).counts = temp;
end
close(loader_bar);

end