function [ PCA_time, PCA_memory, PCA_results ] = ...
    benchmarker_PCA( photon_data, combined_data )
%% PCA Benchmarker
%   By: Niklas Gahm
%   2020/11/16
%
%   This code takes the read in photon data and uses it for PCA FLIM
%   Estimation. This process is benchmarked for time and memory usage. It
%   further outputs PCA's iterative results for qualitative assessment.
% 
%   This is based off of work presented in:
%   https://onlinelibrary.wiley.com/doi/full/10.1002/jbio.201600160
%
%   2020/11/16 - Started
%   2020/11/17 - Finished




%% Initialize Outputs
% Initialize Time Struct
PCA_time = struct;
PCA_time.iterative = struct;
PCA_time.iterative.time(1) = 0;
PCA_time.combined = 0;

% Initialize Memory Struct
PCA_memory = struct;
PCA_memory.iterative = struct;
PCA_memory.iterative.memory(1) = 0;
PCA_memory.combined = 0;

% Initialize Results Struct
PCA_results = struct;
PCA_results.iterative = struct;
PCA_results.iterative.result(1) = 0;
PCA_results.combined = 0;



%% Inform user to be Patient
wait_box = msgbox({'The script is now going into measurement mode.'; ... 
    'No waitbar will be shown to get accurate timing data.'}, ...
    'PCA Benchmarking');



%% Measurement Section for Combined Data
% Memory Usage Estimation
temp_list = whos;
temp_ind = 1;
for i = 1:numel(temp_list)
    if strcmp(temp_list(i).name, 'combined_data')
        temp_ind = i;
        break;
    end
end

% Memory estimate is centered on the largest matrix and the operations done
% to is. This is a rough estimate, not an exact number.
PCA_memory.combined = 3*(temp_list(temp_ind).bytes);


% Start Timing
start_combined = tic;

%%% Calculate Useful Values
img_size = size(combined_data);
num_time_bins = img_size(3);
num_pixels = numel(combined_data(:,:,1));

%%% Perform Noise Correction
normalization_factors = sqrt(squeeze(mean(mean(combined_data, 1),2)));
for i = 1:num_time_bins
    if normalization_factors(i) < 1
        % This case occurs when a time bin is very predominantely zeros. As
        % such since a normalization/division by 0 is impossible, the 
        % normalization for the bin is set to 1 to have no impact on that
        % time bin. 
        normalization_factors(i) = 1;
    end
    combined_data(:,:,i) = combined_data(:,:,i) / normalization_factors(i);
end

%%% Generate Covariance Matrix
% Initialize covariance matrix to be filled in
cov_mat = zeros(num_time_bins);
for i = 1:num_time_bins
    for j = 1:num_time_bins
        cov_mat(i,j) = sum(((combined_data(:,:,i) - ...
            mean(combined_data(:,:,i), 'all')) .* ...
            (combined_data(:,:,j) - ...
            mean(combined_data(:,:,j), 'all'))), 'all');
    end
end
cov_mat = cov_mat / num_pixels;

%%% Get Eigen Vectors and Values via SVD
[~, ~, V] = svd(cov_mat);

%%% Generate Score Values
combined_data = reshape(combined_data, num_pixels, num_time_bins);
score = V' * combined_data';    % score is actually score' here
score = score' * normalization_factors; 
PCA_results.combined = reshape(score, img_size(1), img_size(2));

% End Timing
PCA_time.combined = toc(start_combined);



%% Measurement Section for Iterative Data
cumulative_counts = photon_data(1).counts .* 0;
for i = 1:numel(photon_data)
    % Memory Usage Estimation
    temp_list = whos;
    temp_ind = 1;
    for j = 1:numel(temp_list)
        if strcmp(temp_list(j).name, 'photon_data')
            temp_ind = j;
            break;
        end
    end
    
    PCA_memory.iterative.memory(i) = 3 * (temp_list(temp_ind).bytes / ...
        sum(temp_list(temp_ind).size, 'all'));
    
    % Start Timing
    start_iter = tic;
    
    %%% Build Accumulation Matrix
    cumulative_counts = cumulative_counts + photon_data(i).counts;
    iter_counts = cumulative_counts;
    
    %%% Calculate Useful Values
    img_size = size(cumulative_counts);
    num_time_bins = img_size(3);
    num_pixels = numel(cumulative_counts(:,:,1));

    %%% Perform Noise Correction
    normalization_factors = sqrt(squeeze(mean(mean(iter_counts, 1),2)));
    for j = 1:num_time_bins
        if normalization_factors(j) < 1
            % This case occurs when a time bin is very predominantely
            % zeros. As such since a normalization/division by 0 is
            % impossible, the normalization for the bin is set to 1 to
            % have no impact on that time bin.
            normalization_factors(j) = 1;
        end
        iter_counts(:,:,j) = iter_counts(:,:,j) / normalization_factors(j);
    end
    
    %%% Generate Covariance Matrix
    % Initialize covariance matrix to be filled in
    cov_mat = zeros(num_time_bins);
    for k = 1:num_time_bins
        for j = 1:num_time_bins
            cov_mat(k,j) = sum(((iter_counts(:,:,k) - ...
                mean(iter_counts(:,:,k), 'all')) .* ...
                (iter_counts(:,:,j) - ...
                mean(iter_counts(:,:,j), 'all'))), 'all');
        end
    end
    cov_mat = cov_mat / num_pixels;
    
    %%% Get Eigen Vectors and Values via SVD
    [~, ~, V] = svd(cov_mat);
    
    %%% Generate Score Values
    iter_counts = reshape(iter_counts, num_pixels, num_time_bins);
    score = V' * iter_counts';    % score is actually score' here
    score = score' * normalization_factors;
    PCA_results.iterative(i).result = ...
        reshape(score, img_size(1), img_size(2));
    
    
    % End Timing
    PCA_time.iterative.time(i) = toc(start_iter);
end



%% Cleanup from Benchmarking
if ishandle(wait_box)
    close(wait_box);
end


end