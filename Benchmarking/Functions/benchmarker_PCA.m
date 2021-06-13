function [ PCA_time, PCA_memory, PCA_results ] = ...
    benchmarker_PCA( ...
    photon_data, combined_data, lite_flag, data_cleaning_flag, thresh )
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
PCA_time.iterative(1) = 0;
PCA_time.combined = 0;

% Initialize Memory Struct
PCA_memory = struct;
PCA_memory.iterative(1) = 0;
PCA_memory.combined = 0;

% Initialize Results Struct
PCA_results = struct;
PCA_results.iterative = struct;
PCA_results.iterative.result(1) = 0;
PCA_results.combined = 0;
PCA_results.mid_iter_ind = round(numel(photon_data)/2);



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

% Reshape Combined Data
combined_data = reshape(combined_data, num_pixels, num_time_bins);

% Input Data Cleaning
if data_cleaning_flag == 1
    % Find Max and remove data before max
    [~,I] = max(sum(combined_data(:,1)));
    
    decdata = combined_data( :, I:end );
    timechannels_data = size( decdata, 2 );
    
    if timechannels_data < num_time_bins
        decdata( 1, num_time_bins) = 0;
    elseif  timechannels_data > num_time_bins
        decdata = decdata(:, 1:num_time_bins);
    end
    
    
    % remove offset from data
    try
        data_off = mean(combined_data(:, round(I/3):round(2*I/3)), 2);
    catch
        data_off = mean(combined_data(:, 1), 2);
    end
    data_off = repmat(data_off, [1, size(decdata,2)]);
    decdata = decdata - data_off ;
    decdata(decdata<0) = 0 ;
    
    
    % Threshold (not peak vs peak. Used whole intensity vs max whole intensity)
    maxmax = max(sum(decdata, 2));
    rows_to_remove = any(sum(decdata,2) < (thresh * maxmax), 2);
    combined_data(rows_to_remove, :) = [];
    
end

%%% Perform Noise Correction
normalization_factors = sqrt(squeeze(mean(combined_data, 1)));
for i = 1:num_time_bins
    if normalization_factors(i) < 1
        % This case occurs when a time bin is very predominantely zeros. As
        % such since a normalization/division by 0 is impossible, the 
        % normalization for the bin is set to 1 to have no impact on that
        % time bin. 
        normalization_factors(i) = 1;
    end
    combined_data(:,i) = combined_data(:,i) / normalization_factors(i);
end

%%% Generate Covariance Matrix
% Initialize covariance matrix to be filled in
cov_mat = zeros(num_time_bins);
for i = 1:num_time_bins
    for j = 1:num_time_bins
        cov_mat(i,j) = sum(((combined_data(:, i) - ...
            mean(combined_data(:, i))) .* ...
            (combined_data(:, j) - ...
            mean(combined_data(:, j)))), 'all');
    end
end
cov_mat = cov_mat / num_pixels;

%%% Get Eigen Vectors and Values via SVD
[~, ~, V] = svd(cov_mat);

%%% Generate Score Values
score = V' * combined_data';    % score is actually score' here
score = score' * normalization_factors'; 

comb_score = zeros((img_size(1)*img_size(2)),1);
if data_cleaning_flag == 1
    counter = 0;
    for i = 1:numel(comb_score)
        if~rows_to_remove(i)
            % This case the pixel here has some life time and needs a score
            % assigned to it
            counter = counter + 1;
            comb_score(i) = score(counter);
        end
    end
else
    comb_score = score;
end
PCA_results.combined = reshape(comb_score, img_size(1), img_size(2));

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
    
    PCA_memory.iterative(i) = 3 * (temp_list(temp_ind).bytes / ...
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
    
    % Reshape iter_counts to a vector
    iter_counts = reshape(iter_counts, num_pixels, num_time_bins);
    
    % Input Data Cleaning
    if data_cleaning_flag == 1
        % Find Max and remove data before max
        [~,I] = max(sum(iter_counts(:,1)));
        
        decdata = iter_counts( :, I:end );
        timechannels_data = size( decdata, 2 );
        
        if timechannels_data < num_time_bins
            decdata( 1, num_time_bins) = 0;
        elseif  timechannels_data > num_time_bins
            decdata = decdata(:, 1:num_time_bins);
        end
        
        
        % remove offset from data
        try
            data_off = mean(iter_counts(:, round(I/3):round(2*I/3)), 2);
        catch
            data_off = mean(iter_counts(:, 1), 2);
        end
        data_off = repmat(data_off, [1, size(decdata,2)]);
        decdata = decdata - data_off ;
        decdata(decdata<0) = 0 ;
        
        
        % Threshold (not peak vs peak. Used whole intensity vs max whole intensity)
        maxmax = max(sum(decdata, 2));
        rows_to_remove = any(sum(decdata,2) < (thresh * maxmax), 2);
        iter_counts(rows_to_remove, :) = [];
    end
    
    %%% Perform Noise Correction
    normalization_factors = sqrt(squeeze(mean(iter_counts, 1)));
    for j = 1:num_time_bins
        if normalization_factors(j) < 1
            % This case occurs when a time bin is very predominantely
            % zeros. As such since a normalization/division by 0 is
            % impossible, the normalization for the bin is set to 1 to
            % have no impact on that time bin.
            normalization_factors(j) = 1;
        end
        iter_counts(:,j) = iter_counts(:,j) / normalization_factors(j);
    end
    
    %%% Generate Covariance Matrix
    % Initialize covariance matrix to be filled in
    cov_mat = zeros(num_time_bins);
    for k = 1:num_time_bins
        for j = 1:num_time_bins
            cov_mat(k,j) = sum(((iter_counts(:,k) - ...
                mean(iter_counts(:,k))) .* ...
                (iter_counts(:,j) - ...
                mean(iter_counts(:,j)))), 'all');
        end
    end
    cov_mat = cov_mat / num_pixels;
    
    %%% Get Eigen Vectors and Values via SVD
    [~, ~, V] = svd(cov_mat);
    
    %%% Generate Score Values
    score = V' * iter_counts';    % score is actually score' here
    score = score' * normalization_factors';
    
    % Finish Data Cleaning
    comb_score = zeros((img_size(1)*img_size(2)),1);
    if data_cleaning_flag == 1
        counter = 0;
        for j = 1:numel(comb_score)
            if~rows_to_remove(j)
                % This case the pixel here has some life time and needs a 
                % score assigned to it
                counter = counter + 1;
                comb_score(j) = score(counter);
            end
        end
    else
        comb_score = score;
    end
    
    if lite_flag == 1
        if i == 1
            PCA_results.iterative(1).result = ...
                reshape(comb_score, img_size(1), img_size(2));
        elseif i == PCA_results.mid_iter_ind
            PCA_results.iterative(2).result = ...
                reshape(comb_score, img_size(1), img_size(2));
        elseif i == numel(photon_data)
            PCA_results.iterative(3).result = ...
                reshape(comb_score, img_size(1), img_size(2));
        end
    else
        PCA_results.iterative(i).result = ...
            reshape(comb_score, img_size(1), img_size(2));
    end
    
    % End Timing
    PCA_time.iterative(i) = toc(start_iter);
end



%% Cleanup from Benchmarking
if ishandle(wait_box)
    close(wait_box);
end
end