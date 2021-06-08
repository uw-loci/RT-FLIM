function [ Laguerre_time, Laguerre_memory, Laguerre_results ] = ...
    benchmarker_Laguerre( photon_data, combined_data, lag_degree, ...
    exposure_time, lite_flag )
%% Laguerre Benchmarker
%   By: Niklas Gahm
%   2020/11/16
%
%   This code takes the read in photon data and uses it for Laguerre FLIM
%   Estimation. This process is benchmarked for time and memory usage. It
%   further outputs Laguerre's iterative results for qualitative 
%   assessment.
% 
%   Based on: https://ieeexplore.ieee.org/document/5594996
%             https://numpy.org/doc/stable/reference/generated/numpy.polynomial.laguerre.lagfit.html
%
%   2020/11/16 - Started



%% Initialize Outputs
% Initialize Time Struct
Laguerre_time = struct;
Laguerre_time.iterative(1) = 0;
Laguerre_time.combined = 0;

% Initialize Memory Struct
Laguerre_memory = struct;
Laguerre_memory.iterative(1) = 0;
Laguerre_memory.combined = 0;

% Initialize Results Struct
Laguerre_results = struct;
Laguerre_results.iterative = struct;
Laguerre_results.iterative.result(1) = 0;
Laguerre_results.combined = 0;
Laguerre_results.mid_iter_ind = round(numel(photon_data)/2);



%% Inform user to be Patient
wait_box = msgbox({'The script is now going into measurement mode.'; ... 
    'No waitbar will be shown to get accurate timing data.'}, ...
    'Laguerre Benchmarking');



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
Laguerre_memory.combined = 4*(temp_list(temp_ind).bytes);



% Start Timing
start_combined = tic;

% Reshape Image to a 2D Pixels in Columns and time bins in Rows
pixel_data = zeros(size(combined_data, 3), numel(combined_data(:,:,1)));
for i = 1:size(combined_data,3)
    pixel_data(i,:) = ...
        reshape(combined_data(:,:,i), 1, numel(combined_data(:,:,1)));
end

% Generate Corresponding Time Matrix
bin_times = 1:1:size(combined_data,3);
bin_times = bin_times * exposure_time;

% Fit the Coefficients of the Laguerre Series
coef = py.numpy.polynomial.laguerre.lagfit( ...
    int32(bin_times), int32(pixel_data), int32(lag_degree) );
coef = double(coef);

% The coeficients themselves show the presence and differences of lifetime
% based on this, the resultant score image is the map of the first
% coefficient

% Generate Image
Laguerre_results.combined = reshape(coef(1,:), size(combined_data, 1), ...
    size(combined_data, 2));

% End Timing
Laguerre_time.combined = toc(start_combined);



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
    
    Laguerre_memory.iterative(i) = 3 * (temp_list(temp_ind).bytes / ...
        sum(temp_list(temp_ind).size, 'all'));
    
    % Start Timing
    start_iter = tic;
    
    % Build the Cumulative Counts
    cumulative_counts = cumulative_counts + photon_data(i).counts;
    
    % Collect Image Dimensions for Reshaping
    img_dim_1 = size(cumulative_counts, 1);
    img_dim_2 = size(cumulative_counts, 2);
    img_dim_3 = size(cumulative_counts, 3);
    
    % Reshape Image to a 2D Pixels in Columns and time bins in Rows
    pixel_data = zeros(img_dim_3, (img_dim_1 * img_dim_2));
    for j = 1:img_dim_3
        pixel_data(j,:) = ...
            reshape(cumulative_counts(:,:,j), 1, (img_dim_1 * img_dim_2));
    end
    
    % Generate Corresponding Time Matrix
    bin_times = 1:1:img_dim_3;
    bin_times = bin_times * exposure_time;
    
    % Fit the Coefficients of the Laguerre Series
    coef = py.numpy.polynomial.laguerre.lagfit( ...
        int32(bin_times), int32(pixel_data), int32(lag_degree) );
    coef = double(coef);
    
    % The coeficients themselves show the presence and differences of
    % lifetime based on this, the resultant score image is the map of the
    % first coefficient
    
    % Generate Image
    res_img = reshape(coef(1,:), img_dim_1, img_dim_2);
    
    
    % End Timing
    Laguerre_time.iterative(i) = toc(start_iter);
    
    if lite_flag == 0
        Laguerre_results.iterative(i).result = res_img;
    else
        if i == 1
            Laguerre_results.iterative(1).result = res_img;
        elseif i == round(numel(photon_data)/2)
            Laguerre_results.iterative(2).result = res_img;
        elseif i == numel(photon_data)
            Laguerre_results.iterative(3).result = res_img;
        end
    end
end



%% Cleanup from Benchmarking
if ishandle(wait_box)
    close(wait_box);
end
end