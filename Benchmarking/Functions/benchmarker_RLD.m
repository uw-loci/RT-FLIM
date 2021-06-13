function [ RLD_time, RLD_memory, RLD_results ] = ...
    benchmarker_RLD( photon_data, combined_data, exposure_time, ...
    lite_flag, data_cleaning_flag, thresh )
%% RLD Benchmarker
%   By: Niklas Gahm
%   2020/11/16
%
%   This code takes the read in photon data and uses it for Rapid Lifetime
%   Determination (RLD). This process is benchmarked for time and memory 
%   usage. It further outputs RLD's iterative results for qualitative 
%   assessment.
% 
%   Based on: https://pubs.acs.org/doi/pdf/10.1021/ac00176a007
%
%   2021/06/03 - Started



%% Initialize Outputs
% Initialize Time Struct
RLD_time = struct;
RLD_time.iterative(1) = 0;
RLD_time.combined = 0;

% Initialize Memory Struct
RLD_memory = struct;
RLD_memory.iterative(1) = 0;
RLD_memory.combined = 0;

% Initialize Results Struct
RLD_results = struct;
RLD_results.iterative = struct;
RLD_results.iterative.result(1) = 0;
RLD_results.combined = 0;
RLD_results.mid_iter_ind = round(numel(photon_data)/2);



%% Inform user to be Patient
wait_box = msgbox({'The script is now going into measurement mode.'; ... 
    'No waitbar will be shown to get accurate timing data.'}, ...
    'RLD Benchmarking');



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
RLD_memory.combined = 3 * (temp_list(temp_ind).bytes);


% Start Timing
start_combined = tic;

% Calculate Useful Values
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

% Check how many time bins are in the photon data and generate the two
% component images D_0 and D_1. Also calculates the real time exposure of
% the component images
if rem(num_time_bins, 2) == 0
    D_0 = sum(combined_data( :, 1:(num_time_bins/2)), 2);
    D_1 = sum(combined_data( :, ((num_time_bins/2)+1):end), 2);
    delta_t = exposure_time * (num_time_bins/2);
else
    D_0 = sum(combined_data( :, 1:((num_time_bins-1)/2)), 2);
    D_1 = sum(combined_data( :, (((num_time_bins-1)/2)+1):end), 2);
    delta_t = exposure_time * ((num_time_bins-1)/2);
end

% Calculate tau (Lifetime)
tau = (-1 * delta_t) ./ log( D_1 ./ D_0 );


% % Calculate A (pre-exponential factor)
% A = D_0 ./ (tau .* (1 - (D_1 ./ D_0)));

% Output Combined Lifetime Image
tau_out = zeros((img_size(1)*img_size(2)),1);
if data_cleaning_flag == 1
    counter = 0;
    for i = 1:numel(tau_out)
        if~rows_to_remove(i)
            % This case the pixel here has some life time and needs a score
            % assigned to it
            counter = counter + 1;
            tau_out(i) = tau(counter);
        end
    end
else
    tau_out = tau;
end
RLD_results.combined = reshape(tau_out, img_size(1), img_size(2));


% End Timing
RLD_time.combined = toc(start_combined);



%% Measurement Section for Iterative Data
build_counts = photon_data(1).counts .* 0;
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
    
    RLD_memory.iterative(i) = 3 * (temp_list(temp_ind).bytes / ...
        sum(temp_list(temp_ind).size, 'all'));
    
    
    
    % Start Timing
    start_iter = tic;
    
    % Build the Cumulative Counts
    build_counts = build_counts + photon_data(i).counts;
    cumulative_counts = build_counts;
    
    % Calculate Useful Values
    img_size = size(cumulative_counts);
    num_time_bins = img_size(3);
    num_pixels = numel(cumulative_counts(:,:,1));
    
    % Reshape Combined Data
    cumulative_counts = ...
        reshape(cumulative_counts, num_pixels, num_time_bins);
    
    % Input Data Cleaning
    if data_cleaning_flag == 1
        % Find Max and remove data before max
        [~,I] = max(sum(cumulative_counts(:,1)));
        
        decdata = cumulative_counts( :, I:end );
        timechannels_data = size( decdata, 2 );
        
        if timechannels_data < num_time_bins
            decdata( 1, num_time_bins) = 0;
        elseif  timechannels_data > num_time_bins
            decdata = decdata(:, 1:num_time_bins);
        end
        
        
        % remove offset from data
        try
            data_off = ...
                mean(cumulative_counts(:, round(I/3):round(2*I/3)), 2);
        catch
            data_off = mean(cumulative_counts(:, 1), 2);
        end
        data_off = repmat(data_off, [1, size(decdata,2)]);
        decdata = decdata - data_off ;
        decdata(decdata<0) = 0 ;
        
        
        % Threshold (not peak vs peak. Used whole intensity vs max whole intensity)
        maxmax = max(sum(decdata, 2));
        rows_to_remove = any(sum(decdata,2) < (thresh * maxmax), 2);
        cumulative_counts(rows_to_remove, :) = [];
    end
    
    % Check how many time bins are present and then build D_0 and D_1.
    % Further calculate the delta_t.
    if rem(num_time_bins, 2) == 0
        D_0 = sum(cumulative_counts( :, 1:(num_time_bins/2)), 2);
        D_1 = sum(cumulative_counts( :, ((num_time_bins/2)+1):end), 2);
        delta_t = exposure_time * (num_time_bins/2);
    else
        D_0 = sum(cumulative_counts( :, 1:((num_time_bins-1)/2)), 2);
        D_1 = sum(cumulative_counts( :, (((num_time_bins-1)/2)+1):end), 2);
        delta_t = exposure_time * ((num_time_bins-1)/2);
    end
    
    % Calculate tau (Lifetime)
    tau = (-1 * delta_t) ./ log( D_1 ./ D_0 );
    
    
    % % Calculate A (pre-exponential factor)
    % A = D_0 ./ (tau .* (1 - (D_1 ./ D_0)));
    
    % Output Combined Lifetime Image
    tau_out = zeros((img_size(1)*img_size(2)),1);
    if data_cleaning_flag == 1
        counter = 0;
        for j = 1:numel(tau_out)
            if~rows_to_remove(j)
                % This case the pixel here has some life time and needs a score
                % assigned to it
                counter = counter + 1;
                tau_out(i) = tau(counter);
            end
        end
    else
        tau_out = tau;
    end
    tau_out = reshape(tau_out, img_size(1), img_size(2));
    
    if lite_flag == 0
        RLD_results.iterative(i).result = tau_out;
    else
        if i == 1
            RLD_results.iterative(1).result = tau_out;
        elseif i == round(numel(photon_data)/2)
            RLD_results.iterative(2).result = tau_out;
        elseif i == numel(photon_data)
            RLD_results.iterative(3).result = tau_out;
        end
    end
    
    % End Timing
    RLD_time.iterative(i) = toc(start_iter);
end



%% Cleanup from Benchmarking
if ishandle(wait_box)
    close(wait_box);
end
end