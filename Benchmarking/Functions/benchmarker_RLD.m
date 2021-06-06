function [ RLD_time, RLD_memory, RLD_results ] = ...
    benchmarker_RLD( photon_data, combined_data, exposure_time, lite_flag )
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


% Check how many time bins are in the photon data and generate the two
% component images D_0 and D_1. Also calculates the real time exposure of
% the component images
if rem(size(combined_data, 3), 2) == 0
    D_0 = sum(combined_data( :, :, 1:(size(combined_data, 3)/2)), 3);
    D_1 = sum(combined_data( :, :, ((size(combined_data, 3)/2)+1):end), 3);
    delta_t = exposure_time * (size(combined_data, 3)/2);
else
    D_0 = sum(combined_data( :, :, 1:((size(combined_data, 3)-1)/2)), 3);
    D_1 = sum(combined_data( :, :, (((size(combined_data, 3)-1)/2)+1):end), 3);
    delta_t = exposure_time * ((size(combined_data, 3)-1)/2);
end

% Calculate tau (Lifetime)
tau = (-1 * delta_t) ./ log( D_1 ./ D_0 );


% % Calculate A (pre-exponential factor)
% A = D_0 ./ (tau .* (1 - (D_1 ./ D_0)));


% End Timing
RLD_time.combined = toc(start_combined);

% Output Combined Lifetime Image
RLD_results.combined = tau;



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
    
    RLD_memory.iterative(i) = 3 * (temp_list(temp_ind).bytes / ...
        sum(temp_list(temp_ind).size, 'all'));
    
    
    
    % Start Timing
    start_iter = tic;
    
    % Build the Cumulative Counts
    cumulative_counts = cumulative_counts + photon_data(i).counts;
    
    % Check how many time bins are present and then build D_0 and D_1.
    % Further calculate the delta_t.
    if rem(size(cumulative_counts, 3), 2) == 0
        D_0 = sum(cumulative_counts( ...
            :, :, 1:(size(cumulative_counts, 3)/2)), 3);
        D_1 = sum(cumulative_counts( ...
            :, :, ((size(cumulative_counts, 3)/2)+1):end), 3);
        delta_t = exposure_time * (size(cumulative_counts, 3)/2);
    else
        D_0 = sum(cumulative_counts( ...
            :, :, 1:((size(cumulative_counts, 3)-1)/2)), 3);
        D_1 = sum(cumulative_counts( ...
            :, :, (((size(cumulative_counts, 3)-1)/2)+1):end), 3);
        delta_t = exposure_time * ((size(cumulative_counts, 3)-1)/2);
    end
    
    % Calculate tau (Lifetime)
    tau = (-1 * delta_t) ./ log( D_1 ./ D_0 );
    
    
    % % Calculate A (pre-exponential factor)
    % A = D_0 ./ (tau .* (1 - (D_1 ./ D_0)));
    
    % End Timing
    RLD_time.iterative(i) = toc(start_iter);
    
    if lite_flag == 0 
        RLD_results.iterative(i).result = tau;
    else
        if i == 1
            RLD_results.iterative(1).result = tau;
        elseif i == round(numel(photon_data)/2)
            RLD_results.mid_iter_ind = round(numel(photon_data)/2);
            RLD_results.iterative(2).result = tau;
        elseif i == numel(photon_data)
            RLD_results.iterative(3).result = tau;
        end
    end
end



%% Cleanup from Benchmarking
if ishandle(wait_box)
    close(wait_box);
end
end