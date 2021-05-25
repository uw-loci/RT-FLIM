function [ Laguerre_time, Laguerre_memory, Laguerre_results ] = ...
    benchmarker_Laguerre( photon_data, combined_data, lite_flag )
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



im_size = size(combined_data, 1);
b = construct_laguerre_bases( im_size );
c = regress_classify(b, combined_data);


score = b*c; 










% End Timing
Laguerre_time.combined = toc(start_combined);













%% Measurement Section for Iterative Data

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
    
    
end






end