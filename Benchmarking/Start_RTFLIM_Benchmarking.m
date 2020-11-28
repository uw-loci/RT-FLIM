%% Start Run Time FLIM Benchmarking
%   This is the run file to start the Run Time FLIM benchmarking framework.
% 
%     Copyright (C) 2020 Niklas Gahm
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%   2020/11/19 - Started 
% 
%   To-Do:
%       - Expand to compare an array of time bin sizes
%       - Expand to compare multiple benchmarking files.




%% Setup the Workspace
clear;
format longe; 



%% User Variables

% If Data Visualizer function is on (1) or off (0)
visualizer_flag = 1;

% If it only keeps representative results. 
lite_flag = 1;

% Order of the Dimensions in the Data
data_order = 'TXYS';

% How Large a Time Bin is
% time_bin_size = [1,2,4,8,16,32];
time_bin_size = [16;


%% Benchmark Files
% To add more files, add another element to the benchmark_files cell array 
% with the full file path and file name.

% benchmark_files = {...
%     'D:\LOCI\FLIM\Runtime FLIM Benchmarking\Data_20200317\data_ch2_hig_photons.h5', ... % Low Artifact Benchmark Set High Flux
%     'D:\LOCI\FLIM\Runtime FLIM Benchmarking\Data_20200317\data_ch2_med_photons.h5', ... % Low Artifact Benchmark Set Medium Flux
%     'D:\LOCI\FLIM\Runtime FLIM Benchmarking\Data_20200317\data_ch2_low_photons.h5', ... % Low Artifact Benchmark Set Low Flux
%     'D:\LOCI\FLIM\Runtime FLIM Benchmarking\Data_20191007\data.h5'}; % Benchmark Set with Artifacts

benchmark_files = {...
    'D:\LOCI\FLIM\Runtime FLIM Benchmarking\Data_20200317\data_ch2_hig_photons.h5'};


%% Add Necessary Paths
addpath('Functions');



%% Run Benchmarking
collected_metrics = struct;
for i = 1:numel(benchmark_files)
    [~,collected_metrics(i).name_str,~] = fileparts(benchmark_files{i});
    for j = 1:numel(time_bin_size) 
        collected_metrics(i).name(j).time_bin_size = time_bin_size(j);
        collected_metrics(i).name(j).metrics = ... 
            RTFLIM_Benchmarking_Framework(benchmark_files{i}, ...
            data_order, time_bin_size(j), 0, lite_flag);
        fprintf(['\nBenchmarked ', ...
            strrep(collected_metrics(i).name_str, '_', ' '), ' at ', ...
            num2str(time_bin_size(j)), ' Wide Time Bins\n']);
    end
end



%% Visualize Results
if visualizer_flag == 1
    fprintf('\nVisualizing Results Across Benchmarks and Time Bins`\n');
    RTFLIM_collected_benchmarks_visualizer(collected_metrics);
end


%% Inform User of Completion
fprintf('\n\nAll Benchmarking is Complete.\n\n');