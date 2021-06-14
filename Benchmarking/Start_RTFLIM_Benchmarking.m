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
% 
% 
%   This was written in Matlab R2020B, and is known to use functionality
%   that is unavailable prior to Matlab R2018B.
%
%   2020/11/19 - Started 
% 
%   To-Do:
%       - Rework Visualizer




%% Setup the Workspace
clear;
format longe; 
try 
    pyenv('Version', '3.8');
    % Change '3.8' to whatever version of python is installed on your
    % computer. Most erors arrising from this will be caused by a mismatch
    % between Matlab version and supported python version. See https://www.mathworks.com/content/dam/mathworks/mathworks-dot-com/support/sysreq/files/python-compatibility.pdf
    % for what version works with each other. 
catch
    fprintf('\n\nRunning Python Version:\n');
    pyenv
end



%% User Variables

% If Data Visualizer function is on (1) or off (0)
visualizer_flag = 0;

% If it only keeps representative results. 
lite_flag = 1;

% If basic cleanining of the inpute data should be performed
data_cleaning_flag = 1;
thresh = 0.3;               % threshold [0,1]

% Number of Repetitions to be done to Mitigate Computer Variance
num_reps = 15;

% Order of the Dimensions in the Data
data_order = 'TXYS';

% How Large a Time Bin is [number of individual time gates in one time bin]
time_bin_size = [8,16,32,64];
% time_bin_size = [16];

% Exposure Time of Each Component Image [ns]
exposure_time = 0.040; % 40 ps for LOCI Systems

% The Degree of the Laguerre Series to be Used
lag_degree = 9;



%% Benchmark Files
% To add more files, add another element to the benchmark_files cell array 
% with the full file path and file name.

% benchmark_files = {...
%     'D:\LOCI\RT-FLIM\Runtime FLIM Benchmarking\Data_20200317\data_ch2_hig_photons.h5', ... % Low Artifact Benchmark Set High Flux
%     'D:\LOCI\RT-FLIM\Runtime FLIM Benchmarking\Data_20200317\data_ch2_med_photons.h5', ... % Low Artifact Benchmark Set Medium Flux
%     'D:\LOCI\RT-FLIM\Runtime FLIM Benchmarking\Data_20200317\data_ch2_low_photons.h5', ... % Low Artifact Benchmark Set Low Flux
%     'D:\LOCI\RT-FLIM\Runtime FLIM Benchmarking\Data_20191007\data.h5'}; % Benchmark Set with Artifacts

benchmark_files = {...
    'D:\LOCI\RT-FLIM\Runtime FLIM Benchmarking\Data_20200317\data_ch2_hig_photons.h5', ... % Low Artifact Benchmark Set High Flux
    'D:\LOCI\RT-FLIM\Runtime FLIM Benchmarking\Data_20200317\data_ch2_med_photons.h5', ... % Low Artifact Benchmark Set Medium Flux
    'D:\LOCI\RT-FLIM\Runtime FLIM Benchmarking\Data_20200317\data_ch2_low_photons.h5'}; % Low Artifact Benchmark Set Low Flux

% benchmark_files = {...
%     'D:\LOCI\RT-FLIM\Runtime FLIM Benchmarking\Data_20200317\data_ch2_hig_photons.h5'};


% benchmark_files_processed_FLIM = { ...
%     'D:\LOCI\RT-FLIM\Runtime FLIM Benchmarking\Data_20200317\Processed\data_ch2_hig_photons.tif', ... % Low Artifact Benchmark Set High Flux
%     'D:\LOCI\RT-FLIM\Runtime FLIM Benchmarking\Data_20200317\Processed\data_ch2_med_photons.tif', ... % Low Artifact Benchmark Set Medium Flux
%     'D:\LOCI\RT-FLIM\Runtime FLIM Benchmarking\Data_20200317\Processed\data_ch2_low_photons.tif', ... % Low Artifact Benchmark Set Low Flux
%     ''};  % Benchmark Set with Artifacts

benchmark_files_processed_FLIM = { ...
    'D:\LOCI\RT-FLIM\Runtime FLIM Benchmarking\Data_20200317\Processed\data_ch2_hig_photons.tif', ... % Low Artifact Benchmark Set High Flux
    'D:\LOCI\RT-FLIM\Runtime FLIM Benchmarking\Data_20200317\Processed\data_ch2_med_photons.tif', ... % Low Artifact Benchmark Set Medium Flux
    'D:\LOCI\RT-FLIM\Runtime FLIM Benchmarking\Data_20200317\Processed\data_ch2_low_photons.tif'}; % Low Artifact Benchmark Set Low Flux

% benchmark_files_processed_FLIM = { ...
%     'D:\LOCI\RT-FLIM\Runtime FLIM Benchmarking\Data_20200317\Processed\data_ch2_hig_photons.tif'};



%% System Variables

% Names of the Methods In order they are being benchmarked
method_names = {'NC-PCA', 'Phasor', 'LaGuerre', 'RLD'};

% How many methods are being tested by the framework
num_methods = numel(method_names);


%% Add Necessary Paths
addpath('Functions');
hpath = pwd;
spath = [hpath '\Benchmark_Results_' date];
mkdir(spath);


%% Run Benchmarking
% Initialize Storage Matrix
combined_time = zeros(num_methods, numel(time_bin_size), ...
    numel(benchmark_files), num_reps);
iterative_time = cell(num_methods, numel(time_bin_size), ...
    numel(benchmark_files), num_reps);

combined_memory = combined_time;
iterative_memory = iterative_time; 


for i = 1:num_reps
    collected_metrics = struct;
    for j = 1:numel(benchmark_files)
        [~, name_str, ~]  = fileparts(benchmark_files{j});
        for k = 1:numel(time_bin_size)
            
            temp_metrics = ...
                RTFLIM_Benchmarking_Framework(benchmark_files{j}, ...
                data_order, time_bin_size(k), exposure_time, ...
                lag_degree, 0, lite_flag, data_cleaning_flag, thresh);
            
            % Split Out and the Components 
            results = cell(1,num_methods);
            
            for m = 1:num_methods
                results{m} = getfield(temp_metrics, {m}, 'results');
                
                combined_time(m, k, j, i) = temp_metrics(m).time.combined;
                combined_memory(m, k, j, i) = ...
                    temp_metrics(m).memory.combined;
                
                iterative_time{m, k, j, i} = ...
                    temp_metrics(m).time.iterative;
                iterative_memory{m, k, j, i} = ...
                    temp_metrics(m).memory.iterative;
            end
            
            % Save the Results Images Struct
            save_file_name = [spath '\' name_str ...
                '_bin_size_' num2str(time_bin_size(k)) '_rep_' ...
                num2str(i) '.mat'];
            
            save(save_file_name, 'results');
            
            % Clean Up Memory
            clear temp_metrics results
            
            % Inform User This Round is Complete
            fprintf(['\nBenchmarked ', ...
                strrep(name_str, '_', ' '), ' at ', ...
                num2str(time_bin_size(k)), ' Wide Time Bins\n']);
        end
    end
    fprintf(['\nBenchmark Repetition ' num2str(i) ' Complete\n']);
end

% Estimate Image Accuracy
fprintf('\nImage Accuracy Estimation\n');
accuracy_percentage_matrix = RTFLIM_Accuracy_Estimation( ...
    num_reps, time_bin_size, num_methods, benchmark_files, ...
    benchmark_files_processed_FLIM, spath, data_order);

% Save Resultant Metrics 
save(['statistical_benchmarking_metrics_raw_' date '.mat'], '-v7.3');



%% Visualize Results
% if visualizer_flag == 1
%     fprintf(['\nVisualizing Results Across Benchmarks, Time Bins, ', ...
%         'and Repeats.\n']);
%     load('statistical_benchmarking_metrics_raw.mat');
%     RTFLIM_statistical_benchmarks_visualizer(statistical_metrics);
% end


%% Inform User of Completion
fprintf('\n\nAll Benchmarking is Complete.\n\n');