function [ metrics ] = RTFLIM_Benchmarking_Framework( ...
    benchmark_file, data_order, time_bin_size, visualizer_flag, lite_flag)
%% Runtime FLIM Benchmarking Framework
%   By: Niklas Gahm
%   2020/11/12
%
%   This is a framework that loads up the needed data to benchmark
%   approaches for runtime FLIM. 
% 
%   2020/11/12 - Started 
% 
%   To-Do:
%       - 




%% Initialize Variables
metrics = struct;


%% Navigation Setup
fprintf('\n\nGenerating Paths\n');
home_path = pwd;
[benchmark_file_path, file_name, file_ext] = fileparts(benchmark_file);
benchmark_file_name = [file_name, file_ext];



%% Load in Benchmark Data
fprintf('\nLoading Photon Data\n');
photon_data = img_loader_RTFLIM_Bench(benchmark_file_path, ...
    benchmark_file_name, home_path, data_order);



%% Perform Time-Binning
fprintf('\nTime-Binning Data\n');
[photon_data, time_bin_size] = ...
    photon_time_binning_RTFLIM_Bench(photon_data, time_bin_size);



%% Generate Combined Data Set 
combined_data = photon_data(1).counts;
for i = 2:numel(photon_data)
    combined_data = combined_data + photon_data(i).counts;
end



%% Test NC-PCA
% Based on https://onlinelibrary.wiley.com/doi/full/10.1002/jbio.201600160
fprintf('\nBenchmarking PCA\n');
[ PCA_time, PCA_memory, PCA_results ] = ...
    benchmarker_PCA( photon_data, combined_data, lite_flag );

% Assign outuputs to a useable struct;
metrics(1).method = 'NC-PCA';
metrics(1).time = PCA_time;
metrics(1).memory = PCA_memory;
metrics(1).results = PCA_results;



%% Test Phasor
% Based on https://github.com/PirminLakner/Phasor_FLIM
fprintf('\nBenchmarking Phasor\n');
[ Phasor_time, Phasor_memory, Phasor_results ] = ...
    benchmarker_Phasor( photon_data, combined_data, lite_flag, ...
    time_bin_size );

% Assign outuputs to a useable struct;
metrics(2).method = 'Phasor';
metrics(2).time = Phasor_time;
metrics(2).memory = Phasor_memory;
metrics(2).results = Phasor_results;



% %% Test Laguerre 
% % Based on https://ieeexplore.ieee.org/document/5594996
% fprintf('\nBenchmarking Laguerre\n');
% [ Laguerre_time, Laguerre_memory, Laguerre_results ] = ...
%     benchmarker_Laguerre( photon_data, combined_data, lite_flag );
% 
% % Assign outuputs to a useable struct;
% metrics(3).method = 'Laguerre';
% metrics(3).time = Laguerre_time;
% metrics(3).memory = Laguerre_memory;
% metrics(3).results = Laguerre_results;



%% Visualize Benchmarking Results
if visualizer_flag == 1
    fprintf('\nVisualizing Results\n');
    RTFLIM_benchmark_visualizer( metrics );
end



%% Return to Starting Point
cd(home_path);



%% Confirm Completion
fprintf('\nBenchmarking Complete\n\n\n');

end