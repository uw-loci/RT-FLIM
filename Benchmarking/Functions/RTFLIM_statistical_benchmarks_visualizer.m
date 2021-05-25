function RTFLIM_statistical_benchmarks_visualizer( statistical_metrics )
%% Real Time FLIM Benchmark Visualizer - High Level Statistics
%   By: Niklas Gahm
%   2020/11/20
%
%   This code takes the benchmark metrics and visualizes them in a readily
%   human interpretable manner. It is designed to be able to handle an
%   arbitrary number of benchmarked algorithms, as long as they all conform
%   to the benchmark output standards. This version is specifically for the
%   Real Time FLIM project.
%
%   2020/11/20 - Started
%   2020/11/22 - Finished
%   2021/01/25 - Modified for Statistical Metrics
% 
%   To Do:
%       - Implement accuracy comparison to SPCImage Ground Truth 




%% Gather Metrics and Names
num_runs = numel(statistical_metrics);
num_benchmarks = numel(statistical_metrics(1).metrics);
num_time_bins = numel(statistical_metrics(1).metrics(1).name);
method_names = ...
    {statistical_metrics(1).metrics(1).name(1).metrics.method};
num_methods = numel(method_names);


%% Initialize Formatting Vectors
line_styles = {'-', '--', ':', '-.'}; 
dot_styles = {'o','+','*','d','^', 'p'};
color_styles = {'k', 'b', 'r', 'g', 'y', 'm'};
legend_loc = 'bestoutside'; % 'northwest' Might be a good option too.

% Display Names
% disp_names = cell(1,(num_benchmarks*num_time_bins*num_methods));
% categorized_names = cell(num_benchmarks,1);
% for i = 1:num_benchmarks
%     categorized_names{i} = strrep(collected_metrics(i).name_str, '_', ' ');
%     for j = 1:num_time_bins
%         for k = 1:num_methods
%             ind = k + ((j-1)*num_methods) + ...
%                 ((i-1)*num_methods*num_time_bins);
%             disp_names{ind} = [...
%                 collected_metrics(i).name(j).metrics(k).method ' of ' ...
%                 strrep(collected_metrics(i).name_str, '_', ' ') ...
%                 ' with ' ...
%                 num2str(collected_metrics(i).name(j).time_bin_size) ...
%                 ' Wide Time Bins'];
%         end
%     end
% end
% 
% % Categorization is for bar graph later
% categorized_names_temp = categorical(categorized_names);
% categorized_names = reordercats(categorized_names_temp, categorized_names);



%% Visualize Time Data Matrics
%   This is based on making three levels of data pools then visualizing the
%   average and standard deviations thereof in a human interpretable
%   manner. The first data pool is all iterations across repeats for
%   specific benchmark + number of time gates.

% Generate Data Pools
%   Generate a 5D data set that can be queried in different dimensinos to
%   acquire specific statistics specifically in the dimension order:
%   iteration point, time_gates, method, repeat, benchmark file (This is on
%   the cell level dus to the variability in benchmark file length.)
raw_time_data = cell(1, num_benchmarks);
for i = 1:num_benchmarks
    num_points = numel(...
        statistical_metrics(1).metrics(i).name(1).metrics(1).time.iterative.time);
    benchmark_raw_data = zeros(num_points, num_time_bins, ...
        num_methods, num_runs);
    for j = 1:num_runs
        for k = 1:num_methods
            for m = 1:num_time_bins
                benchmark_raw_data(:,m,k,j) = ...
                    statistical_metrics(j).metrics(i).name(m).metrics(k).time.iterative.time;
            end
        end
    end
    raw_time_data{i} = benchmark_raw_data;
end



   



%% Visualize Result Accuracy Metrics
%   This requires a comparison to the "Gold Standard" ground truth of a
%   processed image. Ideally the processed image from SPCImage 


end