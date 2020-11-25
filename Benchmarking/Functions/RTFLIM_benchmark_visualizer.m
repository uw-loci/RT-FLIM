function RTFLIM_benchmark_visualizer( metrics )
%% Real Time FLIM Benchmark Visualizer
%   By: Niklas Gahm
%   2020/11/16
%
%   This code takes the benchmark metrics and visualizes them in a readily
%   human interpretable manner. It is designed to be able to handle an
%   arbitrary number of benchmarked algorithms, as long as they all conform
%   to the benchmark output standards. This version is specifically for the
%   Real Time FLIM project.
%
%   2020/11/16 - Started
%   2020/11/17 - Finished




%% Gather Metrics and Names
method_names = {metrics.method};
num_methods = numel(metrics);



%% Initialize Formatting Vectors
line_styles = {'-ok', '-+b', '-dg', '-|g'}; 
legend_loc = 'best'; % 'northwest' Might be a good option too.

% Categorization is for bar graph later
categorized_methods = categorical(method_names);
categorized_methods = reordercats(categorized_methods, method_names);



%% Visualize Time Metrics
% Time Metric Global Labels
time_y_label = 'Time [sec]';

% Iterative Time Metric Labels
time_iter_title = 'Algorithm Run Time During Data Acquisition';
time_iter_x_label = 'Iteration';

% Construct Iterative Time Figure
time_iter_fig = figure(1);
clf;
hold on;
for i = 1:num_methods
    plot( 1:1:numel(metrics(i).time.iterative.time), ...
        metrics(i).time.iterative.time, line_styles{i});
end
title(time_iter_title);
xlabel(time_iter_x_label);
ylabel(time_y_label);
legend(method_names, 'Location', legend_loc);
hold off;


% Combined Time Metric Labels
time_comb_title = 'Algorithm Run Time in Fully Combined Data Set';
time_comb_x_label = 'Method';

% Construct Combined Time Figure
time_comb_fig = figure(2);
clf;
hold on;
% Fill in combined time values
combined_time_values = zeros(1, num_methods);
for i = 1:num_methods
    combined_time_values(i) = metrics(i).time.combined;
end
bar(categorized_methods, combined_time_values);
title(time_comb_title);
xlabel(time_comb_x_label);
ylabel(time_y_label);
hold off;



%% Visualize Memory Metrics
% Memory Metric Global Labels
memory_y_label = 'Estimated Memory Usage [GB]';

% Iterative Memory Metric Labels
memory_iter_title = ...
    'Estimated Algorithmic Memory Usage During Data Acquisition';
memory_iter_x_label = 'Iteration';

% Construct Iterative Memory Figure
memory_iter_fig = figure(3);
clf;
hold on;
for i = 1:num_methods
    plot( 1:1:numel(metrics(i).memory.iterative.memory), ...
        (metrics(i).memory.iterative.memory / (2^30)), ...
        line_styles{i});
end
title(memory_iter_title);
xlabel(memory_iter_x_label);
ylabel(memory_y_label);
legend(method_names, 'Location', legend_loc);
hold off;

% Combined Memory Metric Labels
memory_comb_title = ...
    'Estimated Memory Usage of Algorithms in Fully Combined Data Set';
memory_comb_x_label = 'Method';

% Construct Combined Memory Figure
memory_comb_fig = figure(4);
clf;
hold on;
% Fill in combined memory values
combined_memory_values = zeros(1, num_methods);
for i = 1:num_methods
    combined_memory_values(i) = metrics(i).memory.combined  / (2^30);
end
bar(categorized_methods, combined_memory_values);
title(memory_comb_title);
xlabel(memory_comb_x_label);
ylabel(memory_y_label);
hold off;



%% Visualize Results
% Result Metric Global Labels
img_iter_start_title = 'FLIM Estimation after the First Iteration';
img_iter_mid_title = 'FLIM Estimation Halfway Through Iterations';
img_iter_last_title = 'FLIM Estimation after the Last Iteration';
img_comb_title = 'FLIM Estimation from a Fully Combined Data Set';

% Iteratively Construct Qualitative Results Figures for Each Method
results_figure_handles = struct;
for i = 1:num_methods
    
    results_figure_handles(i).handle = figure(4 + i);
    clf;
    hold on;
    sgtitle(method_names{i});
    
    % First Iteration
    subplot(2,2,1);
    imshow(metrics.results.iterative(1).result, [])
    title(img_iter_start_title);
    
    % Middle Iteration
    subplot(2,2,2);
    imshow(metrics.results.iterative(round(numel(...
        metrics.results.iterative)/2)).result, []);
    title(img_iter_mid_title);
    
    % Last Iteration
    subplot(2,2,3);
    imshow(metrics.results.iterative(end).result, [])
    title(img_iter_last_title);
    
    % Combined
    subplot(2,2,4);
    imshow(metrics.results.combined, []);
    title(img_comb_title);
    
    hold off;
end



end

