function RTFLIM_collected_benchmarks_visualizer( collected_metrics )
%% Real Time FLIM Benchmark Visualizer
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




%% Gather Metrics and Names
method_names = {collected_metrics(1).name(1).metrics.method};
num_benchmarks = numel(collected_metrics);
num_time_bins = numel(collected_metrics(1).name);
num_methods = numel(method_names);



%% Initialize Formatting Vectors
line_styles = {'-', '--', ':', '-.'}; 
dot_styles = {'o','+','*','d','^', 'p'};
color_styles = {'k', 'b', 'r', 'g', 'y', 'm'};
legend_loc = 'bestoutside'; % 'northwest' Might be a good option too.

% Display Names
disp_names = cell(1,(num_benchmarks*num_time_bins*num_methods));
categorized_names = cell(num_benchmarks,1);
for i = 1:num_benchmarks
    categorized_names{i} = strrep(collected_metrics(i).name_str, '_', ' ');
    for j = 1:num_time_bins
        for k = 1:num_methods
            ind = k + ((j-1)*num_methods) + ...
                ((i-1)*num_methods*num_time_bins);
            disp_names{ind} = [...
                collected_metrics(i).name(j).metrics(k).method ' of ' ...
                strrep(collected_metrics(i).name_str, '_', ' ') ...
                ' with ' ...
                num2str(collected_metrics(i).name(j).time_bin_size) ...
                ' Wide Time Bins'];
        end
    end
end

% Categorization is for bar graph later
categorized_names_temp = categorical(categorized_names);
categorized_names = reordercats(categorized_names_temp, categorized_names);



%% Visualize Time Data
% Time Metric Global Labels
time_y_label = 'Time [sec]';

% Iterative Time Metric Labels
time_iter_title = 'Algorithm Run Time During Data Acquisition';
time_iter_x_label = 'Iteration';

% Construct Iterative Time Figure
time_iter_fig = figure(1);
clf;
hold on;
for i = 1:num_benchmarks
    for j = 1:num_time_bins
        for k = 1:num_methods
            plot( 1:1:numel(...
                collected_metrics(i).name(j).metrics(k).time.iterative.time),...
                collected_metrics(i).name(j).metrics(k).time.iterative.time,...
                [line_styles{i}, dot_styles{j}, color_styles{k}]);
        end
    end
end
title(time_iter_title);
xlabel(time_iter_x_label);
ylabel(time_y_label);
legend(disp_names, 'Location', legend_loc);
hold off;


% Combined Time Metric Labels
time_comb_title = 'Algorithm Run Time in Fully Combined Data Set';
time_comb_x_label = 'Benchmark Set';

% Construct Combined Time Figure
time_comb_fig = figure(2);
clf;
hold on;
% Fill in combined time values
combined_time_values = zeros(num_benchmarks, (num_time_bins*num_methods));
leg_names_method_time_bin = cell(1, (num_time_bins*num_methods));
for i = 1:num_benchmarks
    for j = 1:num_time_bins
        for k = 1:num_methods
            ind = k + ((j-1) * num_methods);
            combined_time_values(i,ind) = ...
                collected_metrics(i).name(j).metrics(k).time.combined;
            if i == 1
                leg_names_method_time_bin{ind} = [...
                    collected_metrics(i).name(j).metrics(k).method ...
                    ' with ' ...
                    num2str(collected_metrics(i).name(j).time_bin_size) ...
                    ' Wide Time Bins'];
            end
        end
    end
end
bar(categorized_names, combined_time_values);
title(time_comb_title);
xlabel(time_comb_x_label);
ylabel(time_y_label);
legend(leg_names_method_time_bin, 'Location', legend_loc);
hold off;



%% Visualize Memory Data
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
for i = 1:num_benchmarks
    for j = 1:num_time_bins
        for k = 1:num_methods
            plot( 1:1:numel(...
                collected_metrics(i).name(j).metrics(k).memory.iterative.memory),...
                (collected_metrics(i).name(j).metrics(k).memory.iterative.memory / (2^30)),...
                [line_styles{i}, dot_styles{j}, color_styles{k}]);
        end
    end
end
title(memory_iter_title);
xlabel(memory_iter_x_label);
ylabel(memory_y_label);
legend(disp_names, 'Location', legend_loc);
hold off;



% Combined Memory Metric Labels
memory_comb_title = ...
    'Estimated Memory Usage of Algorithms in Fully Combined Data Set';
memory_comb_x_label = 'Benchmark Set';

% Construct Combined Memory Figure
memory_comb_fig = figure(4);
clf;
hold on;
% Fill in combined time values
combined_memory_values = zeros(num_benchmarks, ...
    (num_time_bins*num_methods));
for i = 1:num_benchmarks
    for j = 1:num_time_bins
        for k = 1:num_methods
            ind = k + ((j-1) * num_methods);
            combined_memory_values(i,ind) = ...
                collected_metrics(i).name(j).metrics(k).time.combined;
        end
    end
end
bar(categorized_names, combined_memory_values);
title(memory_comb_title);
xlabel(memory_comb_x_label);
ylabel(memory_y_label);
legend(leg_names_method_time_bin,'Location', legend_loc);
hold off;



%% Visualize Results
% Result Metric Global Labels
img_iter_start_title = 'FLIM Estimation after the First Iteration';
img_iter_mid_title = 'FLIM Estimation Halfway Through Iterations';
img_iter_last_title = 'FLIM Estimation after the Last Iteration';
img_comb_title = 'FLIM Estimation from a Fully Combined Data Set';

% Iteratively Construct Qualitative Results Figures for Each Method
results_figure_handles = struct;
for i = 1:num_benchmarks
    for j = 1:num_time_bins
        ind = j + ((i-1) * num_time_bins);
        
        % First Iteration
        results_figure_handles(ind).handle = figure(4 + ((ind*2)-1));
        clf;
        hold on;
        temp_name = cellstr(categorized_names);
        sgtitle([temp_name{i}, ' with ' ...
            num2str(collected_metrics(i).name(j).time_bin_size), ...
            ' Wide Time Bins On First Iteration']);
        for k = 1:num_methods
            subplot(ceil(num_methods/2),ceil(num_methods/2),k);
            imshow(collected_metrics(i).name(j).metrics.results.iterative(1).result, [])
            title(collected_metrics(i).name(j).metrics.method);
        end
        hold off;
        
        
        % Final Iteration
        results_figure_handles(ind).handle = figure(4 + (ind*2));
        clf;
        hold on;
        temp_name = cellstr(categorized_names);
        sgtitle([temp_name{i}, ' with ' ...
            num2str(collected_metrics(i).name(j).time_bin_size), ...
            ' Wide Time Bins On Final Iteration']);
        for k = 1:num_methods
            subplot(ceil(num_methods/2),ceil(num_methods/2),k);
            imshow(collected_metrics(i).name(j).metrics.results.iterative(end).result, [])
            title(collected_metrics(i).name(j).metrics.method);
        end
        hold off;
    end
end

end