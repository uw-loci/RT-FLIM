function accuracy_percentage_matrix = RTFLIM_Accuracy_Estimation( ...
    num_reps, time_bin_size, num_methods, benchmark_files, ...
    benchmark_files_processed_FLIM, img_path, data_order)
%% RT-FLIM Accuracy Estimation
%   By: Niklas Gahm
%   2021/06/13
%
%   This code loads the resultant images from the benchmarking process and
%   generates an accuracy percentage for each image. Only the results from
%   the iterative process are analyzed since the last cumulative image is
%   equivalent to the combined data image. 
% 
%   2021/06/13 - Started




% Setup useful variables
home_path = pwd;
num_time_bins = numel(time_bin_size);
num_bench_files = numel(benchmark_files);
accuracy_percentage_matrix = cell(1, num_bench_files);
accuracy_waitbar = waitbar(0, 'Image Accuracy Estimation');
counter = 0;

% Iterate on benchmark files
for i_file = 1:num_bench_files
    first_load_flag = 1;
    %% Load Benchmark File
    [benchmark_file_path, name_str, file_ext]  = ...
        fileparts(benchmark_files{i_file});
    benchmark_file_name = [name_str, file_ext];
    fprintf('\nLoading Photon Data\n');
    photon_data = img_loader_RTFLIM_Bench(benchmark_file_path, ...
        benchmark_file_name, home_path, data_order);
    
    %% Generate Combined Data Set
    combined_data = photon_data(1).counts;
    for i = 2:numel(photon_data)
        combined_data = combined_data + photon_data(i).counts;
    end
    
    % Calculate Image size
    img_size = size(combined_data);
    
    %% Load Fitted Benchmark File
    gold_std_img = double(imread(benchmark_files_processed_FLIM{i_file}));
    
    %% Rescale Gold Standard Image to [0,1]
    gold_std_img = gold_std_img - min(gold_std_img, [], 'all');
    gold_std_img = gold_std_img ./ max(gold_std_img, [], 'all');
    
    % Calculate Maximum Possible Error
    max_error = img_size(1) * img_size(2);
    
    % Iterate on Time Bins
    for i_bin = 1:num_time_bins
        % Iterate on Repetitions
        for i_rep = 1:num_reps
            %% Load Img File
            % Get the file name to load
            load_name = [img_path '\' name_str '_bin_size_' ...
                num2str(time_bin_size(i_bin)) '_rep_' num2str(i_rep) ...
                '.mat'];
            
            % Load file
            temp = load(load_name);
            
            % Iterate on Methods
            for i_method = 1:num_methods
                % Get the number of iterative_images
                num_iter_img = numel(temp.results{1}.iterative);
                
                if first_load_flag == 1
                    first_load_flag = 0;
                    accuracy_percentage_matrix{i_file} = zeros(num_reps,...
                        num_iter_img, num_time_bins, num_methods);
                end
                
                % Iterate on Each Iterative Image
                for i_img = 1:num_iter_img
                    %% Calculate Accuracy Percentage for Each Iterative Img
                    % Update the waitbar
                    counter = counter + 1;
                    waitbar( (counter/(num_bench_files * num_time_bins ...
                        * num_methods * num_reps * num_iter_img)), ...
                        accuracy_waitbar);
                    
                    % Get Iterative Image
                    iter_img = ...
                        temp.results{i_method}.iterative(i_img).result;
                    
                    % Rescale Iterative Image to [0,1]
                    iter_img = iter_img - min(iter_img, [], 'all');
                    iter_img = iter_img ./ max(iter_img, [], 'all');
                    
                    % Error Value
                    error_val = sum(abs(gold_std_img - iter_img), 'all');
                    
                    % Scale to percentage 
                    scaled_error = error_val / max_error;
                    accuracy_percentage_matrix{i_file}( i_rep, i_img, ...
                        i_bin, i_method) = (1 - scaled_error) * 100;
                    
                end
            end
        end
    end
end

close(accuracy_waitbar);
end