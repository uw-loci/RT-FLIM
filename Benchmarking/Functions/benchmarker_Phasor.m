function [ Phasor_time, Phasor_memory, Phasor_results ] = ...
    benchmarker_Phasor( ...
    photon_data, combined_data, lite_flag, time_bin_size )
%% Phasor Benchmarker
%   By: Niklas Gahm
%   2020/11/16
%
%   This code takes the read in photon data and uses it for Phasor FLIM
%   Estimation. This process is benchmarked for time and memory usage. It
%   further outputs Phasor's iterative results for qualitative assessment.
% 
%   This code is based off of the code available at:
%   https://github.com/PirminLakner/Phasor_FLIM
%
%   2020/11/16 - Started



%% Initialize Outputs
% Initialize Time Struct
Phasor_time = struct;
Phasor_time.iterative(1) = 0;
Phasor_time.combined = 0;

% Initialize Memory Struct
Phasor_memory = struct;
Phasor_memory.iterative(1) = 0;
Phasor_memory.combined = 0;

% Initialize Results Struct
Phasor_results = struct;
Phasor_results.iterative = struct;
Phasor_results.iterative.result(1) = 0;
Phasor_results.combined = 0;
Phasor_results.mid_iter_ind = round(numel(photon_data)/2);



%% Inform user to be Patient
wait_box = msgbox({'The script is now going into measurement mode.'; ... 
    'No waitbar will be shown to get accurate timing data.'}, ...
    'Phasor Benchmarking');



%% Setting up the Constant Variables Throughout Measurement

% -- set constants
thresh = 0.3;               % threshold [0,1]
harmonic = 2;               % harmonic number (1,2,3,4,...)

shift_bins = 0;             % number of bins shifted to right from max for 
                            % minimizing effect of IRF

freq0 = 7.6E+7;                   % laser frequency, here 80 MHz
delta_t = 48E-12;           % width of one time channel
delta_t = delta_t * time_bin_size;
time_bins = size(photon_data(1).counts,3);
% number of time bins for one period (not 64!). 80 for harmonic no 1

% ------ initial calculations ------ %

freq = harmonic * freq0;        % virtual frequency for higher harmonics
w = 2 * pi * freq ;             % angular frequency



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
Phasor_memory.combined = 6*(temp_list(temp_ind).bytes);



% Start Timing
combined_phasor_timer = tic;

% Get Combined Data and reshape it as needed
E = reshape(combined_data, numel(combined_data(:,:,1)), time_bins);


% Find Max and remove data before max
maxE = sum(E,1);
[~,I] = max(maxE);
decdata = E(:,I+shift_bins:end);
timechannels_data = size(decdata,2);

if timechannels_data < time_bins
    decdata(1,time_bins) = 0;
elseif  timechannels_data > time_bins
    decdata = decdata(:,1:time_bins);
end


% remove offset from data
data_off = mean(E(:,round(I/3):round(2*I/3)),2);
data_off = repmat(data_off,1,size(decdata,2));
decdata = decdata - data_off ;
decdata(decdata<0) = 0 ;


% Threshold (not peak vs peak. Used whole intensity vs max whole intensity)
maxmax = max(sum(decdata,2));
rows_to_remove = any(sum(decdata,2) < (thresh * maxmax), 2);
decdata(rows_to_remove,:) = [];

% calculate G/S-sin/cos-matrix
tb_vec = linspace(1,time_bins,time_bins)'; % time bin indexing vector

Gn_ma = cos(w * delta_t * (tb_vec - 0.5));
Sn_ma = sin(w * delta_t * (tb_vec - 0.5));

% calculate data phasor
Gn = double(decdata) * double(Gn_ma) ;  %take decdata from preprocessing
Sn = double(decdata) * double(Sn_ma) ;
area = sum(decdata, 2) ;

% normalization
G_f = Gn ./ area;
S_f = Sn ./ area;

% data
Z = [G_f,S_f];


% Abbriviated Phasor Since a Score Value can technically get generated at
% this point. This Score value can in turn be used to present an estimated
% view finder

score = zeros(numel(rows_to_remove),1);
counter = 0;
for i = 1:numel(score)
    if ~rows_to_remove(i)
        % This case the pixel here has some life time and needs a score
        % assigned to it
        counter = counter + 1;
        score(i) = Z(counter,1)*Z(counter,2);
    end
    % Otherwise the pixel doesn't have a lifetime.
end
score = score / max(score, [], 'all');
Phasor_results.combined = ...
    reshape(score, size(combined_data,1), size(combined_data,2));

% End Timing
Phasor_time.combined = toc(combined_phasor_timer);



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
    
    % Memory estimate is centered on the largest matrix and the operations
    % done to is. This is a rough estimate, not an exact number.
    Phasor_memory.iterative(i) = 7*(temp_list(temp_ind).bytes);
    
    
    % Start Timing
    iterative_phasor_timer = tic;
    
    % Build the Cumulative Counts
    cumulative_counts = cumulative_counts + photon_data(i).counts;
    
    % Get Cumulative Counts and reshape it as needed
    E = reshape(cumulative_counts, numel(cumulative_counts(:,:,1)), time_bins);
    
    
    % Find Max and remove data before max
    maxE = sum(E,1);
    [~,I] = max(maxE);
    decdata = E(:,I+shift_bins:end);
    timechannels_data = size(decdata,2);
    
    if timechannels_data < time_bins
        decdata(1,time_bins) = 0;
    elseif  timechannels_data > time_bins
        decdata = decdata(:,1:time_bins);
    end
    
    
    % remove offset from data
    data_off = mean(E(:,round(I/3):round(2*I/3)),2);
    data_off = repmat(data_off,1,size(decdata,2));
    decdata = decdata - data_off ;
    decdata(decdata<0) = 0 ;
    
    
    % Threshold (not peak vs peak. Used whole intensity vs max whole intensity)
    maxmax = max(sum(decdata,2));
    rows_to_remove = any(sum(decdata,2) < (thresh * maxmax), 2);
    decdata(rows_to_remove,:) = [];
    
    % calculate G/S-sin/cos-matrix
    tb_vec = linspace(1,time_bins,time_bins)'; % time bin indexing vector
    
    Gn_ma = cos(w * delta_t * (tb_vec - 0.5));
    Sn_ma = sin(w * delta_t * (tb_vec - 0.5));
    
    % calculate data phasor
    Gn = double(decdata) * double(Gn_ma) ;  %take decdata from preprocessing
    Sn = double(decdata) * double(Sn_ma) ;
    area = sum(decdata, 2) ;
    
    % normalization
    G_f = Gn ./ area;
    S_f = Sn ./ area;
    
    % data
    Z = [G_f,S_f];
    
    
    % Abbriviated Phasor Since a Score Value can technically get generated 
    % at this point. This Score value can in turn be used to present an 
    % estimated view finder
    
    score = zeros(numel(rows_to_remove),1);
    counter = 0;
    for j = 1:numel(score)
        if ~rows_to_remove(j)
            % This case the pixel here has some life time and needs a score
            % assigned to it
            counter = counter + 1;
            score(j) = Z(counter,1)*Z(counter,2);
        end
        % Otherwise the pixel doesn't have a lifetime.
    end
    score = score / max(score, [], 'all');
    
    if lite_flag == 1
        if i == 1
            Phasor_results.iterative(1).result = ...
                reshape(score, size(photon_data(i).counts, 1), ...
                size(photon_data(i).counts, 2));
        elseif i == Phasor_results.mid_iter_ind
            Phasor_results.iterative(2).result = ...
                reshape(score, size(photon_data(i).counts, 1), ...
                size(photon_data(i).counts, 2));
        elseif i == numel(photon_data)
            Phasor_results.iterative(3).result = ...
                reshape(score, size(photon_data(i).counts, 1), ...
                size(photon_data(i).counts, 2));
        end
    else
        Phasor_results.iterative(i).result = ...
            reshape(score, size(photon_data(i).counts, 1), ...
            size(photon_data(i).counts, 2));
    end
    
    % End Timing
    Phasor_time.iterative(i) = toc(iterative_phasor_timer);
end



%% Cleanup from Benchmarking
if ishandle(wait_box)
    close(wait_box);
end


end