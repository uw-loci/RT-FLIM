function [ photon_data ] = img_loader_RTFLIM_Bench(...
    benchmark_file_path, benchmark_file_name, home_path, data_order)
%% Image Loader
%   By: Niklas Gahm
%   2020/11/13
%
%   This code programatically loads in hdf5 files and generates the 
%   apropriate structs needed for the remainder of the benchmarking code
%
%   2020/11/13 - Started
%   2020/11/15 - Finished




%% Navigate to Folder Containing the Data
cd(benchmark_file_path);



%% Setup Photon Data Structs and Fill in photon_data.counts
photon_data = struct;

% Use Matlab's native HDF5 reader to import data.
file_info = h5info(benchmark_file_name);
num_datasets = numel(file_info.Datasets);


loader_bar = waitbar((1/num_datasets), 'Loading Photon Data.');

for i = 1:num_datasets
    waitbar((i/num_datasets), loader_bar);
    
    % Read counts data and convert to double
    photon_data(i).counts = double(h5read(benchmark_file_name, ...
        ['/', file_info.Datasets(i).Name]));
    
    % Rearrange dimensions to fit the benchmarking code
    switch data_order
        case 'XYTS'
            % This is what the code is setup for, don't need to do anything
        case 'TXYS'
            % This is the order our FLIM machines commonly use
            photon_data(i).counts = permute(photon_data(i).counts, ...
                [2,3,1]);
        % Added additional orders, just in case.
        case 'YXTS'
            photon_data(i).counts = permute(photon_data(i).counts, ...
                [2,1,3]);
        case 'TYXS'
            photon_data(i).counts = permute(photon_data(i).counts, ...
                [3,2,1]);
        case 'XTYS'
            photon_data(i).counts = permute(photon_data(i).counts, ...
                [1,3,2]);
        case 'YTXS'
            photon_data(i).counts = permute(photon_data(i).counts, ...
                [3,1,2]);
        otherwise
            error('Unsupported order of dimensions.');
    end
end
close(loader_bar);



%% Clean Up Navigation
cd(home_path);

end

