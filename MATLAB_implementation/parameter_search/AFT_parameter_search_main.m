function av_ordermat_output = AFT_parameter_search_main(listing, parameters)

%% calculate parameters for analysis
% image size
im_size_info = imfinfo(fullfile(listing(1).folder, listing(1).name));
im_size = im_size_info.Width; % [px]

% find window size array
min_win_size_px = parameters.min_winsize_px;        % [px]
max_win_size_px = floor(im_size/3); % [px] - this is arbitrary - considering an overlap of 50% this would give 3 windows (vectors) across
interval_win_size_px = parameters.winsize_int_px;   % [px]
win_size_array = min_win_size_px:interval_win_size_px:max_win_size_px;

% winsize has to be odd - make odd if even
cond =  mod(win_size_array,2) == 0; % number is even
win_size_array(cond) = win_size_array(cond) - 1;

% set overlap to value selected by user
overlap = 1 - parameters.overlap_percentage/100;

%% open one file at a time and perform analysis %%

n_files = length(listing);
av_ordermat_output = {};

% for every window size
for window_list = 1:length(win_size_array)
    
    win_size = win_size_array(window_list);	% [px] has to be odd
    
    win_rad = floor(win_size/2);
    win_space = ceil(win_size*overlap);
    
    % find anglemat for every file (vector field)
    av_ordermat = [];
    for file_list = 1:n_files
        
        % file and directory name
        file = listing(file_list).name;
        directory = listing(file_list).folder;
        
        % load images
        im = imread(fullfile(directory, file));
        im = im2double(im);
        
        % call function to calculate vector field
        anglemat = AFT_parameter_search_anglemat(im, win_size, win_rad, win_space);
            
        % for every neighbourhood size
        n_win = length(win_rad+1:win_space:im_size-win_rad);
        max_neighbourhood = floor((n_win-1)/2);
        for neighbourhood_list = 1:max_neighbourhood
            st = neighbourhood_list;	% 2*st+1 window size to calculate order parameter
            
            % output -> rows: each image, cols: neighbourhoods
            av_ordermat(file_list,st) = AFT_parameter_search_ordermat(anglemat,st);
        end
        
    end
    
    % output cell: one cell for each window size
    av_ordermat_output{window_list, 1} = av_ordermat; 
    
end 

%% save order parameter
save(fullfile([listing(1).folder '/output_parameter_search'], 'median_order_parameter_search.mat'), 'av_ordermat_output');

end