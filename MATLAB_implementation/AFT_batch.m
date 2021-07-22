%% load input files [.tif] %%

warning off

% load folder where input files are stored
uiwait(msgbox('Load parent folder'));
parent_d = uigetdir('');

matlab_folder = cd;
cd(parent_d)
listing = dir('*.tif');

% create output folder
mkdir('output')

cd(matlab_folder)

%% user input and set parameters %%

[parameters, listing_masks] = user_input(parent_d);

%% open one file at a time and perform analysis %%

n_files = length(listing);
av_ordermat = zeros(n_files,1);

for file_list = 1:n_files
    
    % file and directory name
    file = listing(file_list).name;
    directory = listing(file_list).folder;
    
    % file and directory name for mask (if local)
    if parameters.mask_method == 1
        file_mask = listing_masks(file_list).name;
        directory_mask = listing_masks(file_list).folder;
        parameters.mask_name = fullfile(directory_mask,file_mask);
    end
   
    % call function
    av_ordermat(file_list,1) = AFT_function(file, directory, parameters);
    
end

% save order parameter
save(fullfile([parent_d '/output'], 'median_order_parameter.mat'), 'av_ordermat');

T = table(av_ordermat);
T.Properties.VariableNames = {'median_order_parameter'};
writetable(T,fullfile([parent_d '/output'], 'median_order_parameter.csv'))

clear; clc