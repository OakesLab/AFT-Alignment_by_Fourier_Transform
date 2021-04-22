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

%% set parameters %%

prompt = {'Window size [px]', ...
    'Window overlap [%]',...
    'Neighbourhood size [vectors]',...
    'Masking method (local = 1, global = 0)',...
    'Display figures (yes = 1, no = 0)'};
prompt_title = 'Parameters';
dims = [1 50];
definput = {'250','50','5','0', '1'};
user_answer = inputdlg(prompt,prompt_title,dims,definput);

parameters.winsize = str2double(user_answer{1,1});  
parameters.overlap = str2double(user_answer{2,1})/100;  
parameters.st = round((str2double(user_answer{3,1})- 1)/2); 
parameters.mask_method = str2double(user_answer{4,1}); % global = 1; local = 2;
parameters.figures = str2double(user_answer{5,1});

parameters.checkpoint = 0; % threshhold sum in each window to do the calculation

% load masks if masking method == local
if parameters.mask_method == 1
    
    uiwait(msgbox('Load folder containing masks'));
    masks_d = uigetdir('');

    matlab_folder = cd;
    cd(masks_d)
    listing_masks = dir('*.tif');

    cd(matlab_folder)
end


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
    av_ordermat(file_list,1) = FFTAlignment(file, directory, parameters);
    
end

save(fullfile([parent_d '/output'], 'average_ordermat.mat'), 'av_ordermat');
clear; clc