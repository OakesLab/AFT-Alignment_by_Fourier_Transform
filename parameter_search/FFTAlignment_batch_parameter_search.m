%% load input files [.tif] %%

warning off

% load folder where input files are stored [1st sample]
uiwait(msgbox('Load parent folder - 1st sample'));
parent_d1 = uigetdir('');

matlab_folder = cd;
cd(parent_d1)
listing1 = dir('*.tif');

% create output folder
mkdir('output_parameter_search')

cd(matlab_folder)

% load folder where input files are stored [2nd sample]
uiwait(msgbox('Load parent folder - 2nd sample'));
parent_d2 = uigetdir('');

cd(parent_d2)
listing2 = dir('*.tif');

% create output folder
mkdir('output_parameter_search')

cd(matlab_folder)

%% set parameters %%

% user input dialog
prompt = {'Minimum window size [px]', ...
    'Window size interval [px]'};
prompt_title = 'Parameters';
dims = [1 50];
definput = {'50','50'};
user_answer = inputdlg(prompt,prompt_title,dims,definput);

% save parameters in [output_parameter_search] folder (for both samples)
parameters_save.min_winsize_px = str2double(user_answer{1,1});
parameters_save.winsize_int_px = str2double(user_answer{2,1});
save(fullfile([parent_d1 '/output_parameter_search'], 'parameters.mat'), 'parameters_save');
save(fullfile([parent_d2 '/output_parameter_search'], 'parameters.mat'), 'parameters_save');

%% perform parameter search on both samples %%
av_ordermat_output1 = FFTAlignment_parameter_search_main(listing1, parameters_save);
av_ordermat_output2 = FFTAlignment_parameter_search_main(listing2, parameters_save);

%% plot results %%
% unpack cells (median values across all images in a sample)

n_windows = length(av_ordermat_output1);
n_neighbourhoods = size(av_ordermat_output1{1,1},2);
av_ordermat1 = NaN(n_windows, n_neighbourhoods);
av_ordermat2 = NaN(n_windows, n_neighbourhoods);

for k = 1:length(av_ordermat_output1)
    
    % for each cell (window size)
    temp1 = av_ordermat_output1{k,1};
    temp2 = av_ordermat_output2{k,1};
    
    % calculate median value across images (obtain one value for each
    % neighbourhood size)
    median_temp1 = median(temp1);
    median_temp2 = median(temp2);

    median_temp1(end+1:n_neighbourhoods) = NaN;
    median_temp2(end+1:n_neighbourhoods) = NaN;
    
    % save output -> row: window sizes, cols: neighbourhood sizes
    av_ordermat1(k,:) = median_temp1;
    av_ordermat2(k,:) = median_temp2;
    
end

% heatmap difference (sample1 - sample2)
figure
imagesc(av_ordermat1-av_ordermat2)
title('Order parameter difference (1st sample - 2nd sample)')
colormap('jet');
caxis([0 1])

xlabel('Neighbourhood [n of vectors]')
ylabel('Window size [px]')

%% check from here
yticks(1:n_windows);
yticklabels([10,20,29,39,49,59,69,79,89,99,109,119,129,139]);

xlabel('Neighbourhood [n of vectors]')
ylabel('Window size [µm]')

% xticks([5 10 15 20 25])
% xticklabels({'15x','30x','45x','60x','75x'})
xticklabels({'15x','30x','45x','60x','75x','90x','105x','120x'})

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 14)
