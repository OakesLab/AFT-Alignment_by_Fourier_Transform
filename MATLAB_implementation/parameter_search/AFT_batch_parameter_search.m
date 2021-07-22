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
    'Window size interval [px]', ...
    'Window overlap (fixed parameter) [%]'};
prompt_title = 'Parameters';
dims = [1 50];
definput = {'50','50','50'};
user_answer = inputdlg(prompt,prompt_title,dims,definput);

% save parameters in [output_parameter_search] folder (for both samples)
parameters_save.min_winsize_px = str2double(user_answer{1,1});
parameters_save.winsize_int_px = str2double(user_answer{2,1});
parameters_save.overlap_percentage = str2double(user_answer{3,1});
save(fullfile([parent_d1 '/output_parameter_search'], 'parameters.mat'), 'parameters_save');
save(fullfile([parent_d2 '/output_parameter_search'], 'parameters.mat'), 'parameters_save');

T = struct2table(parameters_save);
writetable(T,fullfile([parent_d1 '/output_parameter_search'], 'parameters.txt'))
writetable(T,fullfile([parent_d2 '/output_parameter_search'], 'parameters.txt'))

%% perform parameter search on both samples %%
av_ordermat_output1 = AFT_parameter_search_main(listing1, parameters_save);
av_ordermat_output2 = AFT_parameter_search_main(listing2, parameters_save);

%% calculate median values and p-values %%

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

% calculate p-value (Mann-Whitney test)
p_median = NaN(size(av_ordermat1));

for k = 1:length(av_ordermat_output1)
    
    % for each cell (window size)
    temp1 = av_ordermat_output1{k,1};
    temp2 = av_ordermat_output2{k,1};
    
    % Mann-Whitney of all possible combinations of window/neighbourhood size
    for kk = 1:size(temp1,2)
        
        x_median = temp1(:,kk);
        y_median = temp2(:,kk);
        p_median(k,kk) = ranksum(x_median,y_median);
        
    end  
    
end

%% plot results %%
% heatmap difference (sample1 - sample2)
figure
h = imagesc(av_ordermat1-av_ordermat2);
title('Order parameter difference (1st sample - 2nd sample)')
colormap('jet');
colorbar()

% NaNs as black
set(h, 'AlphaData', ~isnan(av_ordermat1-av_ordermat2))
axis on
set(gca, 'Color', 'k')

% labels
xlabel('Neighbourhood (vectors)')
ylabel('Window size (px)')
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 14)

x_tick = 1:n_neighbourhoods;
neighbourhood_list = x_tick*2+1; % neighbourhood = 2*st+1
x_label = string(neighbourhood_list);
x_label = x_label + 'x';

xticks(x_tick);
xticklabels(x_label);
xtickangle(45)
yticks(1:n_windows);
yticklabels(parameters_save.min_winsize_px:...
    parameters_save.winsize_int_px:...
    (n_windows-1)*parameters_save.winsize_int_px+parameters_save.min_winsize_px);

% save
im_out = getframe(gcf);
im_out = im_out.cdata;
imwrite(im_out, fullfile([parent_d1 '/output_parameter_search'], 'parameter_search_difference.tif'));
imwrite(im_out, fullfile([parent_d2 '/output_parameter_search'], 'parameter_search_difference.tif'));
close

% heatmap p-value 
figure
h1 = imagesc(p_median);
title('p-value (Mann-Whitney)');
c_spring = colormap('spring');
c_spring = flipud(c_spring);
colormap(c_spring)
colorbar();

% NaNs as black
set(h1, 'AlphaData', ~isnan(av_ordermat1-av_ordermat2))
axis on
set(gca, 'Color', 'k')

% labels
xlabel('Neighbourhood (vectors)')
ylabel('Window size (px)')
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 14)

x_tick = 1:n_neighbourhoods;
neighbourhood_list = x_tick*2+1; % neighbourhood = 2*st+1
x_label = string(neighbourhood_list);
x_label = x_label + 'x';

xticks(x_tick);
xticklabels(x_label);
xtickangle(45)
yticks(1:n_windows);
yticklabels(parameters_save.min_winsize_px:...
    parameters_save.winsize_int_px:...
    (n_windows-1)*parameters_save.winsize_int_px+parameters_save.min_winsize_px);

% save
im_out = getframe(gcf);
im_out = im_out.cdata;
imwrite(im_out, fullfile([parent_d1 '/output_parameter_search'], 'parameter_search_p-value.tif'));
imwrite(im_out, fullfile([parent_d2 '/output_parameter_search'], 'parameter_search_p-value.tif'));
close

clear; clc