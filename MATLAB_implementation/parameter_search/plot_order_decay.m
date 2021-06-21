
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

cd(matlab_folder)

%% set parameters %%

% user input dialog
prompt = {'Window size to display [px]'};
prompt_title = 'Display decay';
dims = [1 50];
user_answer = inputdlg(prompt,prompt_title,dims);

%% fetch data %%

av_ordermat1 = load(fullfile([parent_d1, '/output_parameter_search'], 'median_order_parameter_search.mat'));
av_ordermat2 = load(fullfile([parent_d2, '/output_parameter_search'], 'median_order_parameter_search.mat'));
parameters = load(fullfile([parent_d2, '/output_parameter_search'],'parameters.mat'));

av_ordermat1 = av_ordermat1.av_ordermat_output;
av_ordermat2 = av_ordermat2.av_ordermat_output;
parameters = parameters.parameters_save;

%% extract order parameter for requested window size %%

winsize = str2double(user_answer{1,1});
array = parameters.min_winsize_px:...
    parameters.winsize_int_px:...
    (length(av_ordermat1)-1)*parameters.winsize_int_px+parameters.min_winsize_px;

cond = winsize == array;
out1 = av_ordermat1{cond,1};
out2 = av_ordermat2{cond,1};

decay1 = nanmean(out1);
decay2 = nanmean(out2);

%% plot and save %%
figure

plot(decay1, 'ko', 'MarkerSize', 10);
hold on
plot(decay2, 'mo', 'MarkerSize', 10);

title_str = ['Order parameter decay for ' num2str(winsize) ' px window size'];
title(title_str)

% labels
xlabel('Neighbourhood (vectors)')
ylabel('Order parameter (-)')
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 14)

n_neighbourhoods = length(decay1);
x_tick = 1:n_neighbourhoods;
neighbourhood_list = x_tick*2+1; % neighbourhood = 2*st+1
x_label = string(neighbourhood_list);
x_label = x_label + 'x';

xticks(x_tick);
xticklabels(x_label);
xtickangle(45)

legend('Sample 1', 'Sample 2', 'box', 'off');

% save figure
im_out = getframe(gcf);
im_out = im_out.cdata;
imwrite(im_out, fullfile([parent_d1 '/output_parameter_search'], ['decay_' num2str(winsize) 'px_window.tif']));
imwrite(im_out, fullfile([parent_d2 '/output_parameter_search'], ['decay_' num2str(winsize) 'px_window.tif']));
close

% save order parameter
save(fullfile([parent_d1 '/output_parameter_search'], ['median_order_parameter_' num2str(winsize) 'px.mat']), 'decay1');
save(fullfile([parent_d2 '/output_parameter_search'], ['median_order_parameter_' num2str(winsize) 'px.mat']), 'decay2');

T1 = table(decay1');
T1.Properties.VariableNames = {['median_order_parameter_sample1_' num2str(winsize) 'px']};
writetable(T1,fullfile([parent_d1 '/output_parameter_search'], ['median_order_parameter_' num2str(winsize) 'px.csv']))
T2 = table(decay2');
T2.Properties.VariableNames = {['median_order_parameter_sample2_' num2str(winsize) 'px']};
writetable(T2,fullfile([parent_d2 '/output_parameter_search'], ['median_order_parameter_' num2str(winsize) 'px.csv']))

clear; clc