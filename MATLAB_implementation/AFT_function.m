%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code and documentation can be found at
% https://github.com/OakesLab/AFT-Alignment_by_Fourier_Transform
%
% This routine uses a vector field of alignment directions using small
% sub-windows in the real space image to calculate an alignment order
% parameter.  Images should be grayscale images.
% Angles determined are oriented as follows:
%
%                            ^  180°
%                            |
%                            |
%                            ------>  90°
%                            |
%                            |
%                            v  0°
%
% Order Parameter value:   0 = Completely random alignment
%                          1 = Perfectly aligned
%
% All rights and permissions belong to:
% Patrick Oakes
% poakes@gmail.com
% 05/15/2015
%
% Citation:
% Cetera et al. Nat Commun 2014; 5:1-12
% http://www.ncbi.nlm.nih.gov/pubmed/25413675
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function av_ordermat = AFT_function(file, directory, parameters)

% load image
im = imread(fullfile(directory, file));
im = im2double(im);

% calculate angle vector field
[anglemat,pc,pr,vc2,ur2] = AFT_anglemat(im, parameters);

% plots 
if parameters.figures == 1
    
    % vector field
    figure
    imshow(im, [])
    hold on
    if exist('pc','var') ~= 0
        quiver(pc,pr,vc2,ur2,0,'y','showarrowhead','off','linewidth',2)
    end
    im_out = getframe(gcf);
    im_out = im_out.cdata;
    imwrite(im_out, fullfile([directory '/output'], ['vectors_' file(1:end-4) '.tif']));
    close
    
    % angle heat map
    figure('Position',[100 100 size(im,2)/3 size(im,1)/3]);
    imagesc(rad2deg(anglemat));
    hsv_nan = [[0,0,0];colormap('hsv')];
    set(gca,'visible','off')
    caxis([0,180]);
    colormap(hsv_nan);
    set(gcf, 'InvertHardCopy', 'off');
    set(gcf, 'Color', [1 1 1]);

    im_out = getframe(gcf);
    im_out = im_out.cdata;
    imwrite(im_out, fullfile([directory '/output'], ['angle_heatmap_' file(1:end-4) '.tif']));
    close

end

% calculate order parameter
av_ordermat = AFT_ordermat(anglemat, parameters);

end