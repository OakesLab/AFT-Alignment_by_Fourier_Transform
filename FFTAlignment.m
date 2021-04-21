%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code can be found at https://github.com/OakesLab/FFT_Alignment
% This routine produces a vectorfield of alignment directions using small
% subwindows in the real space image.  Images should be grayscale images.
% Angles determined are oriented as follows:
%
%                            ^  180?
%                            |
%                            |
%                            ------>  90?
%                            |
%                            |
%                            v  0?
%
% Order Parameter value:   0 = Completely random alignment
%                          1 = Perfectly aligned
%
% Inputs:
%       im          = the image to analyze. Can be any bit-depth. Can input
%                     as a matrix or as a filename
%       winsize     = Size of the sub windows. Must be odd. (Default is 33)
%       overlap     = How much the windows overlap. Value should be between 0
%                     (complete overlap) and 1 (no overlap). (Default is
%                     0.5)
%       st          = order parameter spacing. Considers a window of size
%                     2*st+1 to compare vectors to. Depending on the
%                     overlap size, you'll want to change this. (Default is
%                     2)
%       checkpoint  = threshold value above which to do the caculation. The
%                     sum of the image intensity in the window needs to be
%                     above this value for the routine to calculate a vector
%                     for the given window (Default is 0 - i.e. every
%                     window)
%       mask_method = Determines the method to mask the FFT when
%                     determining the moment calculation. Global (val = 1) uses
%                     a circle for every window. Local (val = 2) uses a
%                     local threshold for each subwindow. (Default = 1)
%       maskname    = Should be a logical mask the same size as im. Will
%                     only analyze regions where the mask is true. An
%                     empty input will analyze the entire image. Can input
%                     as a matrix or a filename
%       figures     = Plot figures. 0 = No; 1 = Yes; (Default = 1)
%
%
% Outputs:
%       Everything is stored in a structure file FFTAlignmentData
%       anglemat     = A matrix of the calculated alignment direction in degrees.
%       ordermat     = A matrix of the calculated order parameter;
%       pos          = row and column positions for each vector
%       vec          = vector direction components [coldisp rowdisp]
%       hist_bins    = bin values for the histogram of angles
%       hist_val     = # of angles in each bin
%       his_val_norm = normalized values of the histogram
%       parameters   = parameters used in the routine
%
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

function av_ordermat = FFTAlignment(file, directory, parameters)

% load image
im = imread(fullfile(directory, file));
im = im2double(im);

% for masking images to only analyze part of it
if parameters.mask_method == 0 % global
    immask = ones(size(im));
else % local
    immask = logical(imread(parameters.mask_name));
end

% unpack parameters for ease of access
winsize = parameters.winsize;
overlap = parameters.overlap;
st = parameters.st;
checkpoint = parameters.checkpoint;

% further parameters and initialisations
if mod(winsize,2) == 0
    winsize = winsize + 1; % make odd
end

im2 = zeros(size(im));

winrad = floor(winsize/2);
winspace = ceil(winsize*overlap);

% filter used to remove edge effects of transform
gauss_filter = fspecial('gaussian',[winsize winsize],winsize/4);

% mask to apply during moments calculation
r = 0.75 * winrad; % Kind of arbitrary choice. Just needs to go to zero near edges

% global masking method
if parameters.mask_method == 0 % global
    mask = fspecial('disk',r);
    mask = vertcat(zeros((winsize-size(mask,1))/2,size(mask,2)),mask,zeros((winsize-size(mask,1))/2,size(mask,2)));
    mask = horzcat(zeros(size(mask,1),(winsize-size(mask,2))/2),mask,zeros(size(mask,1),(winsize-size(mask,2))/2));
    mask = logical(mask);
end

x = repmat(1:winsize,winsize,1);
y = x';

sz = size(im);

% find points to sample with a window size of 2*winrad+1 and spacing of winspace
grid_row = winrad+1:winspace:sz(1)-winrad;
grid_col = winrad+1:winspace:sz(2)-winrad;
numRows = length(grid_row);
numCols = length(grid_col);

% counter
k = 0;

% create a matrix to hold all the data we're going to need
anglemat = zeros(numRows,numCols) * NaN;

% sets string format for verbose progress output
strForm = sprintf('%%.%dd',length(num2str(numRows)));
fprintf('\nStarting row ');

for i = 1:numRows
    
    % Display our count as we go through the image
    procStr = sprintf(strForm,i);
    fprintf(1,[procStr,' of ',num2str(numRows)]);
    
    for j = 1:numCols
        
        if immask(grid_row(i),grid_col(j)) == 1
            
            % Create small window in real space and apply FFT and Gaussian filter
            window = im(grid_row(i)-winrad:grid_row(i)+winrad,grid_col(j)-winrad:grid_col(j)+winrad);
            
            % Create a checkpoint for doing calculation if sum of signal is
            % insignificant (i.e. if it's just blackspace)
            if sum(window(:)) > checkpoint
                
                window_fft = fftshift(fft2(window.*gauss_filter));
                im2(grid_row(i)-winrad:grid_row(i)+winrad,grid_col(j)-winrad:grid_col(j)+winrad) = window_fft;
                k = k + 1;
                if sum(window_fft(:)) ~= 0
                    
                    window_fft = sqrt(window_fft.*conj(window_fft));
                    
                    % local masking method
                    if parameters.mask_method == 1
                        mask = zeros(size(window));
                        pts = window_fft>(2*mean(window_fft(:)));
                        mask(pts) = 1;
                        masklab = bwlabel(mask);
                        stats = regionprops(masklab,'area');
                        stats = [stats.Area];
                        val = find(stats == max(stats));
                        maskthresh = zeros(size(mask));
                        maskthresh(masklab==val) = 1;
                        
                        % Mask the fft for the moment calculations
                        win = log(window_fft).*maskthresh;
                    else
                        win = window_fft.*mask;
                    end
                    
                    % Calculate the various Moments
                    M00 = sum(win(:));
                    M10 = sum(sum(x.*win));
                    M01 = sum(sum(y.*win));
                    M11 = sum(sum(x.*y.*win));
                    M20 = sum(sum(x.*x.*win));
                    M02 = sum(sum(y.*y.*win));
                    
                    % The Center of Mass
                    xave = M10/M00;
                    yave = M01/M00;
                    
                    % Calculate the central moments
                    mu20 = M20/M00-xave^2;
                    mu02 = M02/M00-yave^2;
                    mu11 = M11/M00-xave*yave;
                    
                    % Angle of axis of the least second moment
                    theta = 0.5*atan(2*mu11/(mu20-mu02));
                    
                    % Convert angle to proper orientation for my frame of reference
                    if theta > 0 && theta < pi/4
                        theta = pi/2-theta;
                    end
                    if theta > -pi/4 && theta < 0
                        theta = -1*theta;
                    end
                    
                    % Find points to draw line scans
                    x2 = xave+r*cos(theta);
                    y2 = yave-r*sin(theta);
                    
                    x3 = xave-r*sin(theta);
                    y3 = yave-r*cos(theta);
                    
                    % Interpolate line scans
                    line1 = improfile(log(window_fft),[xave x2],[yave y2],r);
                    line2 = improfile(log(window_fft),[xave x3],[yave y3],r);
                    
                    % Determine which line is the maximum and minimum direction and correct theta accordingly
                    if sum(line1) > sum(line2)
                        maxline = line1;
                        minline = line2;
                        if x2 < xave
                            theta = theta + pi/2;
                        end
                    else
                        maxline = line2;
                        minline = line1;
                        if x3 < xave
                            theta = theta + pi/2;
                        end
                    end
                    
                    anglemat(i,j) = theta;
                    
                    
                else
                    anglemat(i,j) = NaN;
                end
                
                pr(k) = grid_row(i);
                pc(k) = grid_col(j);
                ur(k) = -1*floor(winrad/2)*cos(theta);
                vc(k) = -1*floor(winrad/2)*sin(theta);
                ur2(k) = floor(winrad/2)*cos(theta);
                vc2(k) = floor(winrad/2)*sin(theta);
                ang(k,1) = theta;
            end
        end
    end
    
    % make our verbose output overwrite the previous line
    if i <length(grid_row)
        for jj = 1:(length(procStr)+4+length(num2str(numRows)))
            fprintf(1,'\b');
        end
    end
end

% plot vectorfield
if parameters.figures == 1
    
    figure
    imshow(im, [])
    hold on
    quiver(pc,pr,vc2,ur2,0,'y','showarrowhead','off','linewidth',2)
    im_out = getframe(gcf);
    im_out = im_out.cdata;
    imwrite(im_out, fullfile(directory, ['vectors_' file(1:end-4) '.tif']));
    close
    
end

% Calculate the order parameter
ordermat = NaN(size(anglemat));
for i = st+1:size(anglemat,1)-st
    for j = st+1:size(anglemat,2)-st
        temp = anglemat(i-st:i+st,j-st:j+st);
        temp2 = repmat(anglemat(i,j),2*st+1,2*st+1);
        comp = cos(temp-temp2);
        comp = comp.*comp;
        ordermat(i,j) = 2*(nanmean(comp(:))-.5);
        
        clear temp temp2 comp
    end
end
av_ordermat = nanmedian(ordermat(:));

end
