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

function anglemat = AFT_parameter_search_anglemat(im, winsize, winrad, winspace)

% set further parameters
checkpoint = 0; % threshhold sum in each window to do the calculation
mask_method = 0;	% global = 0; local = 1;

im2 = zeros(size(im));

% mask to apply during moments calculation
r = 0.5 * winrad; 

% global masking method
if mask_method == 0 % global
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

        % Create small window in real space and apply FFT and Gaussian filter
        window = im(grid_row(i)-winrad:grid_row(i)+winrad,grid_col(j)-winrad:grid_col(j)+winrad);
        
        % Create a checkpoint for doing calculation if sum of signal is
        % insignificant (i.e. if it's just blackspace)
        if mean(window(:)) > checkpoint
            
            % separate out the periodic and smooth components
            [im_window_periodic, ~] = periodic_decomposition(window);
            % take the FFT of the periodic component
            window_fft = fftshift(fft2(im_window_periodic));
                
            im2(grid_row(i)-winrad:grid_row(i)+winrad,grid_col(j)-winrad:grid_col(j)+winrad) = window_fft;
            k = k + 1;
            if sum(window_fft(:)) ~= 0
                
                window_fft = sqrt(window_fft.*conj(window_fft));
                win = window_fft.*mask;

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
                theta = 0.5*atan2(2*mu11,(mu20-mu02));
                
                % Convert angle to proper orientation for my frame of reference
                if theta > 0 && theta < pi/2
                    theta = pi-theta;
                end
                if theta > -pi/2 && theta < 0
                    theta = -1*theta;
                end
                
                anglemat(i,j) = theta;
                
            else
                anglemat(i,j) = NaN;
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


end
