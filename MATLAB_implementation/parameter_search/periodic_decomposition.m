function [im_window_periodic, im_window_smooth] = periodic_decomposition(im)

% [Moisan, L (2011). Journal of Mathematical Imaging and Vision, 39(2), 161?179]

% find the number of rows and cols
[N_rows, N_cols] = size(im);

% create a zero matrix the size of the image
v = zeros(size(im));

% fill the edges of V with the difference between the opposite edge of the real image
v(1,:) = im(1,:) - im(end,:);
v(end,:) = -v(1,:);
v(:,1) = v(:,1) + im(:,1) - im(:,end);
v(:,end) = v(:,end) - im(:,1) + im(:,end);

% calculate the frequencies of the image
fx = repmat(cos(2 * pi() * linspace(0,N_cols-1,N_cols) / N_cols),N_rows,1);
fy = repmat(cos(2 * pi() * linspace(0,N_rows-1,N_rows) / N_rows),N_cols,1)';
% set the fx(1,1) to 0 to avoid division by zero
fx(1,1) = 0;
% calculate the smoothed image component
im_window_smooth = real(ifft2(fft2(v).* 0.5./(2 - fx - fy)));
% If you want to calculate the periodic fft directly
% p_fft = fftshift(fft2(actin) - fft2(v) * 0.5 / (2 - fx - fy))
im_window_periodic = im - im_window_smooth;
    
end





  