import numpy as np
from scipy.fft import fft2, fftshift, ifft2, ifftshift, fft, ifft, fftfreq
import skimage.io as io
import cv2                                                     # for filtering vector fields
from skimage.morphology import disk        # morphology operations
import numpy.matlib as matlib
import matplotlib.pyplot as plt


def image_norm(im):
    im_norm = np.sqrt(np.real(im * np.conj(im)))
    return im_norm

def periodic_decomposition(im):
    im = im.astype('float32')
    # find the number of rows and cols
    N_rows, N_cols = im.shape
    # create an zero matrix the size of the image
    v = np.zeros((N_rows,N_cols))
    # fill the edges of V with the difference between the opposite edge of the real image
    v[0,:] = im[0,:] - im[-1,:]
    v[-1,:] = -v[0,:]
    v[:,0] = v[:,0] + im[:,0] - im[:,-1]
    v[:,-1] = v[:,-1] - im[:,0] + im[:,-1]
    # calculate the frequencies of the image
    fx = matlib.repmat(np.cos(2 * np.pi * np.arange(0,N_cols) / N_cols),N_rows,1)
    fy = matlib.repmat(np.cos(2 * np.pi * np.arange(0,N_rows) / N_rows),N_cols,1).T
    # set the fx[0,0] to 0 to avoid division by zero
    fx[0,0] = 0
    # calculate the smoothed image component
    s = np.real(ifft2(fft2(v) * 0.5 / (2 - fx - fy)))
    # If you want to calculate the periodic fft directly
    # p_fft = fftshift(fft2(actin) - fft2(v) * 0.5 / (2 - fx - fy))
    p = im - s

    return p, s

def least_moment(image, xcoords=[], ycoords=[]):
    # get the image shape
    N_rows, N_cols = image.shape

    # check if xcoords and ycoords are passed in the function
    if len(xcoords) == 0:
        # create coordinates for the image
        xcoords, ycoords = np.meshgrid(np.arange(0,N_cols) , np.arange(0,N_rows))

    #calculate the moments
    M00 = np.sum(image)
    M10 = np.sum(image * xcoords)
    M01 = np.sum(image * ycoords)
    M11 = np.sum(image * xcoords * ycoords)
    M20 = np.sum(image * xcoords * xcoords)
    M02 = np.sum(image * ycoords * ycoords)

    # center of mass
    xave = M10 / M00
    yave = M01 / M00

    # calculate the central moments
    mu20 = M20/M00 - xave**2
    mu02 = M02/M00 - yave**2
    mu11 = M11/M00 - xave*yave

    # angle of axis
    theta = 0.5 * np.arctan2((2 * mu11),(mu20 - mu02))

    # multiply by -1 to correct for origin being in top left corner instead of bottom right
    theta = -1 * theta

    # find eigenvectors
    lambda1 = (0.5 * (mu20 + mu02)) + (0.5 * np.sqrt(4 * mu11**2 + (mu20 - mu02)**2))
    lambda2 = (0.5 * (mu20 + mu02)) - (0.5 * np.sqrt(4 * mu11**2 + (mu20 - mu02)**2))
    
    # calculate the eccentricity (e.g. how oblong it is)
    eccentricity = np.sqrt(1 - lambda2/lambda1)

    return theta, eccentricity

def image_local_order(im, window_size = 33, overlap = 0.5, im_mask = None, intensity_thresh = 0, eccentricitiy_thresh = 1, plot_angles=True, plot_eccentricity=True):
    # get the image shape
    N_rows, N_cols = im.shape

    # make window size off if it isn't already
    if window_size % 2 == 0:
        window_size += 1
    
    # define the radius of the window
    radius = int(np.floor((window_size) / 2))
    
    # make a list of the r,c positions for the windows
    rpos = np.arange(radius,N_rows-radius,int(window_size * overlap))
    cpos = np.arange(radius,N_cols-radius,int(window_size * overlap))

    # make a structuring element to filter the mask
    bpass_filter = disk(radius * .5)

    # make window mask
    window_mask = np.zeros((window_size, window_size))
    window_mask[int(np.floor(window_size/2)), int(np.floor(window_size/2))] = 1

    # filter the mask with the structuring element to define the ROI
    window_mask = cv2.filter2D(window_mask, -1, bpass_filter)
    window_mask = np.rint(window_mask) == 1

    # check if there is an image mask
    if im_mask is None:
        im_mask = np.ones_like(im)
    else:
        # make sure the input mask is a boolean
        im_mask = im_mask.astype('bool')

    # make x and y coordinate matrices
    xcoords, ycoords = np.meshgrid(np.arange(0,window_size) , np.arange(0,window_size))

    # make empty list variables
    im_theta, im_ecc = [], []
    x, y, u, v = [], [], [], []
    
    # length of orientation vector
    arrow_length = radius / 2
    
    # loop through each position and measure the local orientation
    for r in rpos:
        for c in cpos:

            # check to see if point is within image mask
            if im_mask[r,c] == True:
                # define the window to analyze
                im_window = im[r-radius:r+radius+1,c-radius:c+radius+1]

                # check that it's above the intensity threshold
                if np.mean(im_window) > intensity_thresh:
                    # separate out the periodic and smooth components
                    im_window_periodic, im_window_smooth = periodic_decomposition(im_window)
                    # take the FFT of the periodic component
                    im_window_fft = fftshift(fft2(im_window_periodic))
                    # find the image norm and mulitply by the mask
                    im_window_fft_norm = image_norm(im_window_fft) * window_mask
                    # calculate the angle and eccentricity of orientation based on the FFT moments
                    theta, eccentricity = least_moment(im_window_fft_norm, xcoords, ycoords)

                    # correct for real space
                    theta = theta + np.pi/2

                    # map everything back to between -pi/2 and pi/2
                    if theta > np.pi/2:
                        theta -= np.pi

                    # filter based on eccentricity
                    if eccentricity > eccentricitiy_thresh:
                        eccentricity = np.nan
                        theta = np.nan

                    # add the values to each list
                    im_theta.append(theta)
                    im_ecc.append(eccentricity)
                    x.append(c)
                    y.append(r)
                    u.append(np.cos(theta) * arrow_length)
                    v.append(np.sin(theta) * arrow_length)
                else:
                    im_theta.append(np.nan)
                    im_ecc.append(np.nan)
                    x.append(c)
                    y.append(r)
                    u.append(np.nan)
                    v.append(np.nan)
            else:
                im_theta.append(np.nan)
                im_ecc.append(np.nan)
                x.append(c)
                y.append(r)
                u.append(np.nan)
                v.append(np.nan)

    # turn all the lists into arrays for simple indexing
    x = np.array(x)
    y = np.array(y)
    u = np.array(u)
    v = np.array(v)
    im_theta = np.reshape(im_theta,(len(rpos),len(cpos)))
    im_ecc = np.reshape(im_ecc,(len(rpos),len(cpos)))

    if plot_angles:
        plt.figure()
        plt.imshow(im_theta, vmin=-np.pi/2, vmax=np.pi/2, cmap='hsv')
        plt.colorbar()
        plt.show()

    if plot_eccentricity:
        plt.figure()
        plt.imshow(im_ecc, vmin=0, vmax=1)
        plt.colorbar()
        plt.show()



    return x, y, u, v, im_theta, im_ecc


def calculate_order_parameter(im_theta, neighborhood_radius):

    order_list = []

    rpos = np.arange(neighborhood_radius,im_theta.shape[0]-neighborhood_radius)
    cpos = np.arange(neighborhood_radius,im_theta.shape[1]-neighborhood_radius)

    for r in rpos:
        for c in cpos:
            search_window = im_theta[r-neighborhood_radius:r+neighborhood_radius+1,c-neighborhood_radius:c+neighborhood_radius+1]
            vector_array = np.ones_like(search_window) * im_theta[r,c]
            order_array = np.cos(search_window - vector_array) ** 2 - 0.5
            order_array = np.delete(order_array, (2*neighborhood_radius + 1)**2 // 2)
            order_list.append(2 * np.nanmean(order_array))

    im_orderparameter = np.nanmedian(order_list)

    return im_orderparameter
