import numpy as np
from scipy.fft import fft2, fftshift, ifft2, ifftshift, fft, ifft, fftfreq
import skimage.io as io
import cv2                                                     # for filtering vector fields
from skimage.morphology import disk        # morphology operations
import numpy.matlib as matlib


def image_norm(im):
    im_norm = np.sqrt(np.real(im * np.conj(im)))
    return im_norm

def periodic_decomposition(im):
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
    
    # correct for real space
    theta = theta + np.pi/2

    return theta, eccentricity

def image_local_order(im, window_size = 33, overlap = 0.5):
    # get the image shape
    N_rows, N_cols = im.shape
    
    # define the radius of the window
    radius = int(np.floor((window_size) / 2))
    
    # make a list of the r,c positions for the windows
    rpos = np.arange(radius,N_rows-radius,int(window_size * overlap))
    cpos = np.arange(radius,N_cols-radius,int(window_size * overlap))

    # make a structuring element to filter the mask
    bpass_filter = disk(radius)

    # make mask
    im_mask = np.zeros((window_size, window_size))
    im_mask[int(np.floor(window_size/2)), int(np.floor(window_size/2))] = 1

    # filter the mask with the structuring element to define the ROI
    im_mask = cv2.filter2D(im_mask, -1, bpass_filter)
    im_mask = np.rint(im_mask) == 1

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
            # define the window to analyze
            im_window = im[r-radius:r+radius+1,c-radius:c+radius+1]
            # separate out the periodic and smooth components
            im_window_periodic, im_window_smooth = periodic_decomposition(im_window)
            # take the FFT of the periodic component
            im_window_fft = fftshift(fft2(im_window_periodic))
            # find the image norm and mulitply by the mask
            im_window_fft_norm = image_norm(im_window_fft) * im_mask
            # calculate the angle and eccentricity of orientation based on the FFT moments
            theta, eccentricity = least_moment(im_window_fft_norm, xcoords, ycoords)
            # add the values to each list
            im_theta.append(theta)
            im_ecc.append(eccentricity)
            x.append(c)
            y.append(r)
            u.append(np.cos(theta) * arrow_length)
            v.append(np.sin(theta) * arrow_length)

    # turn all the lists into arrays for simple indexing
    x = np.array(x)
    y = np.array(y)
    u = np.array(u)
    v = np.array(v)
    im_theta = np.array(im_theta)
    im_ecc = np.array(im_ecc)

    return x, y, u, v, im_theta, im_ecc