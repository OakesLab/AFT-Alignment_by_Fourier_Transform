# AFT - Alignment by Fourier Transform 
## Python implementation

For making images and movies from microscopy images in a Jupyter notebook. Requires the following packages (Version that the approach has been validated with is included in paratheses):

`numpy (1.20.3)`
`scikit-image (0.18.2)`
`matplotlib (3.4.2)`
`scipy (1.7.0)`
`opencv-python (4.1.1.26)`
`pandas (1.3.0)`


# Example Usage

Jupyter notebookd with examples are included in the repository. Example data can be found here (LINK TO DATA): 

The general approach is outlined below.

Import the AFT alignment tools
```python
import AFT_tools as AFT
```

Define your parameters
```python
# AFT parameters

#### required parameters ####
window_size = 65
overlap = 0.5
neighborhood_radius = 5

#### optional parameters ####
intensity_threshold = 0
eccentricity_threshold = 1
im_mask = io.imread('mask').astype('bool')

#### output parameters ####
plot_overlay = True
plot_angles = True
plot_eccentricity = True
save_figures = True
data_save_path = 'output_data/'    
```

Read in your image
```python
im = io.imread('actin/anisotropic/im1.tif')
```

Measure the local alignment
```python
x, y, u, v, im_theta, im_eccentricity = AFT.image_local_order(im, window_size, overlap, save_path = data_save_path,
                                                             plot_overlay=plot_overlay, plot_angles=plot_angles, 
                                                             plot_eccentricity=plot_eccentricity,
                                                             save_figures=data_save_path)
```

Calculate the average order parameter for the image
```python
im_order_parameter = AFT.calculate_order_parameter(im_theta, neighborhood_radius)
print('The order parameter is {0}'.format(im_order_parameter))
```

Additional examples, including doing a parameter search and measuring the order parameter for an image stack can be found in the Jupyter notebooks in the repository
