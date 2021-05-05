# FFT Alignment

Algorithms to calculate fibre alignment from biological images can be found in this repository. Please refer to the documentation below for further information.

This analysis was produced for [Cetera et al. Nat Commun 2014; 5:1-12](http://www.ncbi.nlm.nih.gov/pubmed/25413675). A methods paper is currently in preparation, link and citation will be added in due course.

## FFTAlignment_batch.m
This routine uses a vector field of alignment directions using small sub-windows in the real space image to calculate an alignment order parameter. It calls the functions `FFTAlignment.m` and `user_input.m`.
  * _Input_: a folder containing grayscale images in .tif format
  * _Output_: an array containing an alignment order parameter averaged across each image; analysis parameters; input images overlaid with angle vector field (optional)
  * _Parameters_:
    * `Window size [px]` = Size of the sub-windows in pixel
    * `Window overlap [%]` = Percentage of window overlap
    * `Neighbourhood size [vectors]` = Neighbourhood size on which the order parameter is calculated
  * _Options_:
    * `Save output images` = Option to save vector field overlaid on the original input images
    * `Apply local masking and/or filtering on the images` = Option to apply local masking to part of the input image, filter out blank spaces and/or isotropic regions.
      * `Masking method (local = 1, global = 0)` = Determines the method to mask the FFT when determining the moment calculation. Global (val = 0) uses a circle for every window. Local (val = 1) uses a local threshold for each sub-window. If local, the user will be asked to input a second folder containing logical masks of the regions of each image that require analysis.
      * `Ignore blank spaces (yes = 1, no = 0)` = Option to filter out regions of the image based on their mean intensity value; the user is requested to declare a minimum threshold value between 0 (black) and 255 (white). If not set, all regions in the image are analysed (set to threshold value = 0 by default).
      * `Ignore isotropic regions (yes = 1, no = 0)` = Option to filter out regions of the image based on the isotropy of the signal (coherency of the FFT); the user is requested to declare a minimum threshold value between 0 (isotropic) and 1 (highly aligned). If not set, all regions in the image are analysed (set to threshold value = 0 by default).

Angles determined are oriented as follows:

                            ^  180°
                            |
                            |
                            ------>  90°
                            |
                            |
                            v  0°

Order parameter values:     
  * 0 = Completely random alignment
  * 1 = Perfectly aligned


## FFTAlignment_batch_parameter_search.m
This routine is based on `FFTAlignment_batch.m` and can be used to search for a parameter set where the difference in alignment between two sample populations is more pronounced. The order parameter is calculated for a range of window and neighbourhood sizes; overlap is set to 50% and the masking method is set to global (the whole image is analysed). The median order parameter calculated for each permutation of window and neighbourhood sizes is compared between the two population, either by looking at their difference (sample 1 - sample 2) or the p-value of a non-parameteric statistical comparison (Mann-Whitney test). It calls the functions `FFTAlignment_parameter_search_main.m`, `FFTAlignment_parameter_search_anglemat`, `FFTAlignment_parameter_search_ordermat`.

* _Input_: two separate folders containing grayscale images in .tif format for the 2 samples to be compared
* _Output_: a cell array containing an alignment order parameter averaged across each image for each window and neighbourhood sizes; analysis parameters; heatmaps for order parameter comparison (difference and p-value)
* _Parameters_:
  * `Minimum window size [px]` = Size of the smallest sub-window to be tested (in pixel)
  * `Window size interval [px]` = Interval between subsequent sub-window sizes to be tested (in pixel)
