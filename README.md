# FFT Alignment

Algorithms to calculate fibre alignment from biological images can be found in this repository. Please refer to the documentation below for further information.

This analysis was produced for [Cetera et al. Nat Commun 2014; 5:1-12](http://www.ncbi.nlm.nih.gov/pubmed/25413675). A methods paper is currently in preparation,
link and citation will be added in due course.

## FFTAlignment_batch.m
This routine uses a vector field of alignment directions using small sub-windows in the real space image to calculate an alignment order parameter. It calls the function `FFTAlignment.m`
  * _Input_: a folder containing grayscale images in .tif format
  * _Output_: an array containing an alignment order parameter averaged across each image
  * _Parameters_:
    * `Window size [px]` = Size of the sub-windows in pixel
    * `Window overlap [%]` = Percentage of window overlap
    * `Neighbourhood size [vectors]` = Neighbourhood size on which the order parameter is calculated
    * `Masking method (local = 1, global = 0)` = Determines the method to mask the FFT when determining the moment calculation. Global (val = 0) uses a circle for every window. Local (val = 1) uses a local threshold for each sub-window. If local, the user will be asked to input a second folder containing logical masks of the regions of each image that require analysis.
    * `Display figures (yes = 1, no = 0)` = Option to save vector field overlaid on original input images

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
