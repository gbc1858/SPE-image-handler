# SPE File Handler

## Description

This file is to analyze the .SPE image. Horizontal and vertical RMS beam spot sizes will be calculated separately.

## Usage
For the beam RMS size calculation, there are two methods included. Following is a brief description of how each method 
works. Both methods suggest a similar result.
 
 - #### Method#1
    - **Crop all images.** All beam images will be cropped according to the user provided parameters. 
    ```python
   from spe_image_rms_rms_error import *
   
   SPE_image = DenoiseSPEImage(settings.FILE_START_STRING, settings.FOLDER_PATH, settings.FILE_LIST)
   SPE_image.set_cropped_range(100, 480, 200, 500)
    ```
    - **Generate beam contours.** To outline the beam spot (or to find the ROI) of every frame, contours will be 
    generated based on the user provided contour levels. The contour with the longest path is the beam spot contour. 
    Beam images and the corresponding contours can be displayed using the function of `draw_beam_contour()`.
    ```python
   SPE_image.draw_beam_contour(100, 50)
    ```
    - **Denoise the beam image.** All pixel datapoints outside the beam contour will be set to zero. 
    - **Calculate the background within the contour.** Background will be calculated from the lowest 100 pixel values 
    within the contour. And subtract the calculated background from the copped image file. 
    - **Denoise the beam image (w/o background) the second time.** Use `scipy.ndimage.median_filter()` function.
    - **Calculate the beam RMS sizes.** Beam RMS sizes are calculated from the Gaussian fitting.
    ```python
   x_rms, y_rms, x_std, y_std = SPE_image.get_rms_and_rms_error(contour_method=True)
    ```
 - #### Method#2 (a bit faster)
    - **Crop all images.**
    ```python
   SPE_image.set_cropped_range(100, 480, 200, 500)
    ```
    - **Calculate the background.** Background will be calculated from the upper corner of the image profile. The 
    calculated background will be subtracted from the image file.
    - **Denoise the beam image (w/o background).** Use `scipy.ndimage.median_filter()` function.
    - **Calculate the beam RMS sizes.**
    ```python
   x_rms, y_rms, x_std, y_std = SPE_image.get_rms_and_rms_error(regular_method=True)
    ```
Beam profiles can be saved to one or multiple pdf files using the function of `plot_single_frame()`. Note: the two Boolean 
arguments need to be consistent to the used method.
```python
SPE_image.plot_single_frame(contour_method=True, regular_method=False)    # for method#1
SPE_image.plot_single_frame(contour_method=False, regular_method=True)    # for method#2
```
Run `test.py` for a simple test.

## TODOs
- [ ] Auto detect the cropping region.
- [x] Fix the frame labeling bug.

## Credit
The enclosed `winspec.py` was developed by Anton Loukianov (@antonl). Please visit https://github.com/antonl/pyWinSpec for more 
details.
