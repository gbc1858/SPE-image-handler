# SPE File Handler

## Description

This file is to analyze the .SPE image. Horizontal and vertical RMS beam spot sizes will be calculated separately.

## Usage
There are two methods included. Following is a brief description of how each method works. Both methods suggest a similar 
result.
 
 - Method#1
    - **Crop all images.** All beam images will be cropped according to the user provided parameters. 
    ```python
    SPE_image = DenoiseSPEImage(settings.FILE_START_STRING, settings.FOLDER_PATH, settings.FILE_LIST)
    SPE_image.set_cropped_range(100, 480, 200, 500)
    ```
    - **Generate beam contours.** To outline the beam spot in every frame, contours will be generated based on the user 
    provided contour levels. The contour with the longest path is the beam spot contour.
    - **Denoise the beam image.** All pixel datapoints outside the beam contour will be set to zero. 
    - **Calculate the background within the contour.** Background will be calculated from the lowest 100 pixel values 
    within the contour. And subtract the calculated background from the copped image file. 
    - **Denoise the beam image (w/o background) the second time.** 
    - **Calculate the beam RMS sizes.** Beam RMS sizes are calculated from the Gaussian fitting.
    ```python
    x_rms, y_rms, x_std, y_std = SPE_image.get_rms_and_rms_error(contour_method=True)
    ```
 - Method#2 (a bit faster)
    - **Crop all images.**
    ```python
    SPE_image.set_cropped_range(100, 480, 200, 500)
    ```
    - **Calculate the background.** Background will be calculated from the upper corner of the image profile. The 
    calculated background will be subtracted from the image file.
    - **Denoise the beam image (w/o background).**
    - **Calculate the beam RMS sizes.**
    ```python
    x_rms, y_rms, x_std, y_std = SPE_image.get_rms_and_rms_error(regular_method=True)
    ```


## TODOs
- [ ] Auto detect the cropping region.
