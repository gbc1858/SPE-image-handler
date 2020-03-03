from spe_image_rms_rms_error import *
import settings


SPE_image = DenoiseSPEImage(settings.FILE_START_STRING, settings.FOLDER_PATH, settings.FILE_LIST)
SPE_image.set_cropped_range(100, 480, 100, 500)

# SPE_image.draw_beam_contour(100, 50)
x_rms, y_rms, x_std, y_std = SPE_image.get_rms_and_rms_error(contour_method=False, regular_method=True)
SPE_image.plot_single_frame(contour_method=False, regular_method=True)

print("x RMS size: ", x_rms)
print("x RMS STD: ", x_std)
print("y RMS size: ", y_rms)
print("y RMS STD: ", y_std)

