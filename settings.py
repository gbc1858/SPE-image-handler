import os


# Imaging SOL excel related
IMAGING_SOL_EXCEL_NAME = 'blue_center_ASTRA.xlsx'
IMAGING_SOL_EXCEL_SHEET = 'Sheet1'
IMAGING_SOL_EXCEL_Z_TEMP = 'z_off'
IMAGING_SOL_EXCEL_BZ_TEMP = 'B'

# BSOL file
BSOL_FILE_NAME = 'Blue_solenoid_calibration.xlsx'
BSOL_CURRENT_COLUMN = 'I_clamp_A'
BSOL_B_FIELD_COLUMN = 'Bz_Gauss'

# File Related
FILE_LIST = sorted(os.listdir("/Users/chen/Desktop/github/spe_image_handler/sample_data"))
FOLDER_PATH = "/Users/chen/Desktop/github/spe_image_handler/sample_data/"
FILE_START_STRING = 'set2'
SOLENOID_SETTING_IMG = 'set2_02'

# Solenoid scan related
SOL_SETTINGS = 11
NUM_OF_FRAMES = 20

# Image contour related
CONTOUR_LEVEL = 100