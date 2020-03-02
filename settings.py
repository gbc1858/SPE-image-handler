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
FILE_LIST = sorted(os.listdir("/Users/your_data_path"))
FOLDER_PATH = "/Users/your_data_path/"
FILE_START_STRING = 'set2'

# Solenoid scan related
SOL_SETTINGS = 11
NUM_OF_FRAMES = 20

# Contour related
CONTOUR_LEVEL = 100