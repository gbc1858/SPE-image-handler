import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from numpy import mean, std
from scipy import ndimage as ndi
from scipy.optimize import curve_fit
import scipy.ndimage
import itertools as it
from math import floor

import constant
import settings
from Exceptions import *
from winspec import SpeFile


def quadratic(data, aa, bb, cc):
    return aa * data ** 2 + bb * data + cc


def get_bsol_field(i_sol):
    """
    :param i_sol: current of the blue solenoid [A]
    :return: a magnetic field list of blue solenoid [Gauss]
    """
    bsol_file = pd.read_excel(settings.BSOL_FILE_NAME)
    current = bsol_file[settings.BSOL_CURRENT_COLUMN]
    b_field = bsol_file[settings.BSOL_B_FIELD_COLUMN]
    a, b = np.polyfit(current, b_field, 1)
    return [i * a + b for i in i_sol]


# Define Gaussian function with offset
def gaus(x, a, x0, sigma, c):
    return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2)) + c


def calculate_gauss_fitting_params(interval, intensity_ave_no_bg):
    mean = sum(intensity_ave_no_bg * interval) / sum(intensity_ave_no_bg)
    sigma = np.sqrt(np.absolute(sum(intensity_ave_no_bg * (interval - mean) ** 2) / sum(intensity_ave_no_bg)))
    return curve_fit(gaus, interval, intensity_ave_no_bg, p0=[1., mean, sigma, 0], maxfev=10000000)


def calculate_intensity_avg_no_bg(bg, intensity_ave):
    intensity_ave_no_bg = [i - bg for i in intensity_ave]
    for index in range(len(intensity_ave_no_bg)):
        intensity_ave_no_bg[index] = 0 if intensity_ave_no_bg[index] < 0 else intensity_ave_no_bg[index]
    return intensity_ave_no_bg


def convert_cont_to_current(count: float):
    """
    :param count: solenoid count number
    :return: single number of a current
    """
    if 0 <= count <= 1312:
        count_ls = [0, 608, 672, 736, 800, 864, 928, 992, 1056, 1120, 1184, 1248, 1312]  # blue solenoid count number
        current_ls = [0, 11.7, 13.1, 14.4, 15.6, 17.1, 18.2, 19.6, 20.8, 22, 23.3, 24.7, 26]  # blue solenoid current
        xvals = np.linspace(0, max(count_ls), max(count_ls) + 1)
        yinterp = np.interp(xvals, count_ls, current_ls, 1)
        return yinterp[count]
    else:
        raise CurrentValueError('Exceed the range of the interpolation. Calculation is terminated. '
                                'Please double check the solenoid count number.')


class DenoiseSPEImage:
    def __init__(self, file_start_str, folder_location, file_list):
        self.file_start_str = file_start_str
        self.folder_location = folder_location
        self.file_list = file_list
        self.spe_file_list = []
        self.sol_cont = []
        self.x_rms_all = []
        self.y_rms_all = []
        self.x_rms_std = []
        self.y_rms_std = []
        self.contour_paths_all = []
        self.zoom_in_single_frame = []
        self.x_intensity_no_bg = []
        self.y_intensity_no_bg = []
        self.x_start = None
        self.x_end = None
        self.y_start = None
        self.y_end = None
        self.bg_x_start = None
        self.bg_x_end = None
        self.bg_y_start = None
        self.bg_y_end = None
        self.set_background_range(0, 8, 0, 8)
        for file in self.file_list:
            if file.endswith(".SPE") and file.startswith(self.file_start_str):
                self.spe_file_list.append(file)

    def clear(self):
        self.sol_cont = []
        self.x_rms_all = []
        self.y_rms_all = []
        self.x_rms_std = []
        self.y_rms_std = []

    # def get_all_rms_and_std(self):
    #     return self.x_rms_all, self.y_rms_all, self.x_rms_std, self.y_rms_std

    def set_cropped_range(self, x_start, x_end, y_start, y_end):
        self.x_start = x_start
        self.x_end = x_end
        self.y_start = y_start
        self.y_end = y_end

    def has_cropped_range(self):
        if self.x_start and self.x_end and self.y_end and self.y_start:
            return True
        return False

    def set_background_range(self, bg_x_start, bg_x_end, bg_y_start, bg_y_end):
        self.bg_x_start = bg_x_start
        self.bg_x_end = bg_x_end
        self.bg_y_start = bg_y_start
        self.bg_y_end = bg_y_end

    def has_background_range(self):
        if self.bg_x_start is not None and self.bg_x_end is not None and self.y_start is not None and self.y_end is not None:
            return True
        return False

    def get_background(self, background_arr):
        if not self.has_background_range():
            raise DenoiseBGError('You need to set the background range first.')
        cropped_background = background_arr[self.bg_y_start:self.bg_y_end, self.bg_x_start:self.bg_x_end]
        return np.mean(cropped_background.flatten())

    def has_bg_within_contour(self):
        if self.zoom_in_single_frame:
            return True
        return False

    def get_bg_within_contour(self, current_frame):
        # if not self.has_bg_within_contour():
        #     raise DenoiseFrameAfterContourError("Contours need to be added to the images, and force zero outside the "
        #                                         "contour region.")

        single_frame = self.get_main_beam_contour_and_force_outer_zero(current_frame, settings.CONTOUR_LEVEL)
        single_frame_temp = sorted(single_frame.flatten())
        single_frame_temp[:] = (i for i in single_frame_temp if i != 0)
        return np.mean(single_frame_temp[:100])

    def get_current_all_frame(self, file):
        self.sol_cont.append(file.split('_')[1])
        return SpeFile(self.folder_location + file).data

    def save_pdf(self, file, file_type):
        if file_type == 'all':
            file_end_str = '_all_jet.pdf'
        elif file_type == 'single':
            file_end_str = '_jet_single_frame.pdf'
        else:
            raise ValueError('file_type has to be single or all.')
        pdf = PdfPages(file + file_end_str)
        for fig_index in range(1, plt.gcf().number + 1):
            pdf.savefig(fig_index)
            # plt.close()
        pdf.close()

    def generate_image_and_save(self, zoom_in, x_intensity_ave_no_bg, y_intensity_ave_no_bg, file, file_type,
                                single_frame_counter=None):
        x = np.linspace(0, len(x_intensity_ave_no_bg), len(x_intensity_ave_no_bg))
        y = np.linspace(0, len(y_intensity_ave_no_bg), len(y_intensity_ave_no_bg))
        fig = plt.figure()
        gs = gridspec.GridSpec(4, 4)
        ax_main = plt.subplot(gs[1:4, 0:3])

        ax_x_dist = plt.subplot(gs[0, 0:3], sharex=ax_main)
        ax_y_dist = plt.subplot(gs[1:4, 3], sharey=ax_main)

        ax_main.imshow(zoom_in, aspect='auto', cmap='jet')
        ax_main.axis([0, len(zoom_in[0]), len(zoom_in), 0])
        ax_main.set(xlabel="x pixel", ylabel="y pixel")

        x_intensity_ave_no_bg_normalized = [i * 10 / max(x_intensity_ave_no_bg) for i in x_intensity_ave_no_bg]
        y_intensity_ave_no_bg_normalized = [i * 10 / max(y_intensity_ave_no_bg) for i in y_intensity_ave_no_bg]
        ax_x_dist.bar(x, x_intensity_ave_no_bg_normalized, 1, color='b')
        ax_x_dist.set(ylabel='count')
        ax_x_dist.grid(ls=':', alpha=0.5)
        plt.setp(ax_x_dist.get_xticklabels(), visible=False)

        ax_y_dist.barh(y, y_intensity_ave_no_bg_normalized, 1, color='b')
        ax_y_dist.set(xlabel='count')
        ax_y_dist.grid(ls=':', alpha=0.5)
        plt.setp(ax_y_dist.get_yticklabels(), visible=False)

        # Gaussian fit
        popt, pcov = calculate_gauss_fitting_params(x, x_intensity_ave_no_bg_normalized)
        poptt, pcovv = calculate_gauss_fitting_params(y, y_intensity_ave_no_bg_normalized)

        ax_y_dist.plot(gaus(y, *poptt), y, label='gaussian fit', c='C03', alpha=0.6)
        ax_y_dist.text(5, (self.y_end - self.y_start) / 2 - 60,
                       'y$_{RMS}$=%.4f mm' % (abs(poptt[2]) * constant.pixel_res),
                       rotation=270, fontsize=8)
        ax_x_dist.plot(x, gaus(x, *popt), label='gaussian fit', c='C03', alpha=0.6)
        ax_x_dist.text((self.x_end - self.x_start) / 2 + 60, 5,
                       'x$_{RMS}$=%.4f mm' % (abs(popt[2]) * constant.pixel_res),
                       fontsize=8)

        self.x_rms_all.append(abs(popt[2] * constant.pixel_res))
        self.y_rms_all.append(abs(poptt[2]) * constant.pixel_res)

        if file_type == 'all':
            fig.suptitle(file, fontsize=14)
        elif file_type == 'single':
            if single_frame_counter is None:
                raise ValueError('Single frame counter is needed.')
            fig.suptitle(file + "\n frame#%d" % (single_frame_counter + 1), fontsize=14)
        else:
            raise ValueError('file_type has to be single or all.')
        print("\rSaving the beam profile of %s Frame#%i" % (file, single_frame_counter + 1), end="")
        self.save_pdf(file.split(".")[0], file_type)

    def get_intensity_ave_no_bg(self, frame):
        zoomed_in_current = frame[self.y_start:self.y_end, self.x_start:self.x_end]
        zoomed_in_current = ndi.median_filter(zoomed_in_current, 3)
        bg = self.get_background(zoomed_in_current)
        x_intensity_ave = zoomed_in_current.mean(axis=0).tolist()
        y_intensity_ave = zoomed_in_current.mean(axis=1).tolist()
        return calculate_intensity_avg_no_bg(bg, x_intensity_ave), calculate_intensity_avg_no_bg(bg, y_intensity_ave)

    def get_intensity_ave_using_contour_bg(self, frame):
        zoomed_single_frame = frame[self.y_start:self.y_end, self.x_start:self.x_end]
        bg = self.get_bg_within_contour(frame)
        zoomed_single_frame = np.array(zoomed_single_frame) - bg
        zoomed_single_frame[zoomed_single_frame < 0] = 0

        zoomed_single_frame = ndi.median_filter(zoomed_single_frame, 3)
        x_intensity_ave = zoomed_single_frame.mean(axis=0).tolist()
        y_intensity_ave = zoomed_single_frame.mean(axis=1).tolist()

        return x_intensity_ave, y_intensity_ave, zoomed_single_frame

    def get_rms_plot_all_ave_frames(self):
        """
        Beam image at each solenoid setting is averaged from the 20 frames. Then performing the RMS calculation and
        ploting etc.z
        """
        if not self.has_cropped_range():
            raise DenoiseCroppedError('You need to set the cropped range first.')
        self.clear()

        for file in self.spe_file_list:
            curr_all_frame = self.get_current_all_frame(file)
            results = sum(curr_all_frame[i] for i in range(len(curr_all_frame))) / len(curr_all_frame)

            zoomed_in_results = results[self.y_start:self.y_end, self.x_start:self.x_end]
            background = self.get_background(zoomed_in_results)

            x_intensity_ave = zoomed_in_results.mean(axis=0).tolist()
            x_intensity_ave_no_bg = [i - background for i in x_intensity_ave]
            y_intensity_ave = zoomed_in_results.mean(axis=1).tolist()
            y_intensity_ave_no_bg = [i - background for i in y_intensity_ave]
            self.generate_image_and_save(zoomed_in_results, x_intensity_ave_no_bg, y_intensity_ave_no_bg, file, 'all')

    def get_rms_and_rms_error(self, contour_method=False, regular_method=False):
        if not self.has_cropped_range():
            raise DenoiseCroppedError('You need to set the cropped range first.')
        self.clear()

        if contour_method == True and regular_method == True:
            raise RmsMethodError('Only one method is allowed to be True at once.')

        for file in self.spe_file_list:
            curr_all_frame = self.get_current_all_frame(file)
            x_rms_temp = []
            y_rms_temp = []
            for i in range(len(curr_all_frame)):
                if regular_method:
                    x_intensity_ave_no_bg, y_intensity_ave_no_bg = self.get_intensity_ave_no_bg(curr_all_frame[i])
                if contour_method:
                    x_intensity_ave_no_bg, y_intensity_ave_no_bg, zoomed_single_frame = self.get_intensity_ave_using_contour_bg(curr_all_frame[i])

                x = np.linspace(0, len(x_intensity_ave_no_bg), len(x_intensity_ave_no_bg))
                y = np.linspace(0, len(y_intensity_ave_no_bg), len(y_intensity_ave_no_bg))

                # Gaussian fit
                popt, pcov = calculate_gauss_fitting_params(x, x_intensity_ave_no_bg)
                poptt, pcovv = calculate_gauss_fitting_params(y, y_intensity_ave_no_bg)

                x_rms_temp.append(abs(popt[2] * constant.pixel_res))
                y_rms_temp.append(abs(poptt[2]) * constant.pixel_res)

            self.x_rms_std.append(std(x_rms_temp))
            self.y_rms_std.append(std(y_rms_temp))

            self.x_rms_all.append(mean(x_rms_temp))
            self.y_rms_all.append(mean(y_rms_temp))
        return self.x_rms_all, self.y_rms_all, self.x_rms_std, self.y_rms_std

    def plot_single_frame(self, contour_method=False, regular_method=False):
        if not self.has_cropped_range():
            raise DenoiseCroppedError('You need to set the cropped range first.')
        self.clear()
        
        if contour_method is True and regular_method is True:
            raise RmsMethodError('Only one method is allowed to be True at once.')

        if contour_method is False and regular_method is False:
            raise RmsMethodError('At least one method needs to be used.')
        
        for file in self.spe_file_list:
            curr_all_frame = self.get_current_all_frame(file)
            frame_counter = -1
            plt.close('all')
            for i in range(len(curr_all_frame)):
                frame_counter += 1
                current_frame = curr_all_frame[i]
                zoomed_in_current = current_frame[self.y_start:self.y_end, self.x_start:self.x_end]

                if regular_method:
                    x_intensity_ave_no_bg, y_intensity_ave_no_bg = self.get_intensity_ave_no_bg(curr_all_frame[i])
                    # apply the gaussian filter to de-noise the cropped images.
                    zoomed_in_current = ndi.median_filter(zoomed_in_current, 3)
                if contour_method:
                    x_intensity_ave_no_bg, y_intensity_ave_no_bg, zoomed_in_current = self.get_intensity_ave_using_contour_bg(
                        curr_all_frame[i])

                self.generate_image_and_save(zoomed_in_current, x_intensity_ave_no_bg, y_intensity_ave_no_bg, file,
                                             'single', frame_counter)

    def draw_beam_contour(self, contour_level, lower_diameter_boundary):
        """
        For drawing and imaging plotting purpose. Contours have less fine mesh. The beam/noise with too small contours
        are removed.
        """
        for file in self.spe_file_list:
            curr_all_frame = self.get_current_all_frame(file)
            contr = 0
            for i in range(len(curr_all_frame)):
                contr += 1
                current_frame = curr_all_frame[i]
                zoomed_in_current = current_frame[self.y_start:self.y_end, self.x_start:self.x_end]

                plt.imshow(zoomed_in_current, cmap='jet')
                plt.title(file + "\nFrame#" + str(contr))
                plt.xlabel("x pixel")
                plt.ylabel("y pixel")

                # plot the less fine contour for good looking.
                smooth_results_to_plt = scipy.ndimage.zoom(zoomed_in_current, 0.5)
                x_to_plt = np.linspace(0, len(zoomed_in_current[0]), round(len(zoomed_in_current[0]) * 0.5))
                y_to_plt = np.linspace(0, len(zoomed_in_current), round(len(zoomed_in_current) * 0.5))
                X_to_plt, Y_to_plt = np.meshgrid(x_to_plt, y_to_plt)
                contour_to_plt = plt.contour(X_to_plt, Y_to_plt, smooth_results_to_plt, levels=[contour_level],
                                             colors='r')

                # remove small contours (random noisy spots) for good looking.
                for level in contour_to_plt.collections:
                    for i, path in reversed(list(enumerate(level.get_paths()))):
                        verts = path.vertices
                        diameter = np.max(verts.max(axis=0) - verts.min(axis=0))
                        if diameter < lower_diameter_boundary:
                            del (level.get_paths()[i])
                plt.show()

    def get_main_beam_contour_and_force_outer_zero(self, current_frame, contour_level):
        zoomed_in_current = current_frame[self.y_start:self.y_end, self.x_start:self.x_end]

        smooth_results_to_plt = scipy.ndimage.zoom(zoomed_in_current, 1)
        x = np.linspace(0, len(zoomed_in_current[0]), round(len(zoomed_in_current[0]) * 1))
        y = np.linspace(0, len(zoomed_in_current), round(len(zoomed_in_current) * 1))
        X, Y = np.meshgrid(x, y)
        contour_to_plt = plt.contour(X, Y, smooth_results_to_plt,
                                     levels=[contour_level],
                                     colors='r')
        plt.close('all')

        contour_collection = contour_to_plt.collections[0].get_paths()
        contour_enu = list(contour_collection)
        largest_path = max(contour_enu, key=len).vertices.tolist()

        sorted_by_y = sorted(largest_path, key=lambda tup: tup[1])
        round_sorted_y = [(round(x), round(y)) for x, y in sorted_by_y]

        y_set = {}
        for pointx, pointy in round_sorted_y:
            y_set[pointy] = (y_set[pointy] + [pointx]) if y_set.get(pointy) else [pointx]

        y_cor_list = []
        for y, x_list in y_set.items():
            if y not in y_cor_list:
                y_cor_list.append(y)
            for check_x in it.chain(range(0, min(x_list)), range(max(x_list) + 1, len(zoomed_in_current[0]))):
                zoomed_in_current[y][check_x] = 0

        for enu_x in range(len(zoomed_in_current[0])):
            for row_num in it.chain(range(0, min(y_cor_list)), range(max(y_cor_list), len(zoomed_in_current))):
                zoomed_in_current[row_num][enu_x] = 0

        self.zoom_in_single_frame = zoomed_in_current
        return self.zoom_in_single_frame
