# from gen_matrix import *
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
# from func import *
from numpy.linalg import multi_dot
import random

from spe_image_rms_rms_error import *
import settings


z0 = 0                                  # start of the transport matrix
ze = 254.7 - 0.782 * 100 + 0.5 * 100    # end of the transport matrix
zstep = 0.01                            # step size of the sliced magnetic field
ek = 0.7637                             # kinetic energy of the beam (MeV)
e0 = 0.511                              # electron rest energy (MeV)
e = ek + e0

# constants
c = 299792458
me = 9.10938291e-31
qe = 1.602176565e-19
freq = 1.3161276885390001e9
# rest mass energy (eV)
rest_e = me * c**2 / qe
betagamma = np.sqrt(e**2 / e0**2 - 1)
# ======================================================================================================================
SPE_image = DenoiseSPEImage(settings.FILE_START_STRING, settings.FOLDER_PATH, settings.FILE_LIST)
SPE_image.set_cropped_range(100, 480, 200, 500)
# SPE_image.set_background_range(0, 8, 0, 8)
# SPE_image.draw_beam_contour(100, 50)
x_rms, y_rms, x_std, y_std = SPE_image.get_rms_and_rms_error(contour_method=False, regular_method=True)
# SPE_image.plot_single_frame(contour_method=False, regular_method=True)

# plot_all_ave_frames(file_ls, path, data_set, 100, 480, 200, 500)
print("x RMS size: ", x_rms)
print("x RMS STD: ", x_std)
print("y RMS size: ", y_rms)
print("y RMS STD: ", y_std)
quit()

# blue solenoid count list
start_num = 3
end_num = 3
count_ls = np.linspace(832+32*start_num, 1152-32*end_num, int((1152 - 32*end_num - 832 - 32*start_num)/32 + 1), endpoint=True)
print("Solenoid field (counts) scan list: ", count_ls)
# blue solenoid current list
current_sol = [convert_cont_to_current(int(i)) for i in count_ls]
# blue solenoid B field [T] list
i_sol = [i/10000 for i in get_bsol_field(current_sol)]

R = []
for i in range(len(i_sol)):
    Is = abs(i_sol[i])
    RR = trans_matrix(ek, Is, z0, ze, zstep)
    R.append(RR)

n = len(i_sol)
A = np.zeros((len(i_sol), 3), float).tolist()
# print(A)
counter = -1
for i in range(n):
    counter += 1
    A[counter][0] = R[i][0][0] ** 2
    A[counter][1] = R[i][0][1] ** 2
    A[counter][2] = - 2 * R[i][0][0] * R[i][0][1]
    
# ======================================================================================================================
sigma_11_x = [i**2 for i in x_rms[start_num: len(x_rms)-end_num]]
var_sig11_x = np.array([1 / i for i in x_std[start_num: len(x_rms)-end_num]])     # x weight matrix of the RMS size standard deviation
Awx = np.array(A) * var_sig11_x[:, np.newaxis]
sig_11_xw = np.array(sigma_11_x) * var_sig11_x
twiss_Px, residuals_x, rank_x, cc_x = np.linalg.lstsq(np.array(Awx), np.array(sig_11_xw), rcond=None)
emit_x_un_sq = twiss_Px[0] * twiss_Px[1] - twiss_Px[2] ** 2

x_rms_temp = []
emit_x_un_sq_ls = []
contr = 0
while emit_x_un_sq <= 0:
    contr += 1
    print("\rtrying to modify x RMS list, No.", contr, end="")
    std_x_add = [np.random.uniform(-x_std[i], x_std[i]) for i in range(len(x_std))]
    x_rms_temp = list(np.array(x_rms) + np.array(std_x_add))

    sigma_11_x_temp = [i ** 2 for i in x_rms_temp[start_num: len(x_rms) - end_num]]
    sig_11_xw_temp = np.array(sigma_11_x_temp) * var_sig11_x
    twiss_Px, residuals_x, rank_x, cc_x = np.linalg.lstsq(np.array(Awx), np.array(sig_11_xw_temp), rcond=None)
    emit_x_un_sq = twiss_Px[0] * twiss_Px[1] - twiss_Px[2] ** 2

if x_rms_temp:
    print('\n')
    print("midified x RMS: ", x_rms_temp)
else:
    print('no change on x RMS')

sigma_11_y = [i ** 2 for i in y_rms[start_num: len(y_rms)-end_num]]
var_sig11_y = np.array([1 / i for i in y_std[start_num: len(y_rms)-end_num]])
Awy = np.array(A) * var_sig11_y[:, np.newaxis]
sig_11_yw = np.array(sigma_11_y) * var_sig11_y
twiss_Py, residuals_y, rank_y, cc_y = np.linalg.lstsq(np.array(Awy), np.array(sig_11_yw), rcond=None)
emit_y_un_sq = twiss_Py[0] * twiss_Py[1] - twiss_Py[2] ** 2

y_rms_temp = []
emit_y_un_sq_ls = []
contr = 0
while emit_y_un_sq <= 0:
    contr += 1
    print("\rtrying to modify y RMS list, No.", contr, end="")
    std_y_add = [np.random.uniform(-y_std[i], y_std[i]) for i in range(len(y_std))]
    y_rms_temp = list(np.array(y_rms) + np.array(std_y_add))

    sigma_11_y_temp = [i ** 2 for i in y_rms_temp[start_num: len(y_rms) - end_num]]
    sig_11_yw_temp = np.array(sigma_11_y_temp) * var_sig11_y
    twiss_Py, residuals_y, rank_y, cc_y = np.linalg.lstsq(np.array(Awy), np.array(sig_11_yw_temp), rcond=None)
    emit_y_un_sq = twiss_Py[0] * twiss_Py[1] - twiss_Py[2] ** 2

if y_rms_temp:
    print('\n')
    print("modified y RMS: ", y_rms_temp)
else:
    print('no change on y RMS')

emitx = np.sqrt(emit_x_un_sq) * betagamma
emity = np.sqrt(emit_y_un_sq) * betagamma


# compute the confidence intervals
n_x = len(np.array(sig_11_xw))
k_x = len(twiss_Px)

sigmax = np.sum((np.array(sig_11_xw) - np.dot(np.array(Awx), twiss_Px))**2) / (n_x - k_x)
sigmay = np.sum((np.array(sig_11_yw) - np.dot(np.array(Awy), twiss_Py))**2) / (n_x - k_x)

covar_x = np.linalg.inv(np.dot(np.array(Awx).T, np.array(Awx))) * sigmax   # covariance matrix
covar_y = np.linalg.inv(np.dot(np.array(Awy).T, np.array(Awy))) * sigmay   # covariance matrix
# print(covar_x)
# se = np.sqrt(np.diag(covar_x))         # standard error

# # 写出emt_un对x1,x2,x3的偏导数
pd_1_x = twiss_Px[1] / 2 / (emitx / betagamma)
pd_2_x = twiss_Px[0] / 2 / (emitx / betagamma)
pd_3_x = -twiss_Px[2] / (emitx / betagamma)

pd_1_y = twiss_Py[1] / 2 / (emity / betagamma)
pd_2_y = twiss_Py[0] / 2 / (emity / betagamma)
pd_3_y = -twiss_Py[2] / (emity / betagamma)


print("--------------------------")
print("Calculated x twiss parameters: ", twiss_Px)
print("--------------------------")
print("Calculated y twiss parameters: ", twiss_Py)


emt_un_var_x = multi_dot([[pd_1_x, pd_2_x, pd_3_x], covar_x, np.transpose([pd_1_x, pd_2_x, pd_3_x])])
emt_un_error_x = np.sqrt(emt_un_var_x)
emt_err_x = emt_un_error_x * betagamma

emt_un_var_y = multi_dot([[pd_1_y, pd_2_y, pd_3_y], covar_y, np.transpose([pd_1_y, pd_2_y, pd_3_y])])
emt_un_error_y = np.sqrt(emt_un_var_y)
emt_err_y = emt_un_error_y * betagamma

# =====================================================================================================================
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)

x_new = np.linspace(min(i_sol), max(i_sol), 500)

if x_rms_temp:
    print('x RMS is modified')
    popt, pcov = curve_fit(quadratic, i_sol, sigma_11_x_temp, maxfev=100000000)
    x = np.linspace(min(i_sol), max(i_sol), 1000, endpoint=True)
    ax1.plot(x, quadratic(x, *popt), label=r'fit of $\sigma{_{x}}^2$', c='b', alpha=0.5)

    popt2, pcov2 = curve_fit(quadratic, i_sol, sigma_11_y, maxfev=100000000)
    y = np.linspace(min(i_sol), max(i_sol), 1000, endpoint=True)
    ax1.plot(y, quadratic(y, *popt2), c='r', label=r'fit of $\sigma{_{y}}^2$', alpha=0.5)

    ax1.scatter(i_sol, sigma_11_x_temp, c='b', marker='o', s=30, label=r'$\sigma{_{x}}^2$', alpha=0.9, zorder=10)
    ax1.scatter(i_sol, sigma_11_y, c='r', marker='o', s=30, label=r'$\sigma{_{y}}^2$', alpha=0.9, zorder=10)

elif y_rms_temp:
    print("y RMS is modified.")
    popt, pcov = curve_fit(quadratic, i_sol, sigma_11_x, maxfev=100000000)
    x = np.linspace(min(i_sol), max(i_sol), 1000, endpoint=True)
    ax1.plot(x, quadratic(x, *popt), label=r'fit of $\sigma{_{x}}^2$', c='b', alpha=0.5)

    popt2, pcov2 = curve_fit(quadratic, i_sol, sigma_11_y_temp, maxfev=100000000)
    y = np.linspace(min(i_sol), max(i_sol), 1000, endpoint=True)
    ax1.plot(y, quadratic(y, *popt2), c='r', label=r'fit of $\sigma{_{y}}^2$', alpha=0.5)

    ax1.scatter(i_sol, sigma_11_x, c='b', marker='o', s=30, label=r'$\sigma{_{x}}^2$', alpha=0.9, zorder=10)
    ax1.scatter(i_sol, sigma_11_y_temp, c='r', marker='o', s=30, label=r'$\sigma{_{y}}^2$', alpha=0.9, zorder=10)

elif y_rms_temp and x_rms_temp:
    print('x and y RMS are modified.')
    popt, pcov = curve_fit(quadratic, i_sol, sigma_11_x_temp, maxfev=100000000)
    x = np.linspace(min(i_sol), max(i_sol), 1000, endpoint=True)
    ax1.plot(x, quadratic(x, *popt), label=r'fit of $\sigma{_{x}}^2$', c='b', alpha=0.5)

    popt2, pcov2 = curve_fit(quadratic, i_sol, sigma_11_y_temp, maxfev=100000000)
    y = np.linspace(min(i_sol), max(i_sol), 1000, endpoint=True)
    ax1.plot(y, quadratic(y, *popt2), c='r', label=r'fit of $\sigma{_{y}}^2$', alpha=0.5)

    ax1.scatter(i_sol, sigma_11_x_temp, c='b', marker='o', s=30, label=r'$\sigma{_{x}}^2$', alpha=0.9, zorder=10)
    ax1.scatter(i_sol, sigma_11_y_temp, c='r', marker='o', s=30, label=r'$\sigma{_{y}}^2$', alpha=0.9, zorder=10)

else:
    print('Nothing is modified.')
    popt, pcov = curve_fit(quadratic, i_sol, sigma_11_x, maxfev=100000000)
    x = np.linspace(min(i_sol), max(i_sol), 1000, endpoint=True)
    ax1.plot(x, quadratic(x, *popt), label=r'fit of $\sigma{_{x}}^2$', c='b', alpha=0.5)
    
    popt2, pcov2 = curve_fit(quadratic, i_sol, sigma_11_y, maxfev=100000000)
    y = np.linspace(min(i_sol), max(i_sol), 1000, endpoint=True)
    ax1.plot(y, quadratic(y, *popt2), c='r', label=r'fit of $\sigma{_{y}}^2$',alpha=0.5)
    
    ax1.scatter(i_sol, sigma_11_x, c='b', marker='o', s=30, label=r'$\sigma{_{x}}^2$', alpha=0.9, zorder=10)
    ax1.scatter(i_sol, sigma_11_y, c='r', marker='o', s=30, label=r'$\sigma{_{y}}^2$', alpha=0.9, zorder=10)


plt.grid(alpha=0.5, ls=':')
# plt.title(f'Calculated emittance of {data_set} - 01/09/2020', fontsize=14)
plt.xlabel(r'Solenoid field (T)', fontsize=13)
plt.ylabel(r'$\sigma{_{x}}^2$, $\sigma{_{y}}^2$ (mm$^2$)', fontsize=13)

plt.legend(fancybox=False, framealpha=0.7, edgecolor='k', loc='best')
plt.xlim(min(i_sol), max(i_sol))

# plt.annotate(r'$\varepsilon{_x}$ = %.3f $\pm$ %.3f mm mrad <-> %.3f mm mrad/mm'
#              % (emitx, emt_err_x, emitx / 0.168),
#              ((max(i_sol)-min(i_sol)) / 2 + min(i_sol),
#               (max(sigma_11_x)-min(sigma_11_x))/1.3 + min(sigma_11_x)), horizontalalignment='center')
#
# plt.annotate(r'$\varepsilon{_y}$ = %.3f $\pm$ %.3f mm mrad <-> %.3f mm mrad/mm'
#              % (emity, emt_err_y, emity / 0.132),
#              ((max(i_sol)-min(i_sol)) / 2 + min(i_sol),
#               (max(sigma_11_x)-min(sigma_11_x))/1.6 + min(sigma_11_x)), horizontalalignment='center')
#
# plt.annotate(r'UIC MTE: 266 meV <-> $\varepsilon$ = 0.721 mrad',
#              ((max(i_sol)-min(i_sol)) / 2 + min(i_sol),
#               (max(sigma_11_x)-min(sigma_11_x))/2.0 + min(sigma_11_x)), horizontalalignment='center')


plt.annotate(r'$\varepsilon{_x}$ = %.3f $\pm$ %.3f mm mrad'
             % (emitx, emt_err_x),
             ((max(i_sol)-min(i_sol)) / 2 + min(i_sol),
              (max(sigma_11_x)-min(sigma_11_x))/1.3 + min(sigma_11_x)), horizontalalignment='center')

plt.annotate(r'$\varepsilon{_y}$ = %.3f $\pm$ %.3f mm mrad'
             % (emity, emt_err_y),
             ((max(i_sol)-min(i_sol)) / 2 + min(i_sol),
              (max(sigma_11_x)-min(sigma_11_x))/1.6 + min(sigma_11_x)), horizontalalignment='center')

plt.show()

# #
# plot_ls = ['_03', '_04', '_05', '_06', '_07', '_08']
# # plot_ls = ['_01']
# for i in range(len(plot_ls)):
#     print(f'******** Saving the single frames for set{plot_ls[i]} ********')
#     plot_single_frame(file_ls, path, data_set + plot_ls[i], 100, 480, 150, 500)
#     plt.close('all')
#
# print('Finished.')

plt.scatter(i_sol, x_rms[start_num: len(x_rms)-end_num], c='b', marker='o', s=30, label=r'$\sigma{_{x}}$ with error bars', alpha=0.9, zorder=10)
plt.errorbar(i_sol, x_rms[start_num: len(x_rms)-end_num], yerr=x_std[start_num: len(x_rms)-end_num], fmt='.b', capsize=3.5,
             color='black', elinewidth=0.7)
plt.plot(i_sol, x_rms[start_num: len(x_rms)-end_num], c='b', alpha=0.5)

plt.scatter(i_sol, y_rms[start_num: len(y_rms)-end_num], c='r', marker='o', s=30, label=r'$\sigma{_{y}}$ with error bars', alpha=0.9, zorder=10)
plt.errorbar(i_sol, y_rms[start_num: len(y_rms)-end_num], yerr=y_std[start_num: len(x_rms)-end_num], fmt='.r', capsize=3.5,
             color='black', elinewidth=0.7)
plt.plot(i_sol, y_rms[start_num: len(y_rms)-end_num], c='r', alpha=0.5)

plt.annotate(r'$\varepsilon{_x}$ = %.3f $\pm$ %.3f mm mrad'
             % (emitx, emt_err_x),
             ((max(i_sol)-min(i_sol)) / 2 + min(i_sol),
              (max(x_rms)-min(x_rms))/1.8 + min(x_rms)), horizontalalignment='center')

plt.annotate(r'$\varepsilon{_y}$ = %.3f $\pm$ %.3f mm mrad'
             % (emity, emt_err_y),
             ((max(i_sol)-min(i_sol)) / 2 + min(i_sol),
              (max(x_rms)-min(x_rms))/2 + min(x_rms)), horizontalalignment='center')
plt.grid(alpha=0.5, ls=':')
plt.legend(fancybox=False, framealpha=0.7, edgecolor='k', loc='best')
plt.xlim(0.027, 0.033)

plt.xlabel(r'Solenoid field (T)', fontsize=13)
plt.ylabel(r'$\sigma{_{x}}$, $\sigma{_{y}}$ (mm)', fontsize=13)

plt.show()