
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Program for determining Young's modulus and yield strength according to DIN EN ISO
# 6892-1 and for the adaptation of material data to the consolidation laws
# by Ludwik, Voce, Ghosh and Hockett-Sherby according to the Least Square approach
# ==============================================================================

# !!!Fill first!!!
# !!!Folder name should not start with U, r, v!!!
# to comment: info of test speed in order

# uncomment to use file directory
# 1. SH 4x2 -- 3x10 -- 6x1 -- 3x1/10 --
#file_name = 'path\Excel_data.xlsx'  # path + file name


sheet = "Werte Serie"
_specimen_ID = np.str(np.asarray(pd.read_excel(io=file_name, sheet_name='Parameter Serie', hearder=None, usecols='C', dtype=str, skiprows=[0, 1, 2], nrows=1))[0][0])
_comparison_specimen_index = 2  # to be compare with
_total_specimen = 3
_total_specimen_blank = 1
# Variable and help factors
_optical_yield_strength = 80.0
_yield_strength = 0.002

_width = np.asarray(pd.read_excel(io=file_name, sheet_name='Statistik', header=None, usecols='M', dtype=float, skiprows=[0, 1], nrows=1))[0]  # total width used
_thickness = np.asarray(pd.read_excel(io=file_name, sheet_name='Statistik', header=None, usecols='L', dtype=float, skiprows=[0, 1], nrows=1))[0]

_width_blank = 3.06  # from micrograph
_thickness_blank = 2.19
_edge_radius_blank = 0.74  # from micrograph
_thickness_iso = 1 / 4 * ((_width - _width_blank) + (_thickness - _thickness_blank))

# ==============================================================================

# Data import
rows = 99995  # if error on data import, re-enter value for row. number should exceed existed rows or exactly the same
rows_skip = [0, 1, 2]
matRaw_data = np.ones((_total_specimen, rows, 2))

# Test data - to import only the x1 speed tested specimens
for i in range(0, _total_specimen):
    start = int((i+9)*3)  # change according to specimen test speeds
    stop = int((i+9)*3+1)
    matRaw_data[i] = np.asarray(pd.read_excel(io=file_name, sheet_name=sheet, header=None, usecols=(start, stop), skiprows=rows_skip, nrows=rows))  # being read: Strain [%] and Standard Force [N]

# Comparison data
# comparison_name = 'path\comparison_datasheet.xlsx'
kupferSpringer = np.array(pd.read_excel(io=comparison_name, sheet_name='Tabelle1', header=None, usecols='E:F', skiprows=rows_skip, nrows=105))

# ==============================================================================

# Strain and Tension
# Technical sizes
_nominal_surface = _width * _thickness
_surface_total = _width * _thickness - (4 * (_edge_radius_blank + _thickness_iso) ** 2 - math.pi * (_edge_radius_blank + _thickness_iso) ** 2)
_surface_coeff = _surface_total / _nominal_surface
_surface_blank = _width_blank * _thickness_blank - (4 * _edge_radius_blank ** 2 - math.pi * _edge_radius_blank ** 2)
_surface_iso = _surface_total - _surface_blank
_total_specimen_iso = _total_specimen - _total_specimen_blank
technicalStrain = np.array(matRaw_data[:, :, 0] / 100)
technicalTension = np.array(matRaw_data[:, :, 1] / _nominal_surface)
pullingForce = matRaw_data[:, :, 1]

# True values under the assumption of volume consistency
trueStrain = np.log(technicalStrain + 1)
trueTension = technicalTension * (technicalStrain + 1)
trueStrain = np.nan_to_num(trueStrain)
trueTension = np.nan_to_num(trueTension)
# ==============================================================================

# Surface and pulling force calculations
# Representation of the tensile force curves on the same data points of the technical strain
# Help factors
phi_max = max(technicalStrain[0, 0:_total_specimen])
range_phi = np.arange(0, phi_max, 0.0001)
pullingForce_int = np.zeros((_total_specimen, len(range_phi)))

for i in range(0, _total_specimen):
    XX_help, index_help = np.unique(technicalStrain[i, 0:_total_specimen, ], return_index=True)
    pullingForce_int[i, :] = np.interp(range_phi, XX_help, pullingForce[i, index_help])

sum_pullingForce_mit_iso = np.ones((1, len(pullingForce_int[1, :])))
sum_pullingForce_blank = np.ones((1, len(pullingForce_int[1, :])))

for v in range(0, _total_specimen - _total_specimen_blank):
    sum_pullingForce_mit_iso = sum_pullingForce_mit_iso + pullingForce_int[v, :]

for y in range(_total_specimen_iso, _total_specimen):
    sum_pullingForce_blank = sum_pullingForce_blank + pullingForce_int[y, :]

pullingForce_blank = 1 / _total_specimen_blank * sum_pullingForce_blank
pullingForce_with_iso = 1 / _total_specimen_iso * sum_pullingForce_mit_iso
pullingForce_iso = pullingForce_with_iso - pullingForce_blank  # Force from isolation layer
# ==============================================================================

# Tension calculation with and without Isolation
tension_blank = pullingForce_blank / _surface_blank
tension_iso = pullingForce_iso / _surface_iso
# to get the real technical tension with the right nominal area
technicalTension = technicalTension / _surface_coeff
# ==============================================================================

# E-Modul nach DIN EN ISO 6892-1 Anhang G
# Help factors
R_under_10 = (_optical_yield_strength - technicalTension[:, 0]) * 0.1 + technicalTension[:, 0]
R_over_40 = (_optical_yield_strength - technicalTension[:, 0]) * 0.4 + technicalTension[:, 0]
n_10 = np.ones((len(technicalStrain[:, 1]), 1))
n_40 = np.ones((len(technicalStrain[:, 1]), 1))
eModul_din = np.zeros((len(technicalStrain[:, 1]), 1))

# Ausgleichsgerade
for m in range(0, len(technicalStrain[:, 1])):
    n_10[m] = np.argmax(technicalTension[m, :] > R_under_10[m]) - 1
    n_40[m] = np.argmax(technicalTension[m, :] > R_over_40[m]) - 1
    help10 = int(n_10[m])
    help40 = int(n_40[m])

    magVector_x = technicalTension[m, help10:help40] - technicalTension[m][0]
    obsVector_x = technicalStrain[m, help10:help40]
    interimValue = 0

    for n in range(0, help40 - help10):
        interimValue = interimValue + obsVector_x[n] * obsVector_x[n]

    interimValue_1 = (1 / interimValue) * obsVector_x
    eModul_din[m] = 0

    for o in range(0, help40 - help10):
        eModul_din[m] = eModul_din[m] + interimValue_1[o] * magVector_x[o]

# E_modul for Isolation
# Help factors
#optical_yield_strength_iso = 50
#R_under_10_iso = (optical_yield_strength_iso - tension_iso[0]) * 0.1 + tension_iso[0]
#R_over_40_iso = (optical_yield_strength_iso - tension_iso[0]) * 0.4 + tension_iso[0]

# Equalizing line
#n_10_iso = np.argmax(tension_iso[:] > R_under_10_iso) - 1
#n_40_iso = np.argmax(tension_iso[:] > R_over_40_iso) - 1
#e_modul_iso = (tension_iso[n_40_iso] - tension_iso[n_10_iso]) / (range_phi[n_40_iso].reshape((-1, 1)) - range_phi[n_10_iso].reshape((-1, 1)))

# E-Modul for Blank wire( pure Copper)
# Help factors
#optical_yield_strength_blank = _optical_yield_strength
#R_under_10_blank = (_optical_yield_strength - tension_blank[0]) * 0.1 + tension_blank[0]
#R_over_40_blank = (_optical_yield_strength - tension_blank[0]) * 0.4 + tension_blank[0]

#n_10_blank = np.argmax(tension_blank > R_under_10_blank) - 1
#n_40_blank = np.argmax(tension_blank > R_over_40_blank) - 1
#e_modul_blank = (tension_blank[n_40_blank] - tension_blank[n_10_blank]) / (range_phi[n_40_blank].reshape((-1, 1)) - range_phi[n_10_blank].reshape((-1, 1)))

# ==============================================================================

# Dehngrenze nach DIN EN ISO 6892-1 Anhang G
# Help factors
elastic_measuring_range_DIN = np.ones((len(technicalStrain[:, 1]), 1))
yield_strength_DIN = np.zeros((len(trueStrain[:, 1]), 1))
elastic_range_DIN = np.zeros((len(technicalStrain[:, 1]), 1))

# Calculation
for p in range(0, len(technicalStrain[:, 1])):
    elastic_measuring_range_DIN[p] = np.argmax((technicalStrain[p, :] - technicalTension[p, :] / eModul_din[p]) > _yield_strength)  # R_p_0,2% indices
    yield_strength_DIN[p] = technicalTension[p, int(elastic_measuring_range_DIN[p])]  # assign R_p_0,2
    elastic_range_DIN[p] = technicalStrain[p, int(elastic_measuring_range_DIN[p])]
# ==============================================================================

# Flow curve extraction
temp = len(trueStrain[1, :]) - max(elastic_measuring_range_DIN)
flowCurve = np.zeros((2, int(len(trueStrain[:, 1])), int(temp)))

for q in range(0, len(trueStrain[:, 1])):
    temp = len(trueStrain[1, :]) - (max(elastic_measuring_range_DIN) - elastic_measuring_range_DIN[q])
    flowCurve[0, q, :] = trueStrain[q, int(elastic_measuring_range_DIN[q]):int(temp)] - trueStrain[
        0, int(elastic_measuring_range_DIN[q])]
    flowCurve[1, q, :] = trueTension[q, int(elastic_measuring_range_DIN[q]):int(temp)]
# ==============================================================================

# Flow curve fitting
# Fitting only over relevant measured value range
relevant_measuring_range = np.ones((len(technicalStrain[:, 1]), 1))
minimum_allowable_tension_difference = -0.4
technical_tension_help = np.ones((int(len(technicalTension[:, 1])), int(len(technicalTension[1, :]))))
technical_tension_help[:, 0: (int(len(technicalTension[1, :])) - 1)] = technicalTension[:, 1:int(len(technicalTension[1, :]))]
for s in range(0, len(technicalStrain[:, 1])):
    relevant_measuring_range[s] = np.argmax((technical_tension_help[s, :] - technicalTension[s, :]) < minimum_allowable_tension_difference) - 1

# Database
databaseStrain = flowCurve[0, _comparison_specimen_index-1, 1:int(relevant_measuring_range[_comparison_specimen_index-1] - elastic_measuring_range_DIN[_comparison_specimen_index-1])]
databaseTension = flowCurve[1, _comparison_specimen_index-1, 1:int(relevant_measuring_range[_comparison_specimen_index-1] - elastic_measuring_range_DIN[_comparison_specimen_index-1])]

# clear negative values, if there is any
for i in range(0, 150):
    if databaseStrain[i] < 0:
        last_neg_datStrain = np.where(databaseStrain <= 0)[0][-1]+1
        databaseStrain_alt = databaseStrain[int(last_neg_datStrain):int(len(databaseStrain))]
        databaseTension_alt = databaseTension[int(last_neg_datStrain): int(len(databaseStrain))]
        databaseStrain = databaseStrain_alt
        databaseTension = databaseTension_alt


def fun_ludwik(data_ludwik, *y_ludwik):
    return y_ludwik[0] + y_ludwik[1] * data_ludwik ** y_ludwik[2]


y_ludwik_0 = [115, 300, 0.8]
lower_bound_3 = [0, 0, 0]
[y_ludwik_res, res_ludwik] = curve_fit(fun_ludwik, databaseStrain, databaseTension,p0=y_ludwik_0, check_finite=False, bounds=(lower_bound_3, np.inf), method='trf')


def fun_voce(data_voce=np.asarray, *y_voce):
    return y_voce[0] + y_voce[1] * (1 - np.exp(-data_voce / y_voce[2]) + y_voce[3] * data_voce)


y_voce_0 = [115, 400, 1, 10]
lower_bound_4 = [0, 0, 0, 0]
[y_voce_res, res_voce] = curve_fit(fun_voce, databaseStrain, databaseTension, p0=y_voce_0, check_finite=False, bounds=(lower_bound_4, np.inf), method='trf')


def fun_ghosh(data_ghosh=np.asarray, *y_ghosh):
    return y_ghosh[0] * (y_ghosh[1] + data_ghosh) ** y_ghosh[2] + y_ghosh[3]


y_ghosh_0 = [300, 1, 0.5, 115]
[y_ghosh_res, res_ghosh] = curve_fit(fun_ghosh, databaseStrain, databaseTension, p0=y_ghosh_0, check_finite=False, bounds=(lower_bound_4, np.inf), method='trf')


def fun_hockett_sherby(data_hockett_sherby=np.asarray, *y_hockett_sherby):
    return y_hockett_sherby[0] - (y_hockett_sherby[0] - y_hockett_sherby[1]) * np.exp(
        -y_hockett_sherby[2] * data_hockett_sherby ** y_hockett_sherby[3])


y_hockettSherby_0 = [400, 115, 10, 0.5]
lower_bound_3 = [0, 0, 0]
[y_hockettSherby_res, res_hockettSherby] = curve_fit(fun_hockett_sherby, databaseStrain, databaseTension, p0=y_hockettSherby_0, check_finite=False, bounds=(lower_bound_4, np.inf), method='trf', maxfev=5000)

# ==============================================================================

# Ergebnisvisualisierung
# Datenberechnung
daten_range = np.asarray(np.arange(0, 1, 0.01))
data_ludwik = fun_ludwik(daten_range, *y_ludwik_res)
data_voce = fun_voce(daten_range, *y_voce_res)
data_ghosh = fun_ghosh(daten_range, *y_ghosh_res)
data_hockettSherby = fun_hockett_sherby(daten_range, *y_hockettSherby_res)
flowCurve[flowCurve == 0] = np.nan

# Datenausgabe
# Plot
# first figure: flow curve comparison
fig, ax = plt.subplots(2)
fig.set_size_inches(18.5, 10.5)
springer, = ax[0].plot(kupferSpringer[:, 0], kupferSpringer[:, 1], 'b', label='Springer')
ludwik, = ax[0].plot(daten_range, data_ludwik, 'm', label='Ludwik')
voce, = ax[0].plot(daten_range, data_voce, 'r', label='Voce')
ghosh, = ax[0].plot(daten_range, data_ghosh, 'g', label='Ghosh')
hS, = ax[0].plot(daten_range, data_hockettSherby, 'k', label='Hockett-Sherby')
wbk, = ax[0].plot(flowCurve[0, _comparison_specimen_index-1, :], flowCurve[1, _comparison_specimen_index-1, :], marker='*', markersize=5, markevery=500, label=_specimen_ID)
ax[0].set(xlabel='Umformgrad \u03C6 [-]', ylabel='wahre Spannung $k_f$ [MPa]', title='Fliesskurven ' + ' #' + np.str(_comparison_specimen_index))
ax[0].grid()
leg0 = ax[0].legend(loc='lower right', shadow=None, edgecolor='black')
ax[0].set_xlim([0, 0.4])
ax[0].set_ylim([0, 400])

# create pickable legend
lines = [springer, ludwik, voce, ghosh, hS, wbk]
lined = dict()

for legline, origline in zip(leg0.get_lines(), lines):
    legline.set_picker(5)   # tolerance
    lined[legline] = origline


def onpick(event):
    legline = event.artist
    origline = lined[legline]
    vis = not origline.get_visible()
    origline.set_visible(vis)

    if vis:
        legline.set_alpha(1.0)
    else:
        legline.set_alpha(0.2)
    fig.canvas.draw()


fig.canvas.mpl_connect('pick_event', onpick)
# end pickable legend

# second figure: pulling Force
for i in range(0, _total_specimen):
    ax[1].plot(technicalStrain[i, 0:int(relevant_measuring_range[i])], pullingForce[i, 0:int(relevant_measuring_range[i])], label='P' + np.str(i + 1))
ax[1].set(xlabel='technische Dehnung \u03B5 [%]', ylabel='Zugkraft [N]', title = 'Zugkraefte ')
ax[1].text(0.2, 500,_specimen_ID + '\nKenngroessen nach DIN 6892-1, Anhang G\n=========================\nBreite = ' + np.str(_width) + 'mm\nDicke = ' + np.str(_thickness) +'mm\nE-Modul = ' + np.str(np.around(np.mean(eModul_din)/1000,2)) +' GPa\nDehngrenze = ' + np.str(np.around(np.mean(yield_strength_DIN), 2)) + ' MPa\nelastischer Bereich: 0.00 % .. ' + np.str(np.around(np.mean(elastic_range_DIN) * 100, 2)) + '%\nE-Modul-Ermittlung: ' + np.str(np.around(np.mean(elastic_range_DIN) * 10, 2)) + '% .. ' + np.str(np.around(np.mean(elastic_range_DIN) * 40, 2)) + '% ', bbox=dict(boxstyle='square', fc='none', ec='black'))
ax[1].legend(loc='lower right', shadow=None, edgecolor='blue')
ax[1].set_xlim([0, 0.6])
ax[1].grid()
plt.subplots_adjust(hspace=0.3)
plt.show()


