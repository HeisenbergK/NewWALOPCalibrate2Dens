import csv
import sys
import multiprocessing as mp
from read_zemax import read_zmx, read_zmx_2
from os import walk
from astropy.table import Table, unique
from scipy.optimize import curve_fit
import platform
import warnings
import random
import matplotlib.pyplot as plt
from sympy import symbols
from sympy.core import sympify
from sympy.solvers.solveset import nonlinsolve
from itertools import compress
from math import dist
from models import *

warnings.filterwarnings("ignore")


def testerfunc(testfieldcase, numfields=None, verbose=5, order=2, mixed=True, plim=1.0, ret_model=False):
    print(f"Starting field: {testfieldcase}")
    columnames = ["config", "qu", "xy", "trans"]
    newcolumnames = ['quin', 'trans1', 'trans2', 'trans3', 'trans4', 'quout']

    # get the file names from the folders
    filenames1 = next(walk("Config_1"), (None, None, []))[2]
    filenames2 = next(walk("Config_2"), (None, None, []))[2]
    filenames3 = next(walk("Config_3"), (None, None, []))[2]
    filenames4 = next(walk("Config_4"), (None, None, []))[2]
    unpolfilenames1 = next(walk("Unpol/Config_1"), (None, None, []))[2]
    unpolfilenames2 = next(walk("Unpol/Config_2"), (None, None, []))[2]
    unpolfilenames3 = next(walk("Unpol/Config_3"), (None, None, []))[2]
    unpolfilenames4 = next(walk("Unpol/Config_4"), (None, None, []))[2]

    # make a filename list
    if platform.system == "Windows":
        filelist = [["Config_1\\%s" % i for i in filenames1], ["Config_2\\%s" % i for i in filenames2],
                    ["Config_3\\%s" % i for i in filenames3], ["Config_4\\%s" % i for i in filenames4]]
        unpolfilelist = [["Unpol\\Config_1\\%s" % i for i in unpolfilenames1],
                         ["Unpol\\Config_2\\%s" % i for i in unpolfilenames2],
                         ["Unpol\\Config_3\\%s" % i for i in unpolfilenames3],
                         ["Unpol\\Config_4\\%s" % i for i in unpolfilenames4]]
    else:
        filelist = [["Config_1/%s" % i for i in filenames1], ["Config_2/%s" % i for i in filenames2],
                    ["Config_3/%s" % i for i in filenames3], ["Config_4/%s" % i for i in filenames4]]
        unpolfilelist = [["Unpol/Config_1/%s" % i for i in unpolfilenames1],
                         ["Unpol/Config_2/%s" % i for i in unpolfilenames2],
                         ["Unpol/Config_3/%s" % i for i in unpolfilenames3],
                         ["Unpol/Config_4/%s" % i for i in unpolfilenames4]]

    # make them in a table
    rowlist = []
    for i in range(0, 4):
        for j in range(0, len(filelist[i])):
            q, u, fields, tottransmissions, fieldx, fieldy = read_zmx(filelist[i][j], wavelengthno=7)
            for k in range(0, len(fields)):
                rowlist.append([i, [q, u], [fields[k][0], fields[k][1]], tottransmissions[k]])
        for j in range(0, len(unpolfilelist[i])):
            q, u, fields, tottransmissions, fieldx, fieldy = read_zmx_2(unpolfilelist[i][j], wavelengthno=7)
            for k in range(0, len(fields)):
                rowlist.append([i, [q, u], [fields[k][0], fields[k][1]], tottransmissions[k]])

    master = Table(list(map(list, zip(*rowlist))), names=columnames)

    # get the unique
    uniquefields = unique(master, keys="xy")["xy"].tolist()
    uniquepol = unique(master, keys="qu")["qu"].tolist()
    testfield = uniquefields[testfieldcase]

    if verbose >= 3:
        print("-----------------------------")
        print(f"{len(uniquefields)} unique field positions")
        print("-----------------------------")
        print(f"{len(uniquepol)} unique polarizations")
        print("-----------------------------")

    slave = master[np.all(master["xy"] == testfield, axis=1)]
    slave.remove_column('xy')

    rowlist = []
    for i in range(0, len(uniquepol)):
        current_qu = uniquepol[i]
        trans1 = False
        trans2 = False
        trans3 = False
        trans4 = False
        for j in range(0, len(slave)):
            if slave['qu'][j].tolist() == current_qu:
                if slave['config'][j] == 0:
                    trans1 = slave['trans'][j]
                if slave['config'][j] == 1:
                    trans2 = slave['trans'][j]
                if slave['config'][j] == 2:
                    trans3 = slave['trans'][j]
                if slave['config'][j] == 3:
                    trans4 = slave['trans'][j]
        if not (trans1 and trans2 and trans3 and trans4):
            raise Exception(f"Transmittance not found for {current_qu}")
        quout = [(trans3-trans1)/(trans3+trans1), (trans2-trans4)/(trans2+trans4)]
        rowlist.append([current_qu, trans1, trans2, trans3, trans4, quout])

    newmaster = Table(list(map(list, zip(*rowlist))), names=newcolumnames)

    # fit a model
    xdata = newmaster['quout']
    ydata = newmaster['quin']
    ydata_q = ydata[:, 0]
    ydata_u = ydata[:, 1]
    if order == 2:
        if mixed:
            popt_trans1, pcov_trans1 = curve_fit(modell, ydata, newmaster['trans1'])
            popt_trans2, pcov_trans2 = curve_fit(modell, ydata, newmaster['trans2'])
            popt_trans3, pcov_trans3 = curve_fit(modell, ydata, newmaster['trans3'])
            popt_trans4, pcov_trans4 = curve_fit(modell, ydata, newmaster['trans4'])
            if ret_model:
                return [testfield[0], testfield[1],
                        popt_trans1.tolist(), popt_trans2.tolist(), popt_trans3.tolist(), popt_trans4.tolist(),
                        pcov_trans1.tolist(), pcov_trans2.tolist(), pcov_trans3.tolist(), pcov_trans4.tolist()]
        else:
            popt_trans1, pcov_trans1 = curve_fit(modellnm, ydata, newmaster['trans1'])
            popt_trans2, pcov_trans2 = curve_fit(modellnm, ydata, newmaster['trans2'])
            popt_trans3, pcov_trans3 = curve_fit(modellnm, ydata, newmaster['trans3'])
            popt_trans4, pcov_trans4 = curve_fit(modellnm, ydata, newmaster['trans4'])
            if ret_model:
                return [testfield[0], testfield[1],
                        popt_trans1.tolist(), popt_trans2.tolist(), popt_trans3.tolist(), popt_trans4.tolist(),
                        pcov_trans1.tolist(), pcov_trans2.tolist(), pcov_trans3.tolist(), pcov_trans4.tolist()]
    elif order == 3:
        popt_trans1, pcov_trans1 = curve_fit(modell3, ydata, newmaster['trans1'])
        popt_trans2, pcov_trans2 = curve_fit(modell3, ydata, newmaster['trans2'])
        popt_trans3, pcov_trans3 = curve_fit(modell3, ydata, newmaster['trans3'])
        popt_trans4, pcov_trans4 = curve_fit(modell3, ydata, newmaster['trans4'])
        if ret_model:
            return [testfield[0], testfield[1],
                    popt_trans1.tolist(), popt_trans2.tolist(), popt_trans3.tolist(), popt_trans4.tolist(),
                    pcov_trans1.tolist(), pcov_trans2.tolist(), pcov_trans3.tolist(), pcov_trans4.tolist()]
    else:
        raise Exception("Wrong order")

    rowsel = -1
    zero_row = Table()
    for k in range(0, len(newmaster)):
        if newmaster[k]['quin'][0] == 0 and newmaster[k]['quin'][1] == 0:
            zero_row = newmaster[k]
            rowsel = k
    if rowsel == -1:
        raise Exception("Zero polarization not read")
    else:
        newmaster.remove_row(rowsel)

    random_pols = [random.uniform(0.0, plim) for _ in range(0, numfields)]
    random_norms = np.subtract(1.0, random_pols)
    random_samples = [random.randint(0, len(newmaster)-1) for _ in range(0, numfields)]
    random_qu_in = np.array([np.multiply(random_pols[k],
                                         newmaster['quin'][random_samples[k]]) for k in range(0, numfields)])
    random_ts = [[(random_norms[k] * zero_row['trans1']) + (random_pols[k] * newmaster['trans1'][random_samples[k]]),
                  (random_norms[k] * zero_row['trans2']) + (random_pols[k] * newmaster['trans2'][random_samples[k]]),
                  (random_norms[k] * zero_row['trans3']) + (random_pols[k] * newmaster['trans3'][random_samples[k]]),
                  (random_norms[k] * zero_row['trans4']) + (random_pols[k] * newmaster['trans4'][random_samples[k]])]
                 for k in range(0, numfields)]
    random_quout = [[(k[2]-k[0])/(k[2]+k[0]), (k[1]-k[3])/(k[3]+k[1])] for k in random_ts]

    x, y, z = symbols('x, y, z', real=True)

    if verbose >= 4:
        plt.clf()
        plt.plot(random_qu_in[:, 0], random_qu_in[:, 1], '.')
        plt.plot(ydata_q, ydata_u, '.')
        plt.show()

    if order == 2:
        if mixed:
            random_trans1 = (popt_trans1[0] + (popt_trans1[1] * x) + (popt_trans1[2] * y) + (popt_trans1[3] * x * x) +
                             (popt_trans1[4] * y * y) + (popt_trans1[5] * x * y))
            random_trans2 = (popt_trans2[0] + (popt_trans2[1] * x) + (popt_trans2[2] * y) + (popt_trans2[3] * x * x) +
                             (popt_trans2[4] * y * y) + (popt_trans2[5] * x * y))
            random_trans3 = (popt_trans3[0] + (popt_trans3[1] * x) + (popt_trans3[2] * y) + (popt_trans3[3] * x * x) +
                             (popt_trans3[4] * y * y) + (popt_trans3[5] * x * y))
            random_trans4 = (popt_trans4[0] + (popt_trans4[1] * x) + (popt_trans4[2] * y) + (popt_trans4[3] * x * x) +
                             (popt_trans4[4] * y * y) + (popt_trans4[5] * x * y))
            random_qu_out = np.array([(random_trans3-random_trans1)/(random_trans3+random_trans1),
                                      (random_trans2-random_trans4)/(random_trans2+random_trans4)])
            random_qu_modeled = []
            for k in range(0, numfields):
                print(f"Solving polarization {k+1}/{numfields}")
                random_qu_modeled.append(nonlinsolve([random_qu_out[0] - random_quout[k][0],
                                                      random_qu_out[1] - random_quout[k][1]], [x, y]))
            random_qu_modeled = np.array(random_qu_modeled)
        else:
            random_trans1 = (popt_trans1[0] + (popt_trans1[1] * x) + (popt_trans1[2] * y) + (popt_trans1[3] * x * x) +
                             (popt_trans1[4] * y * y))
            random_trans2 = (popt_trans2[0] + (popt_trans2[1] * x) + (popt_trans2[2] * y) + (popt_trans2[3] * x * x) +
                             (popt_trans2[4] * y * y))
            random_trans3 = (popt_trans3[0] + (popt_trans3[1] * x) + (popt_trans3[2] * y) + (popt_trans3[3] * x * x) +
                             (popt_trans3[4] * y * y))
            random_trans4 = (popt_trans4[0] + (popt_trans4[1] * x) + (popt_trans4[2] * y) + (popt_trans4[3] * x * x) +
                             (popt_trans4[4] * y * y))
            random_qu_out = np.array([(random_trans3-random_trans1)/(random_trans3+random_trans1),
                                      (random_trans2-random_trans4)/(random_trans2+random_trans4)])
            random_qu_modeled = []
            for k in range(0, numfields):
                print(f"Solving polarization {k+1}/{numfields}")
                random_qu_modeled.append(nonlinsolve([random_qu_out[0] - random_quout[k][0],
                                                      random_qu_out[1] - random_quout[k][1]], [x, y]))
            random_qu_modeled = np.array(random_qu_modeled)
    elif order == 3:
        random_trans1 = (popt_trans1[0] + (popt_trans1[1] * x) + (popt_trans1[2] * y) + (popt_trans1[3] * x * x) +
                         (popt_trans1[4] * y * y) + (popt_trans1[5] * x * y) + (popt_trans1[6] * x * x * x) +
                         (popt_trans1[7] * y * y * y) + (popt_trans1[8] * x * x * y) + (popt_trans1[9] * y * y * x))
        random_trans2 = (popt_trans2[0] + (popt_trans2[1] * x) + (popt_trans2[2] * y) + (popt_trans2[3] * x * x) +
                         (popt_trans2[4] * y * y) + (popt_trans2[5] * x * y) + (popt_trans2[6] * x * x * x) +
                         (popt_trans2[7] * y * y * y) + (popt_trans2[8] * x * x * y) + (popt_trans2[9] * y * y * x))
        random_trans3 = (popt_trans3[0] + (popt_trans3[1] * x) + (popt_trans3[2] * y) + (popt_trans3[3] * x * x) +
                         (popt_trans3[4] * y * y) + (popt_trans3[5] * x * y) + (popt_trans3[6] * x * x * x) +
                         (popt_trans3[7] * y * y * y) + (popt_trans3[8] * x * x * y) + (popt_trans3[9] * y * y * x))
        random_trans4 = (popt_trans4[0] + (popt_trans4[1] * x) + (popt_trans4[2] * y) + (popt_trans4[3] * x * x) +
                         (popt_trans4[4] * y * y) + (popt_trans4[5] * x * y) + (popt_trans4[6] * x * x * x) +
                         (popt_trans4[7] * y * y * y) + (popt_trans4[8] * x * x * y) + (popt_trans4[9] * y * y * x))
        random_qu_out = np.array([(random_trans3-random_trans1)/(random_trans3+random_trans1),
                                  (random_trans2-random_trans4)/(random_trans2+random_trans4)])
        random_qu_modeled = []
        for k in range(0, numfields):
            print(f"Solving polarization {k+1}/{numfields}")
            random_qu_modeled.append(nonlinsolve([random_qu_out[0] - random_quout[k][0],
                                                  random_qu_out[1] - random_quout[k][1]], [x, y]))
        random_qu_modeled = np.array(random_qu_modeled)
    else:
        raise Exception("Wrong order")

    bad_apples = 0
    temp = []
    for i in range(0, len(random_qu_modeled)):
        filt = [all([sympify(k).is_real for k in j]) for j in random_qu_modeled[i]]
        random_qu_modeled[i] = list(compress(random_qu_modeled[i], filt))
        try:
            distances = []
            for j in random_qu_modeled[i]:
                distances.append(dist(j, random_qu_in[i]))
            minpos = distances.index(min(distances))
            temp.append(random_qu_modeled[i][minpos])
        except:
            temp.append([np.nan, np.nan])
            bad_apples += 1
    random_qu_modeled = np.array(temp, dtype=np.float64)

    try:
        stdq = round(np.std(random_qu_in[:, 0]-random_qu_modeled[:, 0]), 7)
        stdu = round(np.std(random_qu_in[:, 1]-random_qu_modeled[:, 1]), 7)
        muq = round(np.mean(random_qu_in[:, 0]-random_qu_modeled[:, 0]), 7)
        muu = round(np.mean(random_qu_in[:, 1]-random_qu_modeled[:, 1]), 7)
    except Exception as e:
        print(e)
        stdq = np.nan
        stdu = np.nan
        muq = np.nan
        muu = np.nan

    if verbose >= 2:
        print(f'STD-q: {stdq}')
        print(f'STD-u: {stdu}')
        print(f'Mean-q: {muq}')
        print(f'Mean-u: {muu}')
        print("-----------------------------")

    result = [testfield[0], testfield[1], muq, muu, stdq, stdu, bad_apples]

    return result


if __name__ == '__main__':
    avail_modes = {
        'full': False,
        'model': True
    }
    header_vals = {
        False: ['x', 'y', 'mdq', 'mdu', 'sdq', 'sdu', 'bad_apples'],
        True: ['x', 'y',
               'popt_trans1', 'popt_trans2', 'popt_trans3', 'popt_trans4',
               'pcov_trans1', 'pcov_trans2', 'pcov_trans3', 'pcov_trans4']
    }

    try:
        _order = int(sys.argv[1])
    except:
        _order = 2

    try:
        _mixed = bool(int(sys.argv[2]))
    except:
        _mixed = True

    try:
        _tests = int(sys.argv[3])
    except:
        _tests = 1000

    try:
        _plim = float(sys.argv[4])
    except:
        _plim = 1.0

    try:
        _mode = str(sys.argv[5])
        if _mode not in avail_modes.keys():
            raise Exception("Mode not recognised")
    except:
        _mode = "full"

    fn = f'Results/bw_random{_tests}_order{_order}_mixed{_mixed}_lim{_plim}_mode{_mode}.csv'
    print(fn)

    with open(fn, 'w') as filer:
        write = csv.writer(filer)
        write.writerow(header_vals[avail_modes[_mode]])

    pool = mp.Pool(mp.cpu_count())

    jobs = []
    for looper in range(0, 144):
        job = pool.apply_async(testerfunc, (looper, _tests, 0, _order, _mixed, _plim, avail_modes[_mode]))
        jobs.append(job)

    results = []
    for job in jobs:
        results.append(job.get())

    pool.close()
    pool.join()

    with open(fn, 'a') as f:
        write = csv.writer(f)
        write.writerows(results)
