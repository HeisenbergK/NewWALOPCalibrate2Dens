import csv
import multiprocessing as mp
from read_zemax import read_zmx, read_zmx_2
from os import walk
from astropy.table import Table, unique
from scipy.optimize import curve_fit
import platform
import warnings
import random
import matplotlib.pyplot as plt
from models import *
import argparse
from splash import splash
from polarization import EVPAdeg, pol2cart, trans_to_pol
from os.path import isdir, isfile, join, islink
from os import listdir, unlink, mkdir
from shutil import rmtree

warnings.filterwarnings("ignore")

def testerfunc(
    testfieldcase,
    numfields=None,
    verbose=5,
    order=2,
    mixed=True,
    plim=1.0,
    ret_model=False,
    refit=0,
    noise=0,
    mkplots = False
):
    refitopt = refit
    refittrack = 0
    if refitopt=="auto":
        refit = 1

    columnames = ["config", "qu", "xy", "trans"]
    newcolumnames = ["quin", "evpa", "trans1", "trans2", "trans3", "trans4", "quout"]

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
        filelist = [
            ["Config_1\\%s" % i for i in filenames1],
            ["Config_2\\%s" % i for i in filenames2],
            ["Config_3\\%s" % i for i in filenames3],
            ["Config_4\\%s" % i for i in filenames4],
        ]
        unpolfilelist = [
            ["Unpol\\Config_1\\%s" % i for i in unpolfilenames1],
            ["Unpol\\Config_2\\%s" % i for i in unpolfilenames2],
            ["Unpol\\Config_3\\%s" % i for i in unpolfilenames3],
            ["Unpol\\Config_4\\%s" % i for i in unpolfilenames4],
        ]
    else:
        filelist = [
            ["Config_1/%s" % i for i in filenames1],
            ["Config_2/%s" % i for i in filenames2],
            ["Config_3/%s" % i for i in filenames3],
            ["Config_4/%s" % i for i in filenames4],
        ]
        unpolfilelist = [
            ["Unpol/Config_1/%s" % i for i in unpolfilenames1],
            ["Unpol/Config_2/%s" % i for i in unpolfilenames2],
            ["Unpol/Config_3/%s" % i for i in unpolfilenames3],
            ["Unpol/Config_4/%s" % i for i in unpolfilenames4],
        ]

    # make them in a table
    rowlist = []
    for i in range(0, 4):
        for j in range(0, len(filelist[i])):
            q, u, fields, tottransmissions, fieldx, fieldy = read_zmx(
                filelist[i][j], wavelengthno=7
            )
            print(filelist[i][j],q,u,fields,tottransmissions,fieldx, fieldy)
            for k in range(0, len(fields)):
                rowlist.append(
                    [i, [q, u], [fields[k][0], fields[k][1]], tottransmissions[k]]
                )
        for j in range(0, len(unpolfilelist[i])):
            q, u, fields, tottransmissions, fieldx, fieldy = read_zmx_2(
                unpolfilelist[i][j], wavelengthno=7
            )
            for k in range(0, len(fields)):
                rowlist.append(
                    [i, [q, u], [fields[k][0], fields[k][1]], tottransmissions[k]]
                )

    master = Table(list(map(list, zip(*rowlist))), names=columnames)

    # get the unique
    uniquefields = unique(master, keys="xy")["xy"].tolist()
    uniquepol = unique(master, keys="qu")["qu"].tolist()
    testfield = uniquefields[testfieldcase]
    mkdir(f"Modelplots/{testfield[0]}_{testfield[1]}")

    if verbose >= 3:
        print("-----------------------------")
        print(f"{len(uniquefields)} unique field positions")
        print("-----------------------------")
        print(f"{len(uniquepol)} unique polarizations")
        print("-----------------------------")

    slave = master[np.all(master["xy"] == testfield, axis=1)]
    slave.remove_column("xy")

    rowlist = []
    for i in range(0, len(uniquepol)):
        current_qu = uniquepol[i]
        current_evpa = EVPAdeg(*current_qu)
        trans1 = False
        trans2 = False
        trans3 = False
        trans4 = False
        for j in range(0, len(slave)):
            if slave["qu"][j].tolist() == current_qu:
                if slave["config"][j] == 0:
                    trans1 = slave["trans"][j]
                if slave["config"][j] == 1:
                    trans2 = slave["trans"][j]
                if slave["config"][j] == 2:
                    trans3 = slave["trans"][j]
                if slave["config"][j] == 3:
                    trans4 = slave["trans"][j]
        if not (trans1 and trans2 and trans3 and trans4):
            raise Exception(f"Transmittance not found for {current_qu}")
        quout = trans_to_pol(trans1, trans2, trans3, trans4)
        rowlist.append([current_qu, current_evpa, trans1, trans2, trans3, trans4, quout])

    newmaster = Table(list(map(list, zip(*rowlist))), names=newcolumnames)

    # fit a model
    xdata = newmaster["quout"]
    ydata = newmaster["quin"]
    ydata_q = ydata[:, 0]
    ydata_u = ydata[:, 1]
    if noise:
        noise_q = [np.random.normal(loc=0.0, scale=abs(noise*(i / 100))) for i in ydata_q]
        noise_u = [np.random.normal(loc=0.0, scale=abs(noise*(i / 100))) for i in ydata_u]
        ydata_q = np.add(noise_q, ydata_q)
        ydata_u = np.add(noise_u, ydata_u)
    if order == 2:
        if mixed:
            popt_q, pcov_q = curve_fit(modell, xdata, ydata_q)
            popt_u, pcov_u = curve_fit(modell, xdata, ydata_u)
            popt_trans1, pcov_trans1 = curve_fit(modell, ydata, newmaster["trans1"])
            popt_trans2, pcov_trans2 = curve_fit(modell, ydata, newmaster["trans2"])
            popt_trans3, pcov_trans3 = curve_fit(modell, ydata, newmaster["trans3"])
            popt_trans4, pcov_trans4 = curve_fit(modell, ydata, newmaster["trans4"])
            if mkplots:
                thetarange = np.linspace(0.0, 2.0*np.pi, 1000)
                qu_in_plot = np.array(
                    [
                        [
                            np.cos(2.0 * thetarange[i]),
                            np.sin(2.0 * thetarange[i]),
                        ]
                        for i in range(0, len(thetarange))
                    ]
                )
                plot_trans1 = [
                    modellrun(qu_in_plot[i], *popt_trans1) for i in range(0, len(qu_in_plot))
                ]
                plot_trans2 = [
                    modellrun(qu_in_plot[i], *popt_trans2) for i in range(0, len(qu_in_plot))
                ]
                plot_trans3 = [
                    modellrun(qu_in_plot[i], *popt_trans3) for i in range(0, len(qu_in_plot))
                ]
                plot_trans4 = [
                    modellrun(qu_in_plot[i], *popt_trans4) for i in range(0, len(qu_in_plot))
                ]
                plot_qu_out = np.array(
                    [
                        trans_to_pol(plot_trans1[i], plot_trans2[i], plot_trans3[i], plot_trans4[i])
                        for i in range(0, len(qu_in_plot))
                    ]
                )
                plot_qu_modeled = np.array(
                    [
                        [
                            modellrun(plot_qu_out[i], *popt_q),
                            modellrun(plot_qu_out[i], *popt_u),
                        ]
                        for i in range(0, len(qu_in_plot))
                    ]
                )
                plt.clf()
                plt.plot([EVPAdeg(ydata[i][0], ydata[i][1]) for i in range(0,len(ydata))], newmaster["trans1"], '.', label='Transmittance 1')
                plt.plot([EVPAdeg(ydata[i][0], ydata[i][1]) for i in range(0,len(ydata))], newmaster["trans2"], '.', label='Transmittance 2')
                plt.plot([EVPAdeg(ydata[i][0], ydata[i][1]) for i in range(0,len(ydata))], newmaster["trans3"], '.', label='Transmittance 3')
                plt.plot([EVPAdeg(ydata[i][0], ydata[i][1]) for i in range(0,len(ydata))], newmaster["trans4"], '.', label='Transmittance 4')
                plt.xlabel("EVPA(deg)")
                plt.ylabel("Transmittance")
                plt.legend()
                plt.savefig(f"Modelplots/{testfield[0]}_{testfield[1]}/transvsevpa.png")
                plt.clf()
                plt.plot([EVPAdeg(qu_in_plot[i][0], qu_in_plot[i][1]) for i in range(0,len(qu_in_plot))], plot_qu_out[:,0],'.', label="data-interp")
                plt.plot([EVPAdeg(ydata[i][0], ydata[i][1]) for i in range(0,len(ydata))], xdata[:, 0], '.', label="data")
                plt.xlabel(r'EVPA(deg)')
                plt.ylabel('q-data')
                plt.legend()
                plt.savefig(f"Modelplots/{testfield[0]}_{testfield[1]}/qdata.png")
                plt.clf()
                plt.plot([EVPAdeg(qu_in_plot[i][0], qu_in_plot[i][1]) for i in range(0,len(qu_in_plot))], plot_qu_out[:,1],'.', label="data-interp")
                plt.plot([EVPAdeg(ydata[i][0], ydata[i][1]) for i in range(0,len(ydata))], xdata[:, 1], '.', label="data")
                plt.xlabel(r'EVPA(deg)')
                plt.ylabel('u-data')
                plt.legend()
                plt.savefig(f"Modelplots/{testfield[0]}_{testfield[1]}/udata.png")
                plt.clf()
                plt.plot([EVPAdeg(qu_in_plot[i][0], qu_in_plot[i][1]) for i in range(0,len(qu_in_plot))], plot_qu_modeled[:,0],'.', label="model")
                plt.plot([EVPAdeg(ydata[i][0], ydata[i][1]) for i in range(0,len(ydata))], ydata[:, 0], '.', label="data")
                plt.xlabel(r'EVPA(deg)')
                plt.ylabel('q')
                plt.legend()
                plt.savefig(f"Modelplots/{testfield[0]}_{testfield[1]}/q.png")
                plt.clf()
                plt.plot([EVPAdeg(qu_in_plot[i][0], qu_in_plot[i][1]) for i in range(0,len(qu_in_plot))], plot_qu_modeled[:,1],'.', label="model")
                plt.plot([EVPAdeg(ydata[i][0], ydata[i][1]) for i in range(0,len(ydata))], ydata[:, 1], '.', label="data")
                plt.xlabel(r'EVPA(deg)')
                plt.ylabel('u')
                plt.legend()
                plt.savefig(f"Modelplots/{testfield[0]}_{testfield[1]}/u.png")
                plt.clf()
                plt.plot(qu_in_plot[:,0],plot_qu_modeled[:,0],'.')
                plt.xlabel('q_in')
                plt.ylabel('q_model')
                plt.savefig(f"Modelplots/{testfield[0]}_{testfield[1]}/qvq.png")
                plt.clf()
                plt.plot(qu_in_plot[:,1],plot_qu_modeled[:,1],'.')
                plt.xlabel('u_in')
                plt.ylabel('u_model')
                plt.savefig(f"Modelplots/{testfield[0]}_{testfield[1]}/uvu.png")
                plt.clf()
                plt.plot(plot_qu_out[:,0],plot_qu_modeled[:,0],'.', label="model")
                plt.plot(xdata[:,0],ydata[:,0],'.', label="data")
                plt.xlabel('q_out')
                plt.ylabel('q')
                plt.legend()
                plt.savefig(f"Modelplots/{testfield[0]}_{testfield[1]}/qvqout.png")
                plt.clf()
                plt.plot(plot_qu_out[:,1],plot_qu_modeled[:,1],'.', label="model")
                plt.plot(xdata[:,1],ydata[:,1],'.', label="data")
                plt.xlabel('u_out')
                plt.ylabel('u')
                plt.legend()
                plt.savefig(f"Modelplots/{testfield[0]}_{testfield[1]}/uvuout.png")
                plt.clf()
                plt.plot(plot_qu_out[:,1],plot_qu_modeled[:,0],'.', label="model")
                plt.plot(xdata[:,1],ydata[:,0],'.', label="data")
                plt.xlabel('u_out')
                plt.ylabel('q')
                plt.legend()
                plt.savefig(f"Modelplots/{testfield[0]}_{testfield[1]}/qvuout.png")
                plt.clf()
                plt.plot(plot_qu_out[:,0],plot_qu_modeled[:,1],'.', label="model")
                plt.plot(xdata[:,0],ydata[:,1],'.', label="data")
                plt.xlabel('u_out')
                plt.ylabel('q')
                plt.legend()
                plt.savefig(f"Modelplots/{testfield[0]}_{testfield[1]}/uvqout.png")
            if ret_model:
                return [
                    testfield[0],
                    testfield[1],
                    popt_q.tolist(),
                    popt_u.tolist(),
                    pcov_q.tolist(),
                    pcov_u.tolist(),
                ]
        else:
            popt_q, pcov_q = curve_fit(modellnm, xdata, ydata_q)
            popt_u, pcov_u = curve_fit(modellnm, xdata, ydata_u)
            popt_trans1, pcov_trans1 = curve_fit(modellnm, ydata, newmaster["trans1"])
            popt_trans2, pcov_trans2 = curve_fit(modellnm, ydata, newmaster["trans2"])
            popt_trans3, pcov_trans3 = curve_fit(modellnm, ydata, newmaster["trans3"])
            popt_trans4, pcov_trans4 = curve_fit(modellnm, ydata, newmaster["trans4"])
            if ret_model:
                return [
                    testfield[0],
                    testfield[1],
                    popt_q.tolist(),
                    popt_u.tolist(),
                    pcov_q.tolist(),
                    pcov_u.tolist(),
                ]
    elif order == 3:
        popt_q, pcov_q = curve_fit(modell3, xdata, ydata_q)
        popt_u, pcov_u = curve_fit(modell3, xdata, ydata_u)
        popt_trans1, pcov_trans1 = curve_fit(modell3, ydata, newmaster["trans1"])
        popt_trans2, pcov_trans2 = curve_fit(modell3, ydata, newmaster["trans2"])
        popt_trans3, pcov_trans3 = curve_fit(modell3, ydata, newmaster["trans3"])
        popt_trans4, pcov_trans4 = curve_fit(modell3, ydata, newmaster["trans4"])
        if ret_model:
            return [
                testfield[0],
                testfield[1],
                popt_q.tolist(),
                popt_u.tolist(),
                pcov_q.tolist(),
                pcov_u.tolist(),
            ]
    else:
        raise Exception("Wrong order")

    random_pols = [random.uniform(0.0, plim) for _ in range(0, numfields)]
    random_angles = [np.pi * random.random() for _ in range(0, numfields)]
    random_qu_in = np.array(
        [
            [
                random_pols[i] * np.cos(2.0 * random_angles[i]),
                random_pols[i] * np.sin(2.0 * random_angles[i]),
            ]
            for i in range(0, numfields)
        ]
    )

    if verbose >= 4:
        plt.clf()
        plt.plot(random_qu_in[:, 0], random_qu_in[:, 1], ".")
        plt.plot(ydata_q, ydata_u, ".")
        plt.show()

    if order == 2:
        if mixed:
            random_trans1 = [
                modellrun(random_qu_in[i], *popt_trans1) for i in range(0, numfields)
            ]
            random_trans2 = [
                modellrun(random_qu_in[i], *popt_trans2) for i in range(0, numfields)
            ]
            random_trans3 = [
                modellrun(random_qu_in[i], *popt_trans3) for i in range(0, numfields)
            ]
            random_trans4 = [
                modellrun(random_qu_in[i], *popt_trans4) for i in range(0, numfields)
            ]
            random_qu_out = np.array(
                [
                    trans_to_pol(random_trans1[i], random_trans2[i], random_trans3[i], random_trans4[i])
                    for i in range(0, numfields)
                ]
            )
            random_qu_modeled = np.array(
                [
                    [
                        modellrun(random_qu_out[i], *popt_q),
                        modellrun(random_qu_out[i], *popt_u),
                    ]
                    for i in range(0, numfields)
                ]
            )
        else:
            random_trans1 = [
                modellnmrun(random_qu_in[i], *popt_trans1) for i in range(0, numfields)
            ]
            random_trans2 = [
                modellnmrun(random_qu_in[i], *popt_trans2) for i in range(0, numfields)
            ]
            random_trans3 = [
                modellnmrun(random_qu_in[i], *popt_trans3) for i in range(0, numfields)
            ]
            random_trans4 = [
                modellnmrun(random_qu_in[i], *popt_trans4) for i in range(0, numfields)
            ]
            random_qu_out = np.array(
                [
                    trans_to_pol(random_trans1[i], random_trans2[i], random_trans3[i], random_trans4[i])
                    for i in range(0, numfields)
                ]
            )
            random_qu_modeled = np.array(
                [
                    [
                        modellnmrun(random_qu_out[i], *popt_q),
                        modellnmrun(random_qu_out[i], *popt_u),
                    ]
                    for i in range(0, numfields)
                ]
            )
    elif order == 3:
        random_trans1 = [
            modell3run(random_qu_in[i], *popt_trans1) for i in range(0, numfields)
        ]
        random_trans2 = [
            modell3run(random_qu_in[i], *popt_trans2) for i in range(0, numfields)
        ]
        random_trans3 = [
            modell3run(random_qu_in[i], *popt_trans3) for i in range(0, numfields)
        ]
        random_trans4 = [
            modell3run(random_qu_in[i], *popt_trans4) for i in range(0, numfields)
        ]
        random_qu_out = np.array(
            [
                trans_to_pol(random_trans1[i], random_trans2[i], random_trans3[i], random_trans4[i])
                for i in range(0, numfields)
            ]
        )
        random_qu_modeled = np.array(
            [
                [
                    modell3run(random_qu_out[i], *popt_q),
                    modell3run(random_qu_out[i], *popt_u),
                ]
                for i in range(0, numfields)
            ]
        )
    else:
        raise Exception("Wrong order")

    qdiff = random_qu_in[:, 0] - random_qu_modeled[:, 0]
    udiff = random_qu_in[:, 1] - random_qu_modeled[:, 1]

    if verbose >= 2:
        print(f"STD-q: {round(float(np.std(qdiff)), 7)}")
        print(f"STD-u: {round(float(np.std(udiff)), 7)}")
        print(f"Mean-q: {round(float(np.mean(qdiff)), 7)}")
        print(f"Mean-u: {round(float(np.mean(udiff)), 7)}")
        print(f"Covariance Matrix :\n\t\t\t{np.cov(qdiff, udiff)}")
        print("-----------------------------")

    result = [
        testfield[0],
        testfield[1],
        round(float(np.mean(qdiff)), 7),
        round(float(np.mean(udiff)), 7),
        round(float(np.std(qdiff)), 7),
        round(float(np.std(udiff)), 7),
        0,
        np.cov(qdiff, udiff),
    ]

    if not refit:
        return result

    qcorr = popt_q[0] - result[2]
    ucorr = popt_u[0] - result[3]
    laststdq = np.std(qdiff)
    laststdu = np.std(udiff)
    first_time = True
    # refit a model
    while refit:
        if not first_time:
            qcorr = qcorr - result[2]
            ucorr = ucorr - result[3]
        first_time = False
        if refitopt!="auto":
            refit -= 1
        xdata = newmaster["quout"]
        xdatanewq = []
        xdatanewu = []
        for i in range(0, len(xdata)):
            xdatanewq.append(np.append(xdata[i], qcorr))
            xdatanewu.append(np.append(xdata[i], ucorr))
        xdatanewq = np.array(xdatanewq)
        xdatanewu = np.array(xdatanewu)
        ydata = newmaster["quin"]
        ydata_q = ydata[:, 0]
        ydata_u = ydata[:, 1]
        if order == 2:
            if mixed:
                popt_q, pcov_q = curve_fit(modellconst, xdatanewq, ydata_q, method='dogbox')
                popt_u, pcov_u = curve_fit(modellconst, xdatanewu, ydata_u, method='dogbox')
                popt_trans1, pcov_trans1 = curve_fit(modell, ydata, newmaster["trans1"])
                popt_trans2, pcov_trans2 = curve_fit(modell, ydata, newmaster["trans2"])
                popt_trans3, pcov_trans3 = curve_fit(modell, ydata, newmaster["trans3"])
                popt_trans4, pcov_trans4 = curve_fit(modell, ydata, newmaster["trans4"])
                if ret_model:
                    return [
                        testfield[0],
                        testfield[1],
                        popt_q.tolist(),
                        popt_u.tolist(),
                        pcov_q.tolist(),
                        pcov_u.tolist(),
                    ]
            else:
                popt_q, pcov_q = curve_fit(modellnmconst, xdatanewq, ydata_q)
                popt_u, pcov_u = curve_fit(modellnmconst, xdatanewu, ydata_u)
                popt_trans1, pcov_trans1 = curve_fit(
                    modellnm, ydata, newmaster["trans1"]
                )
                popt_trans2, pcov_trans2 = curve_fit(
                    modellnm, ydata, newmaster["trans2"]
                )
                popt_trans3, pcov_trans3 = curve_fit(
                    modellnm, ydata, newmaster["trans3"]
                )
                popt_trans4, pcov_trans4 = curve_fit(
                    modellnm, ydata, newmaster["trans4"]
                )
                if ret_model:
                    return [
                        testfield[0],
                        testfield[1],
                        popt_q.tolist(),
                        popt_u.tolist(),
                        pcov_q.tolist(),
                        pcov_u.tolist(),
                    ]
        elif order == 3:
            popt_q, pcov_q = curve_fit(modell3const, xdatanewq, ydata_q)
            popt_u, pcov_u = curve_fit(modell3const, xdatanewu, ydata_u)
            popt_trans1, pcov_trans1 = curve_fit(modell3, ydata, newmaster["trans1"])
            popt_trans2, pcov_trans2 = curve_fit(modell3, ydata, newmaster["trans2"])
            popt_trans3, pcov_trans3 = curve_fit(modell3, ydata, newmaster["trans3"])
            popt_trans4, pcov_trans4 = curve_fit(modell3, ydata, newmaster["trans4"])
            if ret_model:
                return [
                    testfield[0],
                    testfield[1],
                    popt_q.tolist(),
                    popt_u.tolist(),
                    pcov_q.tolist(),
                    pcov_u.tolist(),
                ]
        else:
            raise Exception("Wrong order")

        random_pols = [random.uniform(0.0, plim) for _ in range(0, numfields)]
        random_angles = [np.pi * random.random() for _ in range(0, numfields)]
        random_qu_in = np.array(
            [
                [
                    random_pols[i] * np.cos(2.0 * random_angles[i]),
                    random_pols[i] * np.sin(2.0 * random_angles[i]),
                ]
                for i in range(0, numfields)
            ]
        )

        if verbose >= 4:
            plt.clf()
            plt.plot(random_qu_in[:, 0], random_qu_in[:, 1], ".")
            plt.plot(ydata_q, ydata_u, ".")
            plt.show()

        if order == 2:
            if mixed:
                random_trans1 = [
                    modellrun(random_qu_in[i], *popt_trans1)
                    for i in range(0, numfields)
                ]
                random_trans2 = [
                    modellrun(random_qu_in[i], *popt_trans2)
                    for i in range(0, numfields)
                ]
                random_trans3 = [
                    modellrun(random_qu_in[i], *popt_trans3)
                    for i in range(0, numfields)
                ]
                random_trans4 = [
                    modellrun(random_qu_in[i], *popt_trans4)
                    for i in range(0, numfields)
                ]
                random_qu_out = np.array(
                    [
                        trans_to_pol(random_trans1[i], random_trans2[i], random_trans3[i], random_trans4[i])
                        for i in range(0, numfields)
                    ]
                )
                random_qu_modeled = np.array(
                    [
                        [
                            modellrunconst([*random_qu_out[i], qcorr], *popt_q),
                            modellrunconst([*random_qu_out[i], ucorr], *popt_u),
                        ]
                        for i in range(0, numfields)
                    ]
                )
            else:
                random_trans1 = [
                    modellnmrun(random_qu_in[i], *popt_trans1)
                    for i in range(0, numfields)
                ]
                random_trans2 = [
                    modellnmrun(random_qu_in[i], *popt_trans2)
                    for i in range(0, numfields)
                ]
                random_trans3 = [
                    modellnmrun(random_qu_in[i], *popt_trans3)
                    for i in range(0, numfields)
                ]
                random_trans4 = [
                    modellnmrun(random_qu_in[i], *popt_trans4)
                    for i in range(0, numfields)
                ]
                random_qu_out = np.array(
                    [
                        trans_to_pol(random_trans1[i], random_trans2[i], random_trans3[i], random_trans4[i])
                        for i in range(0, numfields)
                    ]
                )
                random_qu_modeled = np.array(
                    [
                        [
                            modellnmrunconst([*random_qu_out[i], qcorr], *popt_q),
                            modellnmrunconst([*random_qu_out[i], ucorr], *popt_u),
                        ]
                        for i in range(0, numfields)
                    ]
                )
        elif order == 3:
            random_trans1 = [
                modell3run(random_qu_in[i], *popt_trans1) for i in range(0, numfields)
            ]
            random_trans2 = [
                modell3run(random_qu_in[i], *popt_trans2) for i in range(0, numfields)
            ]
            random_trans3 = [
                modell3run(random_qu_in[i], *popt_trans3) for i in range(0, numfields)
            ]
            random_trans4 = [
                modell3run(random_qu_in[i], *popt_trans4) for i in range(0, numfields)
            ]
            random_qu_out = np.array(
                [
                    trans_to_pol(random_trans1[i], random_trans2[i], random_trans3[i], random_trans4[i])
                    for i in range(0, numfields)
                ]
            )
            random_qu_modeled = np.array(
                [
                    [
                        modell3runconst([*random_qu_out[i], qcorr], *popt_q),
                        modell3runconst([*random_qu_out[i], ucorr], *popt_u),
                    ]
                    for i in range(0, numfields)
                ]
            )
        else:
            raise Exception("Wrong order")

        qdiff = random_qu_in[:, 0] - random_qu_modeled[:, 0]
        udiff = random_qu_in[:, 1] - random_qu_modeled[:, 1]

        if verbose >= 2:
            print(f"STD-q: {round(float(np.std(qdiff)), 7)}")
            print(f"STD-u: {round(float(np.std(udiff)), 7)}")
            print(f"Mean-q: {round(float(np.mean(qdiff)), 7)}")
            print(f"Mean-u: {round(float(np.mean(udiff)), 7)}")
            print(f"Covariance Matrix :\n\t\t\t{np.cov(qdiff, udiff)}")
            print("-----------------------------")
        
        qdiffimpr = laststdq - np.std(qdiff)
        udiffimpr = laststdu - np.std(udiff)
        totimpr = qdiffimpr + udiffimpr

        if totimpr>0 and refitopt=="auto":
            result = [
                testfield[0],
                testfield[1],
                round(float(np.mean(qdiff)), 7),
                round(float(np.mean(udiff)), 7),
                round(float(np.std(qdiff)), 7),
                round(float(np.std(udiff)), 7),
                0,
                np.cov(qdiff, udiff),
            ]
            refittrack+=1
        elif totimpr<=0 and refitopt=="auto":
            refit=0
        else:
            result = [
                testfield[0],
                testfield[1],
                round(float(np.mean(qdiff)), 7),
                round(float(np.mean(udiff)), 7),
                round(float(np.std(qdiff)), 7),
                round(float(np.std(udiff)), 7),
                0,
                np.cov(qdiff, udiff),
            ]

    return result


if __name__ == "__main__":
    print(chr(27) + "[2J")

    if isdir("Modelplots"):
        for filename in listdir("Modelplots"):
            file_path = join("Modelplots", filename)
            if isfile(file_path) or islink(file_path):
                unlink(file_path)
            elif isdir(file_path):
                rmtree(file_path)
        print("Removed all contents of Modelplots")
    else:
        mkdir("Modelplots")

    splash()
    avail_modes = {"full": False, "model": True}
    header_vals = {
        False: ["x", "y", "mdq", "mdu", "sdq", "sdu", "bad_apples"],
        True: ["x", "y", "popt_q", "popt_u", "pcov_q", "pcov_u"],
    }
    header_cov = ["x", "y", "covariance"]

    parser = argparse.ArgumentParser("Model the WALOP instrument - forward method")
    parser.add_argument("order", type=int, default=2, help="Fit order", choices=[2, 3])
    parser.add_argument(
        "mixed",
        type=int,
        default=1,
        help="Whether to use the mixed term with 2nd order fit",
        choices=[0, 1],
    )
    parser.add_argument("tests", type=int, default=100, help="Number of tests to run")
    parser.add_argument(
        "polarization_limit",
        type=float,
        default=1.0,
        help="Polarization limit for the tests in the range [0.001, 1.0]",
    )
    parser.add_argument(
        "refit",
        default=0,
        help="Whether to refit and how many times, freezing the constant terms",
        choices=["0","1","2","3","4","5","6","7","8","9","10","auto"],
        metavar="[0-10/auto]",
    )
    parser.add_argument(
        "noise",
        type=float,
        default=0.0,
        help="Noise for the model in %%",
    )
    parser.add_argument(
        "--model",
        help="Add this flag to only save the model parameters instead of the model std",
        action="store_true",
    )
    parser.add_argument(
        "--plots",
        help="Add this flag to save the plots of modelled vs input",
        action="store_true",
    )
    args = parser.parse_args()

    _order = args.order
    _mixed = bool(args.mixed)
    _tests = args.tests
    _plim = args.polarization_limit
    _noise = args.noise
    if _noise==0.0:
        _noise = False
    if args.refit=="auto":
        _refit="auto"
    else:
        _refit = int(args.refit)
    if _refit==1 or _refit==0:
        _refit=bool(_refit)
    if args.model:
        _mode = "model"
    else:
        _mode = "full"
    if args.plots:
        _mkplots = True
    else:
        _mkplots = False

    if _plim < 0.001 or _plim > 1.0:
        print("Polarization limit cannot be outside the range [0.001, 1.0]")
        exit()

    print(f"Running with the following parameters:")
    print(f"Order:\t{_order}")
    print(f"Mixed:\t{_mixed}")
    print(f"Tests:\t{_tests}")
    print(f"Polarization Limit:\t{_plim}")
    print(f"Refit:\t{_refit}")
    print(f'Noise:\t{_noise}%')
    print(f"Mode:\t{_mode}")
    print(f"Make Plots:\t{_mkplots}")
    goon = input("Press enter to continue, Ctrl+c to abort")

    fn = f"Results/fw_random{_tests}_order{_order}_mixed{_mixed}_lim{_plim}_mode{_mode}_refit{_refit}_noise{_noise}.csv"
    fn2 = f"Results/fw_random{_tests}_order{_order}_mixed{_mixed}_lim{_plim}_mode{_mode}_refit{_refit}_noise{_noise}_covs.csv"
    print(fn)

    with open(fn, "w") as filer:
        write = csv.writer(filer)
        write.writerow(header_vals[avail_modes[_mode]])

    if not avail_modes[_mode]:
        with open(fn2, "w") as filer:
            write = csv.writer(filer)
            write.writerow(header_cov)

    pool = mp.Pool(mp.cpu_count())

    jobs = []
    for looper in range(0, 576):
        job = pool.apply_async(
            testerfunc,
            (looper, _tests, 0, _order, _mixed, _plim, avail_modes[_mode], _refit, _noise, _mkplots),
        )
        jobs.append(job)

    results = []
    covs = []
    for job in jobs:
        _result = job.get()
        _covariance = _result.pop()
        _covariance = [_result[0], _result[1], _covariance]
        results.append(_result)
        covs.append(_covariance)

    pool.close()
    pool.join()

    with open(fn, "a") as f:
        write = csv.writer(f)
        write.writerows(results)

    if not avail_modes[_mode]:
        with open(fn2, "a") as f:
            write = csv.writer(f)
            write.writerows(covs)
