from read_zemax import read_zmx, read_zmx_2
from os import walk
from astropy.table import Table, unique, vstack
import platform
import warnings
from models import *
from splash import splash
from polarization import EVPAdeg, trans_to_pol
from tqdm import tqdm
import multiprocessing as mp

warnings.filterwarnings("ignore")


def testerfunc(testfieldcase, verbose=5):
    columnames = ["config", "qu", "xy", "trans"]
    newcolumnames = [
        "x",
        "y",
        "qin",
        "uin",
        "evpa",
        "trans1",
        "trans2",
        "trans3",
        "trans4",
        "qout",
        "uout",
    ]

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
            if verbose > 2:
                print(filelist[i][j], q, u, fields, tottransmissions, fieldx, fieldy)
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
        rowlist.append(
            [
                testfield[0],
                testfield[1],
                current_qu[0],
                current_qu[1],
                current_evpa,
                trans1,
                trans2,
                trans3,
                trans4,
                quout[0],
                quout[1],
            ]
        )

    newmaster = Table(list(map(list, zip(*rowlist))), names=newcolumnames)
    return newmaster


if __name__ == "__main__":
    print(chr(27) + "[2J")

    splash()

    print(f"Running the extractor")
    goon = input("Press enter to continue, Ctrl+c to abort")

    jobs = []

    pool = mp.Pool(mp.cpu_count())

    for looper in range(576):
        job = pool.apply_async(
            testerfunc,
            (looper, 0),
        )
        jobs.append(job)

    pool.close()
    pool.join()

    tabs = []

    for job in jobs:
        _tab = job.get()
        tabs.append(_tab)

    mastertab = vstack([i for i in tabs])

    mastertab.write("Results/extracted.ecsv", format="ascii")
