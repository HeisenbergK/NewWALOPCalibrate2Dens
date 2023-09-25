import numpy as np
from polarization import qu


def read_zmx(filename, wavelengthno=7):
    with open(filename, encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()[7:9]

    lines = [i.strip("\n") for i in lines]

    xfield = float(lines[0].split()[2])
    yfield = float(lines[1].split()[2])

    q, u = qu(xfield, yfield)

    fields = np.full([12, 2], -9999.0)
    tottransmissions = np.full(12, -9999.0)

    with open(filename, encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()[16:]
    lines = [line.strip("\n") for line in lines]
    for i in range(0, 12):
        for j in range(0, wavelengthno + 3):
            linetoexam = i * (wavelengthno + 3) + j
            # print("%d %d %d %s" %(i,j,linetoexam,lines[linetoexam]))
            if j == 0:
                # print(lines[linetoexam].split())
                coord1 = float(lines[linetoexam].split()[3].strip(','))
                coord2 = float(lines[linetoexam].split()[4])
                fields[i][0] = coord1
                fields[i][1] = coord2
            elif j == 8:
                # print(lines[linetoexam].split())
                tottransmissions[i] = (float(lines[linetoexam].split()[3]))
            elif j == 9:
                # print(lines[linetoexam].split())
                pass
            else:
                # print(lines[linetoexam].split())
                pass
    return q, u, fields, tottransmissions, xfield, yfield


def read_zmx_2(filename, wavelengthno):
    wavelengths = np.zeros(wavelengthno)
    fields = np.full([12, 2], -9999.0)
    tottransmissions = np.full(12, -9999.0)
    transmissions = np.full([12, wavelengthno], -9999.0)

    with open(filename, encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()[13:]
    lines = [line.strip("\n") for line in lines]
    for i in range(0, 12):
        for j in range(0, wavelengthno + 3):
            linetoexam = i * (wavelengthno + 3) + j
            # print("%d %d %d %s" %(i,j,linetoexam,lines[linetoexam]))
            if j == 0:
                # print(lines[linetoexam].split())
                coord1 = float(lines[linetoexam].split()[3].strip(','))
                coord2 = float(lines[linetoexam].split()[4])
                fields[i][0] = coord1
                fields[i][1] = coord2
            elif j == 8:
                # print(lines[linetoexam].split())
                tottransmissions[i] = (float(lines[linetoexam].split()[3]))
            elif j == 9:
                # print(lines[linetoexam].split())
                pass
            else:
                # print(lines[linetoexam].split())
                wav = float(lines[linetoexam].split()[2].strip(":"))
                trans = float(lines[linetoexam].split()[3])
                if i == 0:
                    wavelengths[j - 1] = wav
                transmissions[i][j - 1] = trans
    return 0, 0, fields, tottransmissions, 0, 0


if __name__ == '__main__':
    _wavelengths, _fields, _tottransmissions, _transmissions = read_zmx("Unpol/Config_1/transmission_field1_unpol.txt",
                                                                        wavelengthno=7)
    print(_fields)
    print(_tottransmissions)
    print(_wavelengths)
    print(_transmissions)
    print("-----------------------------------------------------------")

    _q, _u, _fields, _tottransmissions, _fieldx, _fieldy = read_zmx_2("Unpol/Config_1/transmission_field1_unpol.txt",
                                                                      wavelengthno=7)
    print(_q)
    print(_u)
    print(_fields)
    print(_tottransmissions)
    print(_fieldx)
    print(_fieldy)
