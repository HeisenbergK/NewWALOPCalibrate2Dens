import numpy as np
import csv
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from splash import splash
import argparse
from models import modellrun
from astropy.table import Table
from astropy.io import ascii
import matplotlib.pyplot as plt

print(chr(27) + "[2J")
splash()

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
    choices=["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "auto"],
    metavar="[0-10/auto]",
)
parser.add_argument(
    "noise",
    type=float,
    default=0.0,
    help="Noise for the model in %%",
)
parser.add_argument(
    "satlim",
    type=float,
    default=1000,
    help="Saturation limit",
)
args = parser.parse_args()

_order = args.order
_mixed = bool(args.mixed)
_tests = args.tests
_plim = args.polarization_limit
_noise = args.noise
_satlim = args.satlim
if _noise == 0.0:
    _noise = False
if args.refit == "auto":
    _refit = "auto"
else:
    _refit = int(args.refit)
if _refit == 1 or _refit == 0:
    _refit = bool(_refit)

if _plim < 0.001 or _plim > 1.0:
    print("Polarization limit cannot be outside the range [0.001, 1.0]")
    exit()

if _satlim < 0:
    _satlim = 0.0 - _satlim

print(f"Running with the following parameters:")
print(f"Order:\t{_order}")
print(f"Mixed:\t{_mixed}")
print(f"Tests:\t{_tests}")
print(f"Polarization Limit:\t{_plim}")
print(f"Refit:\t{_refit}")
print(f"Noise:\t{_noise}%")
print(f"Saturation Limit:\t{_satlim}")
goon = input("Press enter to continue, Ctrl+c to abort")

filenam = f"Results/fw_random{_tests}_order{_order}_mixed{_mixed}_lim{_plim}_modemodel_refit{_refit}_noise{_noise}.csv"

master = []
with open(filenam) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=",")
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            print()
            line_count += 1
        else:
            master.append(
                [
                    float(row[0]),
                    float(row[1]),
                    [float(i) for i in row[2].split("[")[1].split("]")[0].split(", ")],
                    [float(i) for i in row[3].split("[")[1].split("]")[0].split(", ")],
                ]
            )
            line_count += 1
    print(f"Processed {line_count} lines.")

extracted = ascii.read("Results/extracted.ecsv")
extracted["qin"].format = "%.16f"
extracted["uin"].format = "%.16f"

unx = np.sort(np.unique(np.array(master).T[0]))
uny = np.sort(np.unique(np.array(master).T[1]))
popt_q_0 = np.empty((len(unx), len(uny)))
popt_q_q = np.empty((len(unx), len(uny)))
popt_q_u = np.empty((len(unx), len(uny)))
popt_q_q2 = np.empty((len(unx), len(uny)))
popt_q_u2 = np.empty((len(unx), len(uny)))
popt_q_qu = np.empty((len(unx), len(uny)))
popt_q_q3 = np.empty((len(unx), len(uny)))
popt_q_u3 = np.empty((len(unx), len(uny)))
popt_q_q2u = np.empty((len(unx), len(uny)))
popt_q_qu2 = np.empty((len(unx), len(uny)))
popt_u_0 = np.empty((len(unx), len(uny)))
popt_u_q = np.empty((len(unx), len(uny)))
popt_u_u = np.empty((len(unx), len(uny)))
popt_u_q2 = np.empty((len(unx), len(uny)))
popt_u_u2 = np.empty((len(unx), len(uny)))
popt_u_qu = np.empty((len(unx), len(uny)))
popt_u_q3 = np.empty((len(unx), len(uny)))
popt_u_u3 = np.empty((len(unx), len(uny)))
popt_u_q2u = np.empty((len(unx), len(uny)))
popt_u_qu2 = np.empty((len(unx), len(uny)))
popt_q_0[:] = np.nan
popt_q_q[:] = np.nan
popt_q_u[:] = np.nan
popt_q_q2[:] = np.nan
popt_q_u2[:] = np.nan
popt_q_qu[:] = np.nan
popt_q_q3[:] = np.nan
popt_q_u3[:] = np.nan
popt_q_q2u[:] = np.nan
popt_q_qu2[:] = np.nan
popt_u_0[:] = np.nan
popt_u_q[:] = np.nan
popt_u_u[:] = np.nan
popt_u_q2[:] = np.nan
popt_u_u2[:] = np.nan
popt_u_qu[:] = np.nan
popt_u_q3[:] = np.nan
popt_u_u3[:] = np.nan
popt_u_q2u[:] = np.nan
popt_u_qu2[:] = np.nan

q0q = np.empty((len(unx), len(uny)))
q1q = np.empty((len(unx), len(uny)))
q0u = np.empty((len(unx), len(uny)))
q1u = np.empty((len(unx), len(uny)))
q0q[:] = np.nan
q1q[:] = np.nan
q0u[:] = np.nan
q1u[:] = np.nan

for i in range(0, len(master)):
    curx = master[i][0]
    cury = master[i][1]
    indx = np.where(unx == curx)[0]
    indy = np.where(uny == cury)[0]
    partextracted1q0u = extracted[
        np.logical_and(
            np.logical_and(extracted["x"] == curx, extracted["y"] == cury),
            np.logical_and(extracted["qin"] == 1, extracted["uin"] == 0),
        )
    ]
    partextracted0q1u = extracted[
        np.logical_and(
            np.logical_and(extracted["x"] == curx, extracted["y"] == cury),
            np.logical_and(
                extracted["qin"] < 0.0000000000000002, extracted["uin"] == 1
            ),
        )
    ]

    try:
        popt_q_0[indy, indx] = master[i][2][0]
        popt_u_0[indy, indx] = master[i][3][0]

        popt_q_q[indy, indx] = master[i][2][1]
        popt_u_q[indy, indx] = master[i][3][1]

        popt_q_u[indy, indx] = master[i][2][2]
        popt_u_u[indy, indx] = master[i][3][2]

        popt_q_q2[indy, indx] = master[i][2][3]
        popt_u_q2[indy, indx] = master[i][3][3]

        popt_q_u2[indy, indx] = master[i][2][4]
        popt_u_u2[indy, indx] = master[i][3][4]

        popt_q_qu[indy, indx] = master[i][2][5]
        popt_u_qu[indy, indx] = master[i][3][5]

        q0qval = modellrun(
            [partextracted0q1u["qin"][0], partextracted0q1u["uin"][0]],
            popt_q_0[indy, indx],
            popt_q_q[indy, indx],
            popt_q_u[indy, indx],
            popt_q_q2[indy, indx],
            popt_q_u2[indy, indx],
            popt_q_qu[indy, indx],
        )

        if q0qval < 0:
            q0q[indy, indx] = abs(max(q0qval, -_satlim))
        else:
            q0q[indy, indx] = min(q0qval, _satlim)

        q1qval = modellrun(
            [partextracted1q0u["qin"][0], partextracted1q0u["uin"][0]],
            popt_q_0[indy, indx],
            popt_q_q[indy, indx],
            popt_q_u[indy, indx],
            popt_q_q2[indy, indx],
            popt_q_u2[indy, indx],
            popt_q_qu[indy, indx],
        )

        if q1qval < 0:
            q1q[indy, indx] = abs(max(q1qval, -_satlim))
        else:
            q1q[indy, indx] = min(q1qval, _satlim)

        q0uval = modellrun(
            [partextracted0q1u["qin"][0], partextracted0q1u["uin"][0]],
            popt_u_0[indy, indx],
            popt_u_q[indy, indx],
            popt_u_u[indy, indx],
            popt_u_q2[indy, indx],
            popt_u_u2[indy, indx],
            popt_u_qu[indy, indx],
        )

        if q0uval < 0:
            q0u[indy, indx] = abs(max(q0uval, -_satlim))
        else:
            q0u[indy, indx] = min(q0uval, _satlim)

        q1uval = modellrun(
            [partextracted1q0u["qin"][0], partextracted1q0u["uin"][0]],
            popt_u_0[indy, indx],
            popt_u_q[indy, indx],
            popt_u_u[indy, indx],
            popt_u_q2[indy, indx],
            popt_u_u2[indy, indx],
            popt_u_qu[indy, indx],
        )

        if q1uval < 0:
            q1u[indy, indx] = abs(max(q1uval, -_satlim))
        else:
            q1u[indy, indx] = min(q1uval, _satlim)

    except Exception as e:
        print(e)


fig = make_subplots(
    rows=2,
    cols=2,
    subplot_titles=(
        f"q for (q,u)in=(1,0)",
        f"q for (q,u)in=(0,1)",
        f"u for (q,u)in=(1,0)",
        f"u for (q,u)in=(0,1)",
    ),
    column_widths=[0.5, 0.5],
)

fig.add_trace(
    go.Heatmap(
        z=q1q,
        x=unx,  # horizontal axis
        y=uny,
        colorbar=dict(exponentformat="power", len=0.4, x=0.45, y=0.8),
        zsmooth="best",
    ),
    row=1,
    col=1,
)

fig.add_trace(
    go.Heatmap(
        z=q0q,
        x=unx,  # horizontal axis
        y=uny,
        colorbar=dict(exponentformat="power", len=0.4, x=1, y=0.8),
        zsmooth="best",
    ),
    row=1,
    col=2,
)

fig.add_trace(
    go.Heatmap(
        z=q1u,
        x=unx,  # horizontal axis
        y=uny,
        colorbar=dict(exponentformat="power", len=0.4, x=0.45, y=0.2),
        zsmooth="best",
    ),
    row=2,
    col=1,
)

fig.add_trace(
    go.Heatmap(
        z=q0u,
        x=unx,  # horizontal axis
        y=uny,
        colorbar=dict(exponentformat="power", len=0.4, x=1, y=0.2),
        zsmooth="best",
    ),
    row=2,
    col=2,
)

fig["layout"]["yaxis1"].update(scaleanchor="x1", scaleratio=1)
fig["layout"]["yaxis2"].update(scaleanchor="x2", scaleratio=1)
fig["layout"]["yaxis3"].update(scaleanchor="x3", scaleratio=1)
fig["layout"]["yaxis4"].update(scaleanchor="x4", scaleratio=1)
fig["layout"]["xaxis1"].update(range=[-0.25, 0.25])
fig["layout"]["yaxis2"].update(range=[-0.25, 0.25])
fig["layout"]["yaxis3"].update(range=[-0.25, 0.25])
fig["layout"]["yaxis4"].update(range=[-0.25, 0.25])
# fig.show()

malakia = [i for i in range(0,27,4)]

real_x=np.array([np.round(((i/24)*0.5)-0.25, 2) for i in malakia])
real_y=np.array([np.round(((i/24)*0.5)-0.25, 2) for i in malakia])

plt.clf()
plt.imshow(q0q, interpolation="quadric", extent=[0, 24, 0, 24])
plt.title(r"$q_c$ for $\left(q_i,u_i\right)=\left(0,1\right)$")
plt.xlabel("x(deg)")
plt.ylabel("y(deg)")
plt.gca().invert_yaxis()
plt.gca().set_xticks(malakia)
plt.gca().set_xticklabels(real_x)
plt.gca().set_yticks(malakia)
plt.gca().set_yticklabels(real_y)
plt.colorbar()
plt.tight_layout()
plt.savefig("Images/q0q.svg")

plt.clf()
plt.imshow(q0u, interpolation="quadric", extent=[0, 24, 0, 24])
plt.title(r"$u_c$ for $\left(q_i,u_i\right)=\left(0,1\right)$")
plt.xlabel("x(deg)")
plt.ylabel("y(deg)")
plt.gca().invert_yaxis()
plt.gca().set_xticks(malakia)
plt.gca().set_xticklabels(real_x)
plt.gca().set_yticks(malakia)
plt.gca().set_yticklabels(real_y)
plt.colorbar()
plt.tight_layout()
plt.savefig("Images/q0u.svg")

plt.clf()
plt.imshow(q1q, interpolation="quadric", extent=[0, 24, 0, 24])
plt.title(r"$q_c$ for $\left(q_i,u_i\right)=\left(1,0\right)$")
plt.xlabel("x(deg)")
plt.ylabel("y(deg)")
plt.gca().invert_yaxis()
plt.gca().set_xticks(malakia)
plt.gca().set_xticklabels(real_x)
plt.gca().set_yticks(malakia)
plt.gca().set_yticklabels(real_y)
plt.colorbar()
plt.tight_layout()
plt.savefig("Images/q1q.svg")

plt.clf()
plt.imshow(q1u, interpolation="quadric", extent=[0, 24, 0, 24])
plt.title(r"$u_c$ for $\left(q_i,u_i\right)=\left(1,0\right)$")
plt.xlabel("x(deg)")
plt.ylabel("y(deg)")
plt.gca().invert_yaxis()
plt.gca().set_xticks(malakia)
plt.gca().set_xticklabels(real_x)
plt.gca().set_yticks(malakia)
plt.gca().set_yticklabels(real_y)
plt.colorbar()
plt.tight_layout()
plt.savefig("Images/q1u.svg")