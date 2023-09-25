import numpy as np
import csv
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from splash import splash
import argparse

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
    choices=["0","1","2","3","4","5","6","7","8","9","10","auto"],
    metavar="[0-10/auto]",
)
parser.add_argument(
    "noise",
    type=float,
    default=0.0,
    help="Noise for the model in %%",
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
goon = input("Press enter to continue, Ctrl+c to abort")

filenam = f"Results/fw_random{_tests}_order{_order}_mixed{_mixed}_lim{_plim}_modemodel_refit{_refit}_noise{_noise}.csv"

master = []
with open(filenam) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            print()
            line_count += 1
        else:
            master.append([float(row[0]), float(row[1]),
                           [float(i) for i in row[2].split('[')[1].split(']')[0].split(', ')],
                           [float(i) for i in row[3].split('[')[1].split(']')[0].split(', ')]])
            line_count += 1
    print(f'Processed {line_count} lines.')

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

for i in range(0, len(master)):
    curx = master[i][0]
    cury = master[i][1]
    indx = np.where(unx == curx)[0]
    indy = np.where(uny == cury)[0]
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

        popt_q_q3[indy, indx] = master[i][2][6]
        popt_u_q3[indy, indx] = master[i][3][6]

        popt_q_u3[indy, indx] = master[i][2][7]
        popt_u_u3[indy, indx] = master[i][3][7]

        popt_q_q2u[indy, indx] = master[i][2][8]
        popt_u_q2u[indy, indx] = master[i][3][8]

        popt_q_qu2[indy, indx] = master[i][2][9]
        popt_u_qu2[indy, indx] = master[i][3][9]
    except:
        pass

# plt.clf()
# fig = plt.figure(figsize=(8, 8))
# gs = mpl.gridspec.GridSpec(2, 2, height_ratios=[2, 2], width_ratios=[2, 2])
# ax1 = fig.add_subplot(gs[0, 0])
# ax2 = fig.add_subplot(gs[1, 0])
# ax3 = fig.add_subplot(gs[0, 1])
# ax4 = fig.add_subplot(gs[1, 1])
# cont1 = ax1.contourf(unx, uny, popt_q_0, 20)
# cont2 = ax2.contourf(unx, uny, popt_u_0, 20)
# cont3 = ax3.contourf(unx, uny, popt_q_q, 20)
# cont4 = ax4.contourf(unx, uny, popt_u_q, 20)
# plt.tick_params(which='both', top=False, right=False)
# plt.show()

mean_popt_q_0 = round(np.median(popt_q_0),3)
mean_popt_q_q = round(np.median(popt_q_q),3)
mean_popt_q_u = round(np.median(popt_q_u),3)
mean_popt_q_q2 = round(np.median(popt_q_q2),3)
mean_popt_q_u2 = round(np.median(popt_q_u2),3)
mean_popt_q_qu = round(np.median(popt_q_qu),3)
mean_popt_q_q3 = round(np.median(popt_q_q3),3)
mean_popt_q_u3 = round(np.median(popt_q_u3),3)
mean_popt_q_q2u = round(np.median(popt_q_q2u),3)
mean_popt_q_qu2 = round(np.median(popt_q_qu2),3)
mean_popt_u_0 = round(np.median(popt_u_0),3)
mean_popt_u_q = round(np.median(popt_u_q),3)
mean_popt_u_u = round(np.median(popt_u_u),3)
mean_popt_u_q2 = round(np.median(popt_u_q2),3)
mean_popt_u_u2 = round(np.median(popt_u_u2),3)
mean_popt_u_qu = round(np.median(popt_u_qu),3)
mean_popt_u_q3 = round(np.median(popt_u_q3),3)
mean_popt_u_u3 = round(np.median(popt_u_u3),3)
mean_popt_u_q2u = round(np.median(popt_u_q2u),3)
mean_popt_u_qu2 = round(np.median(popt_u_qu2),3)

std_popt_q_0 = round(np.std(popt_q_0),3)
std_popt_q_q = round(np.std(popt_q_q),3)
std_popt_q_u = round(np.std(popt_q_u),3)
std_popt_q_q2 = round(np.std(popt_q_q2),3)
std_popt_q_u2 = round(np.std(popt_q_u2),3)
std_popt_q_qu = round(np.std(popt_q_qu),3)
std_popt_q_q3 = round(np.std(popt_q_q3),3)
std_popt_q_u3 = round(np.std(popt_q_u3),3)
std_popt_q_q2u = round(np.std(popt_q_q2u),3)
std_popt_q_qu2 = round(np.std(popt_q_qu2),3)
std_popt_u_0 = round(np.std(popt_u_0),3)
std_popt_u_q = round(np.std(popt_u_q),3)
std_popt_u_u = round(np.std(popt_u_u),3)
std_popt_u_q2 = round(np.std(popt_u_q2),3)
std_popt_u_u2 = round(np.std(popt_u_u2),3)
std_popt_u_qu = round(np.std(popt_u_qu),3)
std_popt_u_q3 = round(np.std(popt_u_q3),3)
std_popt_u_u3 = round(np.std(popt_u_u3),3)
std_popt_u_q2u = round(np.std(popt_u_q2u),3)
std_popt_u_qu2 = round(np.std(popt_u_qu2),3)



fig = make_subplots(
    rows=2, cols=6,
    subplot_titles=(f"q dep 0: {mean_popt_q_0}/{std_popt_q_0}", f"q dep q: {mean_popt_q_q}/{std_popt_q_q}", 
                    f"q dep u: {mean_popt_q_u}/{std_popt_q_u}", f"q dep q^2: {mean_popt_q_q2}/{std_popt_q_q2}", 
                    f"q dep u^2: {mean_popt_q_u2}/{std_popt_q_u2}", f"q dep qu: {mean_popt_q_qu}/{std_popt_q_qu}",
                    f"u dep 0: {mean_popt_u_0}/{std_popt_u_0}", f"u dep q: {mean_popt_u_q}/{std_popt_u_q}", 
                    f"u dep u: {mean_popt_u_u}/{std_popt_u_u}", f"u dep q^2: {mean_popt_u_q2}/{std_popt_u_q2}", 
                    f"u dep u^2: {mean_popt_u_u2}/{std_popt_u_u2}", f"u dep qu: {mean_popt_u_qu}/{std_popt_u_qu}"),
    column_widths=[0.5, 0.5, 0.5, 0.5, 0.5, 0.5])

fig.add_trace(go.Contour(
    z=popt_q_0,
    x=unx,  # horizontal axis
    y=uny,
    line_smoothing=0.95,
    contours=dict(
        size=0.01,
    )),
    row=1, col=1)

fig.add_trace(go.Contour(
    z=popt_q_q,
    x=unx,  # horizontal axis
    y=uny,
    line_smoothing=0.95,
    contours=dict(
        size=0.01,
    )),
    row=1, col=2)

fig.add_trace(go.Contour(
    z=popt_q_u,
    x=unx,  # horizontal axis
    y=uny,
    line_smoothing=0.95,
    contours=dict(
        size=0.01,
    )),
    row=1, col=3)

fig.add_trace(go.Contour(
    z=popt_q_q2,
    x=unx,  # horizontal axis
    y=uny,
    line_smoothing=0.95,
    contours=dict(
        size=0.01,
    )),
    row=1, col=4)

fig.add_trace(go.Contour(
    z=popt_q_u2,
    x=unx,  # horizontal axis
    y=uny,
    line_smoothing=0.95,
    contours=dict(
        size=0.01,
    )),
    row=1, col=5)

fig.add_trace(go.Contour(
    z=popt_q_qu,
    x=unx,  # horizontal axis
    y=uny,
    line_smoothing=0.95,
    contours=dict(
        size=0.01,
    )),
    row=1, col=6)

fig.add_trace(go.Contour(
    z=popt_u_0,
    x=unx,  # horizontal axis
    y=uny,
    line_smoothing=0.95,
    contours=dict(
        size=0.01,
    )),
    row=2, col=1)

fig.add_trace(go.Contour(
    z=popt_u_q,
    x=unx,  # horizontal axis
    y=uny,
    line_smoothing=0.95,
    contours=dict(
        size=0.01,
    )),
    row=2, col=2)

fig.add_trace(go.Contour(
    z=popt_u_u,
    x=unx,  # horizontal axis
    y=uny,
    line_smoothing=0.95,
    contours=dict(
        size=0.01,
    )),
    row=2, col=3)

fig.add_trace(go.Contour(
    z=popt_u_q2,
    x=unx,  # horizontal axis
    y=uny,
    line_smoothing=0.95,
    contours=dict(
        size=0.01,
    )),
    row=2, col=4)

fig.add_trace(go.Contour(
    z=popt_u_u2,
    x=unx,  # horizontal axis
    y=uny,
    line_smoothing=0.95,
    contours=dict(
        size=0.01,
    )),
    row=2, col=5)

fig.add_trace(go.Contour(
    z=popt_u_qu,
    x=unx,  # horizontal axis
    y=uny,
    line_smoothing=0.95,
    contours=dict(
        size=0.01,
    )),
    row=2, col=6)

fig.show()
