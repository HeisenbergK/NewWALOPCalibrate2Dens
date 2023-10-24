import numpy as np
import csv
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import argparse
from splash import splashplot

print(chr(27) + "[2J")
splashplot()

parser = argparse.ArgumentParser(
    "Model the WALOP instrument - forward method")
parser.add_argument("order", type=int, default=2,
                    help="Fit order", choices=[2, 3])
parser.add_argument("mixed", type=int, default=1,
                    help="Whether to use the mixed term with 2nd order fit", choices=[0, 1])
parser.add_argument("tests", type=int, default=100,
                    help="Number of tests to run")
parser.add_argument("polarization_limit", type=float,
                    default=1.0, help="Polarization limit for the tests in the range [0.001, 1.0]")
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

order = args.order
mixed = bool(args.mixed)
tests = args.tests
plim = args.polarization_limit
if args.refit=="auto":
    refit="auto"
else:
    refit = int(args.refit)

if refit==1 or refit==0:
    refit=bool(refit)

noise = args.noise
if noise==0.0:
    noise = False

if plim < 0.001 or plim > 1.0:
    print("Polarization limit cannot be outside the range [0.001, 1.0]")
    exit()

filenam = f'Results/fw_random{tests}_order{order}_mixed{mixed}_lim{plim}_modefull_refit{refit}_noise{noise}.csv'

master = []
with open(filenam) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            print()
            line_count += 1
        else:
            master.append([float(row[0]), float(row[1]), float(
                row[2]), float(row[3]), float(row[4]), float(row[5]), float(row[6])])
            line_count += 1
    print(f"Plotting the following parameters:")
    print(f"Order:\t{order}")
    print(f"Mixed:\t{mixed}")
    print(f"Tests:\t{tests}")
    print(f"Polarization Limit:\t{plim}")
    print(f"Refit:\t{refit}")
    print(f'Noise:\t{noise}%')
    goon = input("Press enter to continue, Ctrl+c to abort")

unx = np.sort(np.unique(np.array(master).T[0]))
uny = np.sort(np.unique(np.array(master).T[1]))
mdq = np.empty((len(unx), len(uny)))
mdu = np.empty((len(unx), len(uny)))
sdq = np.empty((len(unx), len(uny)))
sdu = np.empty((len(unx), len(uny)))
mdq[:] = np.nan
mdu[:] = np.nan
sdq[:] = np.nan
sdu[:] = np.nan


for i in range(0, len(master)):
    if not master[i][6] > 0:
        curx = master[i][0]
        cury = master[i][1]
        indx = np.where(unx == curx)[0]
        indy = np.where(uny == cury)[0]
        mdq[indy, indx] = master[i][2]
        mdu[indy, indx] = master[i][3]
        if master[i][4] < 0.002 and master[i][5] < 0.002:
            sdq[indy, indx] = master[i][4]
            sdu[indy, indx] = master[i][5]

fig = make_subplots(
    rows=2, cols=2,
    subplot_titles=("Inverted q residual", "Inverted u residual",
                    "Inverted q std", "Inverted u std"),
    column_widths=[0.5, 0.5])

fig.add_trace(go.Contour(
    z=mdq,
    x=unx,  # horizontal axis
    y=uny,
    colorbar=dict(exponentformat='power', len=0.4, x=0.45, y=0.8),
    line_smoothing=0.95,
    contours=dict(
        size=0.01,
    )),
    row=1, col=1)

fig.add_trace(go.Contour(
    z=mdu,
    x=unx,  # horizontal axis
    y=uny,
    colorbar=dict(exponentformat='power', len=0.4, x=1, y=0.8),
    line_smoothing=0.95,
    contours=dict(
        size=0.01,
    )),
    row=1, col=2)

fig.add_trace(go.Contour(
    z=sdq,
    x=unx,  # horizontal axis
    y=uny,
    colorbar=dict(exponentformat='power', len=0.4, x=0.45, y=0.2),
    line_smoothing=0.95,
    contours=dict(
        size=0.01,
    )),
    row=2, col=1)

fig.add_trace(go.Contour(
    z=sdu,
    x=unx,  # horizontal axis
    y=uny,
    colorbar=dict(exponentformat='power', len=0.4, x=1, y=0.2),
    line_smoothing=0.95,
    contours=dict(
        size=0.01,
    )),
    row=2, col=2)

fig['layout'].update(
    title=f'Partially (and zero) Polarized: order: {order}, mixed: {mixed}, teststars: {tests}, polarization limit: {plim}, refit: {refit}, noise:{noise}%'
)
fig['layout']['yaxis1'].update(
    scaleanchor="x1",
    scaleratio=1
)
fig['layout']['yaxis2'].update(
    scaleanchor="x2",
    scaleratio=1
)
fig['layout']['yaxis3'].update(
    scaleanchor="x3",
    scaleratio=1
)
fig['layout']['yaxis4'].update(
    scaleanchor="x4",
    scaleratio=1
)
fig['layout']['xaxis1'].update(
    range=[-0.25, 0.25]
)
fig['layout']['yaxis2'].update(
    range=[-0.25, 0.25]
)
fig['layout']['yaxis3'].update(
    range=[-0.25, 0.25]
)
fig['layout']['yaxis4'].update(
    range=[-0.25, 0.25]
)
fig.show()
