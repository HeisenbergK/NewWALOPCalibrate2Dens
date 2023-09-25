import numpy as np
import csv
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import sys

try:
    order = int(sys.argv[1])
except:
    order = 2

try:
    mixed = bool(int(sys.argv[2]))
except:
    mixed = True

try:
    tests = int(sys.argv[3])
except:
    tests = 1000

try:
    plim = float(sys.argv[4])
except:
    plim = 1.0

filenam = f'Results/bw_random{tests}_order{order}_mixed{mixed}_lim{plim}_modefull.csv'

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
    print(f'Processed {line_count} lines.')

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
    title=f'Partially (and zero) Polarized: order: {order}, mixed: {mixed}, teststars: {tests}, polarization limit: {plim}'
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
