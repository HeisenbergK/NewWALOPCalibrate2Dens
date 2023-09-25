import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import sys
from astropy.io import ascii
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

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

filenam = f'Results/fw_random{tests}_order{order}_mixed{mixed}_lim{plim}_modefull_covs.csv'
limit = 1

master = ascii.read(filenam)
covariances = []
for i in range(0, len(master)):
    temp = [float(j) for j in master['covariance'][i].replace('[', '').replace(']', '').replace('\n', '').split()]
    covariances.append(temp)
master['covariance'] = covariances
# print(master.pprint_all())

unx = np.sort(np.unique(np.array(master['x'])))
uny = np.sort(np.unique(np.array(master['y'])))
cov11 = np.empty((len(unx), len(uny)))
cov12 = np.empty((len(unx), len(uny)))
cov21 = np.empty((len(unx), len(uny)))
cov22 = np.empty((len(unx), len(uny)))
cov11[:] = np.nan
cov12[:] = np.nan
cov21[:] = np.nan
cov22[:] = np.nan

for i in range(0, len(master)):
    curx = float(master['x'][i])
    cury = float(master['y'][i])
    indx = np.where(unx == curx)[0][0]
    indy = np.where(uny == cury)[0][0]
    if master['covariance'][i][0] <= limit and master['covariance'][i][3] <= limit:
        cov11[indy, indx] = float(master['covariance'][i][0])
        cov12[indy, indx] = float(master['covariance'][i][1])
        cov21[indy, indx] = float(master['covariance'][i][2])
        cov22[indy, indx] = float(master['covariance'][i][3])

# plt.clf()
# fig = plt.figure(figsize=(8, 8))
# gs = mpl.gridspec.GridSpec(2, 2, height_ratios=[2, 2], width_ratios=[2, 2])
# ax1 = fig.add_subplot(gs[0, 0])
# ax2 = fig.add_subplot(gs[1, 0])
# ax3 = fig.add_subplot(gs[0, 1])
# ax4 = fig.add_subplot(gs[1, 1])
# cont1 = ax1.contourf(unx, uny, cov11, 20)
# cont2 = ax2.contourf(unx, uny, cov21, 20)
# cont3 = ax3.contourf(unx, uny, cov12, 20)
# cont4 = ax4.contourf(unx, uny, cov22, 20)
# plt.tick_params(which='both', top=False, right=False)
# plt.show()


fig = make_subplots(
    rows=2, cols=2,
    subplot_titles=("Cqq", "Cqu",
                    "Cuq", "Cuu"),
    column_widths=[0.5, 0.5])

fig.add_trace(go.Contour(
    z=cov11,
    x=unx,  # horizontal axis
    y=uny,
    colorbar=dict(exponentformat='power', len=0.4, x=0.45, y=0.8),
    line_smoothing=0.95,
    contours=dict(
        size=0.01,
    )),
    row=1, col=1)

fig.add_trace(go.Contour(
    z=cov12,
    x=unx,  # horizontal axis
    y=uny,
    colorbar=dict(exponentformat='power', len=0.4, x=1, y=0.8),
    line_smoothing=0.95,
    contours=dict(
        size=0.01,
    )),
    row=1, col=2)

fig.add_trace(go.Contour(
    z=cov21,
    x=unx,  # horizontal axis
    y=uny,
    colorbar=dict(exponentformat='power', len=0.4, x=0.45, y=0.2),
    line_smoothing=0.95,
    contours=dict(
        size=0.01,
    )),
    row=2, col=1)

fig.add_trace(go.Contour(
    z=cov22,
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