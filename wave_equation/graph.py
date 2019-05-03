import glob
import os
import shutil

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D


def load(item):
    with open(item) as f:
        time = f.readline()
    matr = pd.read_csv(item, sep=' ', skiprows=1, header=None)
    return (np.array(matr), time)


def norm(item):
    X, Y = np.arange(
        0, 1, 1.0 / item.shape[1]), np.arange(0, 1, 1.0 / item.shape[0])
    X, Y = np.meshgrid(X, Y)
    return X, Y, item


def update_plot(frm, data, plot):
    plot[0].remove()
    ax.set_title(str(data[frm][1]))
    plot[0] = ax.plot_surface(*norm(data[frm][0]),  cmap=cm.coolwarm)


if __name__ == "__main__":
    FOLDERNAME = 'res/{}.txt'
    GIFPATH = 'plots/result.gif'
    SAVE = True
    DELAY = 10  # default 200 as milliseconds

    N = len(glob.glob(FOLDERNAME.format('*')))
    data = list(map(load, [FOLDERNAME.format(x) for x in range(N)]))

    fig=plt.figure()
    ax=fig.add_subplot(111, projection='3d')
    plot=[ax.plot_surface(*norm(data[0][0]), cmap=cm.Blues)]
    ax.set_zlim(0, 0.5)
    gif=animation.FuncAnimation(fig, update_plot, len(
        data), fargs=(data, plot), interval=DELAY)

    plt.show()

    if SAVE:
        if os.path.exists("plots"):
            shutil.rmtree("plots")
        os.mkdir("plots")

        gif.save(GIFPATH, writer='pillow', fps=len(data))
