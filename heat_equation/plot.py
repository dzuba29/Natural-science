import glob
import math
import os

import imageio
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import argparse


def surface(data, file, time):
    Z = np.array(data)
    X = np.arange(0, 1, 1.0 / Z.shape[1])
    Y = np.arange(0, 1, 1.0 / Z.shape[0])
    X2D, Y2D = np.meshgrid(X, Y)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot_surface(X2D, Y2D, Z, cmap=cm.coolwarm)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('t')
    ax.set_zlim(0, 1)
    plt.title(str(time), pad=0.5)
    plt.subplots_adjust(wspace=0.5, hspace=0.6)
    plt.savefig('{}.png'.format(file))
    plt.clf()
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--n', type=int, default=100)
    args = parser.parse_args()

    res = glob.glob('res/*')
    print(len(res))

    if (len(res) == 0):
        raise NotImplementedError

    for i in range(0, args.n):
        p = i/(args.n-1) * 100
        print(f"{p}%")
        with open(res[i]) as f:
            time = f.readline()
            data = pd.read_csv(res[i], sep=' ', skiprows=1, header=None)
            surface(data, 'plots/'+str(i), time)

    n = len(glob.glob("plots/*.png"))
    images = []
    for idx,item in enumerate(sorted(n)):
        images.append(imageio.imread('plots/'+str(idx)+'.png'))
    imageio.mimsave('plots/movie.gif', images, duration=0.5)
