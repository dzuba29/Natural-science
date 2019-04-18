# T = 1.0
# nT = 10000.0
# def test(t, h1, h2):
#     d = t/pow(h1, 2) + t/pow(h2, 2)
#     return d <= 0.5


# for x in range(10, 1000, 10):
#     if test(T/nT, 1./x, 1./x):
#         print("!!!", x)

import glob
import os
import shutil

import imageio
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D


def surface(data, file, time):
    Z = np.array(data)
    X = np.arange(0, 1, 1.0 / Z.shape[1])
    Y = np.arange(0, 1, 1.0 / Z.shape[0])
    X2D, Y2D = np.meshgrid(X, Y)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot_surface(X2D, Y2D, Z)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('t')
    ax.set_zlim(0, 1)
    plt.title(str(time), pad=0.5)
    plt.subplots_adjust(wspace=0.5, hspace=0.6)
    plt.savefig('{}.png'.format(file))
    plt.clf()
    plt.close()

    
n = len(glob.glob("plots/*.png"))
images = []
for i in range(0, n):
	images.append(imageio.imread('plots/'+str(i)+'.png'))
imageio.mimsave('plots/movie.gif', images, duration=0.5)
