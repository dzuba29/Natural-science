import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import argparse


def surface(data, file):
    Z = np.array(data)
    X = np.arange(0, 1, 1.0 / Z.shape[1])
    Y = np.arange(0, 1, 1.0 / Z.shape[0])
    X2D, Y2D = np.meshgrid(X, Y)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot_surface(X2D, Y2D, Z)
    plt.subplots_adjust(wspace=0.5, hspace=0.6)
    plt.savefig('{}.pdf'.format(file))


data = pd.read_csv("result.txt", header=None, sep=' ')
surface(data, "surface")
