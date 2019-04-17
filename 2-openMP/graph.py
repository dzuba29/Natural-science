import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
import imageio

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
	ax.set_zlim(0, 2)
	plt.title(str(time), pad=0.5)
	plt.subplots_adjust(wspace=0.5, hspace=0.6)
	plt.savefig('{}.png'.format(file))

n = 10
images = []
for i in range(0, n+1):
	with open('res/'+str(i)+'.txt') as f:
		time = f.readline()
		data = pd.read_csv('res/'+str(i)+'.txt', sep = ' ', skiprows=1, header = None)
		surface(data,'plots/'+str(i), time)
	images.append(imageio.imread('plots/'+str(i)+'.png'))

imageio.mimsave('plots/movie.gif', images, duration=0.3)



# png_dir = 'plots/'

# for file_name in os.listdir(png_dir):
#     if file_name.endswith('.png'):
#         file_path = os.path.join(png_dir, file_name)
#         images.append(imageio.imread(file_path))
# imageio.mimsave('plots/movie.gif', images, duration=0.3)