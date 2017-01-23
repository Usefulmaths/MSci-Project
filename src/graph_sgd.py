import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

data = pd.read_csv('sgd_toffoli.txt')

X = data['index']
y = data['value']

while True:

	data = pd.read_csv('sgd_toffoli2.txt')

	X = np.array(data['index'])
	y = np.array(data['value'])

	plt.ion()
	plt.plot(X, y, color='r')
	plt.pause(2)