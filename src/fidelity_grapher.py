import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

while True:
	data = pd.read_csv('QHA_output.txt', names=["index", "value"])

	X = np.array(data.index)
	y = np.array(data.value)

	y = pd.rolling_mean(data.value, window = 50)

	plt.ion()
	plt.xlabel('Number of iterations')
	plt.ylabel('Likelihood Fidelity')
	plt.title('SGD likelihood for QHA')
	plt.plot(X, y, color='r')
	
	plt.pause(2)
