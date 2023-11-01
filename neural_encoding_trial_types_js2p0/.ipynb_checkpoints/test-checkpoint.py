import numpy as np
from preprocessing import GaussianSmoothing
import matplotlib.pyplot as plt


x = np.zeros((90, 100))
for i in range(90):
    x[i, i:(i + 10)] = np.ones(10)


my = GaussianSmoothing(x, sigma=5, length=20, axis='column')
conv_x = my.conv()