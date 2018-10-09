#!/usr/bin/env pythob
import numpy as np
import matplotlib.pyplot as plt

theta = np.arange(0., 2.*np.pi,0.01)

fig = plt.figure(1)
plt.plot(theta, 1.*np.cos(theta) - 1.*np.sin(theta))
plt.show()


