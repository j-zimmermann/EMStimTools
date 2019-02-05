import numpy as np
import matplotlib.pyplot as plt


def Tupaj_func(V, a, d, x):
    return V / (2. * np.log(a / (d - a))) * (((x - d / 2.) / ((x - d / 2.)**2)) - (((x + d / 2.)) / ((x + d / 2.)**2)))


# voltage
V = 1.0
# radius of rod
a = 0.0005
# spacing between electrodes
d = 0.022
xmin = -0.01025
xmax = 0.01025

xvalues = np.linspace(xmin, xmax, 50)
result = np.empty(len(xvalues))
for i in range(len(xvalues)):
    result[i] = Tupaj_func(V, a, d, xvalues[i])


# shift for distance
shift = xmin + 0.5 * ((xmax - xmin) - d)
print(xmin)
print(xmin + ((xmax - xmin) - d))
print(shift)
xvalues = xvalues - shift
# plot and output
plt.plot(xvalues, result)
plt.show()
np.savetxt('field_tupaj.txt', np.transpose([xvalues, result]), delimiter='\t', header='position\tAnalytical')
