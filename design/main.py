import numpy as np
import matplotlib.pyplot as plt


x = np.linspace(0, 2*np.pi, 100)
y = np.sin(x)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x,y)
fig.tight_layout()
plt.show()
