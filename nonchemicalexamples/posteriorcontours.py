import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
matplotlib.use('Agg')
plt.rcParams.update({'font.size': 22})
numpoints=1000
QoI_0 = np.linspace(-4.0, 4.0, numpoints)
QoI_1 = np.linspace(-4.0, 4.0, numpoints)
QoI_0, QoI_1 = np.meshgrid(QoI_0, QoI_1)
contours = 1/np.sqrt(2*np.pi)*np.exp(-(QoI_0 - 0.0)**2)*1/np.sqrt(2*np.pi)*np.exp(-(QoI_1 - 0.0)**2)
likelihood_probability = plt.contour(QoI_0, QoI_1, contours,levels=np.linspace(0,0.2,21), cmap=plt.get_cmap('Reds'))
plt.xlabel(r'$QoI_0$')
plt.ylabel(r'$QoI_1$')
plt.ylim(-2,2)
plt.xlim(-2,2)
plt.colorbar()
plt.tight_layout()
plt.savefig('likelihood_probability.png', dpi=220)
