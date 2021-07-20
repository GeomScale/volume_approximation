import numpy as np
import matplotlib.pyplot as plt 

def get_ellipsoid_vals(A, c, x1, x2):
    X = np.vstack([x1, x2]).T
    vals = [(x-c).T @ A @ (x-c) for x in X]
    return np.array(vals)
 

# plot data
filename = "ellipsoid.txt"
data = np.genfromtxt(filename, delimiter=' ')
plt.plot(data[:, 0], data[:, 1], 'x')

# plot ellipsoid
# Plot (x-c)' A (x-c) = 1
# A = L @ L.T
c = np.array([2, 2])
L = np.array([[0.5, 0.0], [1.5, 1.0]])
A = L @ (L.T)

tol = 0.5
x1_vec, x2_vec = np.meshgrid(np.linspace(data[:, 0].min() - tol, data[:, 0].max() + tol, 100), np.linspace(data[:, 1].min() - tol, data[:, 1].max() + tol, 100))
vals = get_ellipsoid_vals(A, c, x1_vec.ravel(), x2_vec.ravel())
vals = vals.reshape(x1_vec.shape)
plt.contour(x1_vec, x2_vec, vals, levels=[1], cmap="Greys_r")

# show all plots
plt.show()