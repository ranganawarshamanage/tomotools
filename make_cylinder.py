# create a cylinder and output as a mrc file

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def create_cylinder_array(radius, height, length, resolution=100):
    # Set up a grid
    x = np.linspace(-length/2, length/2, resolution)
    y = np.linspace(-radius, radius, resolution)
    z = np.linspace(0, height, resolution)
    x, y, z = np.meshgrid(x, y, z)

    # Create a 3D array representing the cylinder
    cylinder_array = np.zeros_like(x, dtype=int)
    cylinder_mask = x**2 + y**2 <= radius**2
    cylinder_array[cylinder_mask] = 1

    return cylinder_array

# Example usage with a cylinder of radius 0.2, height 1.0, and length 1.0
radius = 0.2
height = 1.0
length = 1.0
cylinder_array = create_cylinder_array(radius, height, length)

# Display the numpy array
print(cylinder_array)

# Plot the cylinder using matplotlib for visualization
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.voxels(cylinder_array, edgecolor='k')

ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_zlabel('Z-axis')
ax.set_title('Cylinder as Numpy Array')

plt.show()

