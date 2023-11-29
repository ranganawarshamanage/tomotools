import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import starfile
import pandas as pd
import argparse


parser = argparse.ArgumentParser(description='')
parser.add_argument('--particlefile1', type=str, required=True, help='Input particles.star file')
parser.add_argument('--particlefile2', type=str, required=True, help='Input particles.star file')
parser.add_argument('--distance_cutoff', type=float, required=True, help='Specify the distance cutoff in Angstroms')
parser.add_argument('--pixsize', type=float, required=True, help='Specify the pixel size of the tilt series')
args = parser.parse_args()


def euler_rotmat(angles):
    phi, theta, psi = angles
    R_x = np.array([[1, 0, 0],
                    [0, np.cos(phi), -np.sin(phi)],
                    [0, np.sin(phi), np.cos(phi)]])
   
    R_y = np.array([[np.cos(theta), 0, np.sin(theta)],
                    [0, 1, 0],
                    [-np.sin(theta), 0, np.cos(theta)]])
   
    R_z = np.array([[np.cos(psi), -np.sin(psi), 0],
                    [np.sin(psi), np.cos(psi), 0],
                    [0, 0, 1]])
   
    return R_z @ R_y @ R_x



distance_cutoff = args.distance_cutoff / args.pixsize

bluedf = starfile.read(args.particlefile1)
bluelist = bluedf["rlnClassNumber"]
xyz_blues = []
ang_blues = []
for i, _ in enumerate(bluelist):
    xyz_blues.append(
        [bluedf["rlnCoordinateX"][i], bluedf["rlnCoordinateY"][i], bluedf["rlnCoordinateZ"][i]])
    ang_blues.append(
        [bluedf["rlnAngleRot"][i], bluedf["rlnAngleTilt"][i], bluedf["rlnAnglePsi"][i]])

reddf = starfile.read(args.particlefile2)
redlist = reddf["rlnClassNumber"]
xyz_reds = []
for i, _ in enumerate(redlist):
    xyz_reds.append([reddf["rlnCoordinateX"][i], reddf["rlnCoordinateY"][i], reddf["rlnCoordinateZ"][i]])

blues = np.array(xyz_blues)
reds = np.array(xyz_reds)
angs = np.array(ang_blues)

# Create a figure and a 3D axis
fig = plt.figure()

selected_vecs = []
selected_angles = []

for r in reds:
    vecs = blues - r
    magnitudes = np.linalg.norm(vecs, axis=1)
    index_of_min_magnitude = np.argmin(magnitudes)
    shortest_vec = vecs[index_of_min_magnitude]
    shortest_vec_ang = angs[index_of_min_magnitude]
    if np.amin(magnitudes) <= distance_cutoff:
        selected_vecs.append(shortest_vec)
        selected_angles.append(shortest_vec_ang)

selected_vecs = np.array(selected_vecs)
selected_angles = np.deg2rad(np.array(selected_angles))

# sanity check
# selected_vecs = reds
# selected_angles = angs


# ================
# First plot- vecs before rotation 
# ================
ax = fig.add_subplot(1,2,1, projection='3d')
for v, angle in zip(selected_vecs, selected_angles):
    norm_v = v / np.linalg.norm(v)
    ax.plot([0, norm_v[0]], [0, norm_v[1]], [0, norm_v[2]], color='r')

    # Set plot limits
    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.set_zlim([-1, 1])

# ================
# Second plot- vecs after rotation 
# ================
ax = fig.add_subplot(1,2,2, projection='3d')
for v, angle in zip(selected_vecs, selected_angles):
    rm = euler_rotmat(angle)
    v = rm @ v
    norm_v = v / np.linalg.norm(v)
    ax.plot([0, norm_v[0]], [0, norm_v[1]], [0, norm_v[2]], color='k')

# Set axis labels
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_zlabel('Z-axis')

# Set plot limits
ax.set_xlim([-1, 1])
ax.set_ylim([-1, 1])
ax.set_zlim([-1, 1])

# Add a legend
ax.legend()

# Show the plot
plt.show()