import numpy as np
import starfile
import pandas as pd
import argparse
import traceback

parser = argparse.ArgumentParser(description='')
parser.add_argument('--particles', type=str, required=True, help='Input particles.star file')
parser.add_argument('--tomoname', type=str, required=True, help='Specify tomogram name')
#parser.add_argument('--binfactor', type=float, required=True, help='Specify the bin factor')
args = parser.parse_args()


def find_points_within_distance(a_points, b_points, c_points, distance_cutoff):
    # Convert point arrays to numpy arrays
    a_points = np.array(a_points)
    b_points = np.array(b_points)
    c_points = np.array(c_points)

    # Calculate distances between each point of type c and points of type a
    distances_a = np.linalg.norm(a_points[:, np.newaxis, :] - c_points, axis=2)

    """ for i in range(2):
        for j in range(3):
            v = a_points[j] - c_points[i]
            print(np.sqrt(np.sum(v*v))) """

    # Find indices where distances are within the cutoff for points of type a
    indices_a_within_cutoff = np.any(distances_a < distance_cutoff, axis=1)

    # Calculate distances between each point of type c and points of type b
    distances_b = np.linalg.norm(b_points[:, np.newaxis, :] - c_points, axis=2)

    # Find indices where distances are within the cutoff for points of type b
    indices_b_within_cutoff = np.any(distances_b < distance_cutoff, axis=1)

    # Extract points of type a and b within the cutoff
    points_a_within_cutoff = a_points[indices_a_within_cutoff]
    points_b_within_cutoff = b_points[indices_b_within_cutoff]

    return points_a_within_cutoff, points_b_within_cutoff


# Example usage:
#a_points = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
#b_points = np.array([[10, 11, 12], [13, 14, 15], [16, 17, 18]])
#c_points = np.array([[19, 20, 21], [22, 23, 24]])

#distance_cutoff = 20.0

#points_a_within_cutoff, points_b_within_cutoff = find_points_within_distance(a_points, b_points, c_points, distance_cutoff)

#print("Points of type a within cutoff:", points_a_within_cutoff)
#print("Points of type b within cutoff:", points_b_within_cutoff)


# Reading from particles.star file into classes
df = starfile.read(args.particles)
classlist = df["rlnClassNumber"]
xyz_cls1 = []
xyz_cls2 = []
xyz_cls4 = []
cls_number = []
for i, _ in enumerate(classlist):
    if df["rlnTomoName"][i] == args.tomoname:
        if df["rlnClassNumber"][i] == 1:
            cls_number.append(1)
            xyz_cls1.append([df["rlnCoordinateX"][i], df["rlnCoordinateY"][i], df["rlnCoordinateZ"][i]])
        if df["rlnClassNumber"][i] == 2:
            cls_number.append(2)
            xyz_cls2.append([df["rlnCoordinateX"][i], df["rlnCoordinateY"][i], df["rlnCoordinateZ"][i]])
        if df["rlnClassNumber"][i] == 4:
            cls_number.append(4)
            xyz_cls4.append([df["rlnCoordinateX"][i], df["rlnCoordinateY"][i], df["rlnCoordinateZ"][i]])

distance_cutoff = 1000. # Angstroms

points_a_within_cutoff, points_b_within_cutoff = find_points_within_distance(
    a_points=np.array(xyz_cls1, dtype=float), 
    b_points=np.array(xyz_cls4, dtype=float), 
    c_points=np.array(xyz_cls2, dtype=float), 
    distance_cutoff=distance_cutoff)

tma = []
if 0 < len(points_a_within_cutoff): 
    xa = points_a_within_cutoff[:, 0]
    ya = points_a_within_cutoff[:, 1]
    za = points_a_within_cutoff[:, 2]
    tma = [args.tomoname for _ in range(len(xa))]

tmb = []
if 0 < len(points_b_within_cutoff):
    xb = points_b_within_cutoff[:, 0]
    yb = points_b_within_cutoff[:, 1]
    zb = points_b_within_cutoff[:, 2]
    tmb = [args.tomoname for _ in range(len(xb))]

# Output those particles into a star file
try:
    if 0 < len(tma) and 0 < len(tmb): 
        tomoname = tma + tmb
        x = xa.tolist() + xb.tolist()
        y = ya.tolist() + yb.tolist()
        z = za.tolist() + zb.tolist()
    elif 0 < len(tma) and not (0 < len(tmb)):
        tomoname = tma
        x = xa.tolist()
        y = ya.tolist()
        z = za.tolist()
    elif not (0 < len(tma)) and 0 < len(tmb):
        tomoname = tmb
        x = xb.tolist()
        y = yb.tolist()
        z = zb.tolist()
except:
    print(traceback.exit())

df2 = pd.DataFrame()
df2["rlnTomoName"] = tomoname
df2["rlnCoordinateX"] = x
df2["rlnCoordinateY"] = y
df2["rlnCoordinateZ"] = z
df2["rlnClassNumber"] = cls_number

outfilename = 'neighbor_particles_%s.star' %args.tomoname
starfile.write(df2, outfilename, overwrite=True)