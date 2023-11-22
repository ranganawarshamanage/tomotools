import numpy as np
import starfile
import pandas as pd
import argparse
import traceback

parser = argparse.ArgumentParser(description='')
parser.add_argument('--particles', type=str, required=True, help='Input particles.star file')
#parser.add_argument('--tomoname', type=str, required=True, help='Specify tomogram name')
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

    # Combine indices to find c points near both a and b
    #indices_near_both_a_b = np.logical_or(indices_a_within_cutoff, indices_b_within_cutoff)

    # Extract c points near both a and b
    # c_points_near_a = c_points[indices_a_within_cutoff]
    # c_points_near_b = c_points[indices_b_within_cutoff]
    # c_points_near_a_b = np.concatenate((c_points_near_a, c_points_near_b), axis=0)

    # return points_a_within_cutoff, points_b_within_cutoff#, c_points_near_a_b
    return indices_a_within_cutoff, indices_b_within_cutoff


def select_by_tomoname(df, tname):
    classlist = df["rlnClassNumber"]
    xyz_cls1 = []
    xyz_cls2 = []
    xyz_cls4 = []
    rtp_cls1 = []
    rtp_cls2 = []
    rtp_cls4 = []

    for i, _ in enumerate(classlist):
        if df["rlnTomoName"][i] == tname:
            if df["rlnClassNumber"][i] == 1:
                xyz_cls1.append([df["rlnCoordinateX"][i], df["rlnCoordinateY"][i], df["rlnCoordinateZ"][i]])
                rtp_cls1.append([df["rlnAngleRot"][i], df["rlnAngleTilt"][i], df["rlnAnglePsi"][i]])
            if df["rlnClassNumber"][i] == 2:
                xyz_cls2.append([df["rlnCoordinateX"][i], df["rlnCoordinateY"][i], df["rlnCoordinateZ"][i]])
                rtp_cls2.append([df["rlnAngleRot"][i], df["rlnAngleTilt"][i], df["rlnAnglePsi"][i]])
            if df["rlnClassNumber"][i] == 4:
                xyz_cls4.append([df["rlnCoordinateX"][i], df["rlnCoordinateY"][i], df["rlnCoordinateZ"][i]])
                rtp_cls4.append([df["rlnAngleRot"][i], df["rlnAngleTilt"][i], df["rlnAnglePsi"][i]])

    # output class 2 points for sanity check

    dfx = pd.DataFrame()
    xx = [row[0] for row in xyz_cls2]
    dfx["rlnTomoName"] = [tname for _ in range(len(xx))]
    dfx["rlnCoordinateX"] = [row[0] for row in xyz_cls2]
    dfx["rlnCoordinateY"] = [row[1] for row in xyz_cls2]
    dfx["rlnCoordinateZ"] = [row[2] for row in xyz_cls2]
    dfx["rlnClassNumber"] = [2 for _ in range(len(xx))]

    outfilename = 'cls2_particles_%s.star' % tname
    starfile.write(dfx, outfilename, overwrite=True)

    distance_cutoff = 500. # Angstroms

    a_points = np.array(xyz_cls1, dtype=float)
    b_points = np.array(xyz_cls4, dtype=float)
    c_points = np.array(xyz_cls2, dtype=float)

    a_points_rtp = np.array(rtp_cls1, dtype=float)
    b_points_rtp = np.array(rtp_cls4, dtype=float)


    indices_a_within_cutoff, indices_b_within_cutoff = find_points_within_distance(
        a_points=a_points, 
        b_points=b_points, 
        c_points=c_points, 
        distance_cutoff=distance_cutoff)

    points_a_within_cutoff = a_points[indices_a_within_cutoff]
    points_b_within_cutoff = b_points[indices_b_within_cutoff]

    a_points_rtp_within_cutoff = a_points_rtp[indices_a_within_cutoff]
    b_points_rtp_within_cutoff = b_points_rtp[indices_b_within_cutoff]



    tma = []
    if 0 < len(points_a_within_cutoff): 
        xa = points_a_within_cutoff[:, 0]
        ya = points_a_within_cutoff[:, 1]
        za = points_a_within_cutoff[:, 2]
        ra = a_points_rtp_within_cutoff[:, 0]
        ta = a_points_rtp_within_cutoff[:, 1]
        pa = a_points_rtp_within_cutoff[:, 2]
        tma = [tname for _ in range(len(xa))]
        cls1_number = [1 for _ in range(len(xa))]

    tmb = []
    if 0 < len(points_b_within_cutoff):
        xb = points_b_within_cutoff[:, 0]
        yb = points_b_within_cutoff[:, 1]
        zb = points_b_within_cutoff[:, 2]
        rb = b_points_rtp_within_cutoff[:, 0]
        tb = b_points_rtp_within_cutoff[:, 1]
        pb = b_points_rtp_within_cutoff[:, 2]
        tmb = [tname for _ in range(len(xb))]
        cls4_number = [4 for _ in range(len(xb))]

    # Output those particles into a star file
    try:
        if 0 < len(tma) and 0 < len(tmb): 
            tomoname = tma + tmb
            cls_number = cls1_number + cls4_number
            x = xa.tolist() + xb.tolist()
            y = ya.tolist() + yb.tolist()
            z = za.tolist() + zb.tolist()
            rot = ra.tolist() + rb.tolist()
            tilt = ta.tolist() + tb.tolist()
            psi = pa.tolist() + pb.tolist()
        elif 0 < len(tma) and not (0 < len(tmb)):
            tomoname = tma
            cls_number = cls1_number
            x = xa.tolist()
            y = ya.tolist()
            z = za.tolist()
            rot = ra.tolist() 
            tilt = ta.tolist() 
            psi = pa.tolist() 
        elif not (0 < len(tma)) and 0 < len(tmb):
            tomoname = tmb
            cls_number = cls4_number
            x = xb.tolist()
            y = yb.tolist()
            z = zb.tolist()
            rot = rb.tolist()
            tilt = tb.tolist()
            psi = pb.tolist()
    except:
        print(traceback.exit())

    return x, y, z, rot, tilt, psi, cls_number, tomoname


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

uniq_tnames = set(df["rlnTomoName"])
xlist = []
ylist = []
zlist = []
rotlist = []
tiltlist = []
psilist = []
clslist = []
tomolist = []
for tname in uniq_tnames:
    x, y, z, rot, tilt, psi, cls_number, tomoname = select_by_tomoname(df, tname)

    df2 = pd.DataFrame()
    df2["rlnTomoName"] = tomoname
    df2["rlnCoordinateX"] = x
    df2["rlnCoordinateY"] = y
    df2["rlnCoordinateZ"] = z
    df2["rlnAngleRot"] = rot
    df2["rlnAngleTilt"] = tilt
    df2["rlnAnglePsi"] = psi
    df2["rlnClassNumber"] = cls_number

    outfilename = 'neighbor_particles_%s.star' % tname
    starfile.write(df2, outfilename, overwrite=True)

    xlist += x
    ylist += y
    zlist += z
    rotlist += rot
    tiltlist += tilt
    psilist += psi
    clslist += cls_number
    tomolist += tomoname

# Output all particles
df3 = pd.DataFrame()
df3["rlnTomoName"] = tomolist
df3["rlnCoordinateX"] = xlist
df3["rlnCoordinateY"] = ylist
df3["rlnCoordinateZ"] = zlist
df3["rlnAngleRot"] = rotlist
df3["rlnAngleTilt"] = tiltlist
df3["rlnAnglePsi"] = psilist
df3["rlnClassNumber"] = clslist

outfilename = 'neighbor_particles.star'
starfile.write(df3, outfilename, overwrite=True)