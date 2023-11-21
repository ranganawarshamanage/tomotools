import starfile
import pandas as pd
import argparse
import napari

particles1 =     '/cephfs/ranganaw/processing/G3_all/cls2_particles_Grid_G3_T4.star'
particles4 = '/cephfs/ranganaw/processing/G3_all/neighbor_particles_Grid_G3_T4.star'
df1 = starfile.read(particles1)
df4 = starfile.read(particles4)

x1 = df1["rlnCoordinateX"]
y1 = df1["rlnCoordinateY"]
z1 = df1["rlnCoordinateZ"]

x4 = df4["rlnCoordinateX"]
y4 = df4["rlnCoordinateY"]
z4 = df4["rlnCoordinateZ"]

viewer = napari.Viewer(ndisplay=3)

bf = 10 # bin factor

particle1_crds = [[xi / bf, yi / bf, zi / bf] for xi, yi, zi in zip(x1, y1, z1)]
particle4_crds = [[xi / bf, yi / bf, zi / bf] for xi, yi, zi in zip(x4, y4, z4)]

particles11 = viewer.add_points(particle1_crds, face_color='red')
particles44 = viewer.add_points(particle4_crds, face_color='blue')

napari.run()