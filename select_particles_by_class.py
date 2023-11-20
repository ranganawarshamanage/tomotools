import starfile
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('--particles', type=str, required=True, help='Input particles.star file')
parser.add_argument('--class_number', type=int, required=True, help='Specify class number')
parser.add_argument('--binfactor', type=float, required=True, help='Specify the bin factor')
args = parser.parse_args()

path = './'

particles = args.particles
class_number = args.class_number
binfactor = args.binfactor
outfname = "_class%d_particles.star" % class_number

df = starfile.read(particles)

x = []
y = []
z = []
tname = []
classlist = df["rlnClassNumber"]
for i, clas in enumerate(classlist):
    if df["rlnClassNumber"][i] == class_number:
        tname.append(df["rlnTomoName"][i])
        x.append(df["rlnCoordinateX"][i] // binfactor)
        y.append(df["rlnCoordinateY"][i] // binfactor)
        z.append(df["rlnCoordinateZ"][i] // binfactor)

df2 = pd.DataFrame()
df2["rlnTomoName"] = tname #[tomoname for _ in range(len(x))]
df2["rlnCoordinateX"] = [ix for ix in x]
df2["rlnCoordinateY"] = [iy for iy in y]
df2["rlnCoordinateZ"] = [iz for iz in z]
df2["rlnClassNumber"] = [class_number for _ in range(len(x))]

starfile.write(df2, path + tname[0]+outfname, overwrite=True)
