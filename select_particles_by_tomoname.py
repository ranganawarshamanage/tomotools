import starfile
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('--particles', type=str, required=True, help='Input particles.star file')
parser.add_argument('--tomoname', type=str, required=True, help='Specify tomogram name')
parser.add_argument('--binfactor', type=float, required=True, help='Specify the bin factor')
args = parser.parse_args()

path = './'
#particles = "/ceph/scheres_grp/ranganaw/TomoTutorial/PickTomo/job014/particles.star"
#outfname = 'TS_54_particles.star'
#binfactor = 7.4

particles = args.particles
binfactor = args.binfactor
tomoname = args.tomoname
outfname = tomoname + "_particles.star"

df = starfile.read(particles)

x = []
y = []
z = []
classname = []
tnamelist = df["rlnTomoName"]
for i, tname in enumerate(tnamelist):
    if df["rlnTomoName"][i] == tomoname:
        x.append(df["rlnCoordinateX"][i] // binfactor)
        y.append(df["rlnCoordinateY"][i] // binfactor)
        z.append(df["rlnCoordinateZ"][i] // binfactor)
        classname.append(df["rlnClassNumber"][i])

df2 = pd.DataFrame()
df2["rlnTomoName"] = [tomoname for _ in range(len(x))]
df2["rlnCoordinateX"] = [ix for ix in x]
df2["rlnCoordinateY"] = [iy for iy in y]
df2["rlnCoordinateZ"] = [iz for iz in z]
df2["rlnClassNumber"] = classname

starfile.write(df2, path + outfname, overwrite=True)
