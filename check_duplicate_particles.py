# Check for duplicate particles in relion particles.star

import starfile
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('--particles', type=str, required=True, help='Input particles.star file')
args = parser.parse_args()



df = starfile.read(args.particles)


x = df["rlnCoordinateX"]
y = df["rlnCoordinateY"]
z = df["rlnCoordinateZ"]

print(len(x), len(list(set(x))))