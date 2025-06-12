import numpy as np
import argparse
import napari
import iotools
from extract_particle_info import (
    get_particle_info,
    extract_centeredcoodinates_by_tomoname,
)


parser = argparse.ArgumentParser(description="")
parser.add_argument(
    "--tomogram", type=str, required=True, help="Specify the tomogram"
)
parser.add_argument(
    "--particlefile1", type=str, required=True, help="Input particle file 1"
)
parser.add_argument(
    "--particlefile2", type=str, required=False, default=None,
    help="Input particle file 2"
)
parser.add_argument(
    "--block",
    type=str,
    required=False,
    default="particles",
    help="Block name in the particlefile",
)
parser.add_argument(
    "--tomogram_name",
    type=str,
    required=True,
    help="Give the tomogram name to select from particle file",
)
args = parser.parse_args()


def get_map(mapname):
    m1 = iotools.Map(mapname)
    m1.read()
    print(m1.workarr.shape)
    return m1.workarr, m1.workcell


def get_particle_coordinates(
    particle_file, block_name, tomogram_name, picle_size, offset
):

    df = get_particle_info(
        particle_file=particle_file, block_name=block_name
    )

    # Centered coordinates
    coords_angstroms = extract_centeredcoodinates_by_tomoname(
        df, tomogram_name
    )

    coords_angstroms = \
        np.array(coords_angstroms, dtype=np.float32) + np.array(offset)
    coords_pxl = coords_angstroms / pixel_size

    return coords_pxl


if __name__ == "__main__":

    arr, cell = get_map(args.tomogram)
    pixel_size = cell[0] / arr.shape[0]
    offset = [x / 2 for x in cell[:3]]

    coords_file1 = get_particle_coordinates(
        particle_file=args.particlefile1,
        block_name=args.block,
        tomogram_name=args.tomogram_name,
        picle_size=pixel_size,
        offset=offset,
    )
    if args.particlefile2:
        coords_file2 = get_particle_coordinates(
            particle_file=args.particlefile2,
            block_name=args.block,
            tomogram_name=args.tomogram_name,
            picle_size=pixel_size,
            offset=offset,
        )

    viewer = napari.Viewer(ndisplay=3)
    new_layer = viewer.add_image(arr)
    particles1 = viewer.add_points(
        coords_file1.tolist(), face_color="green", size=5
    )

    if args.particlefile2:
        particles2 = viewer.add_points(
            coords_file2.tolist(), face_color="red", size=5
        )

    napari.run()
