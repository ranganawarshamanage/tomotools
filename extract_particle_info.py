import gemmi
import pandas as pd


def get_particle_info(particle_file, block_name, columns=None):
    # Load the STAR file using gemmi
    star = gemmi.cif.read_file(particle_file)

    # Find the block by name (e.g., 'data_particles')
    block = star.find_block(block_name)

    # Check if the block exists
    if block is None:
        print("Block 'data_particles' not found.")
    else:
        loop = block[0].loop
        cols = loop.tags[0:] if columns is None else columns
        data = []
        for row in block.find(cols):
            data.append([gemmi.cif.as_string(x) for x in row])
        df = pd.DataFrame(data, columns=cols)
        print(df)
        return df


def extract_centeredcoodinates_by_tomoname(df, tomoname):
    x1, y1, z1 = [], [], []
    for i, name in enumerate(df["_rlnTomoName"]):
        if name == tomoname:
            x1.append(float(df["_rlnCenteredCoordinateXAngst"][i]))
            y1.append(float(df["_rlnCenteredCoordinateYAngst"][i]))
            z1.append(float(df["_rlnCenteredCoordinateZAngst"][i]))

    particle_crd_list = [[xi, yi, zi] for xi, yi, zi in zip(x1, y1, z1)]

    return particle_crd_list


def extract_origincoodinates_by_tomoname(df, tomoname):
    x1, y1, z1 = [], [], []
    for i, name in enumerate(df["_rlnTomoName"]):
        if name == tomoname:
            x1.append(float(df["_rlnOriginXAngst"][i]))
            y1.append(float(df["_rlnOriginYAngst"][i]))
            z1.append(float(df["_rlnOriginZAngst"][i]))

    particle_crd_list = [[xi, yi, zi] for xi, yi, zi in zip(x1, y1, z1)]

    return particle_crd_list


def extract_subtomoangles_by_tomoname(df, tomoname):
    subRot, subTilt, subPsi = [], [], []
    for i, name in enumerate(df["_rlnTomoName"]):
        if name == tomoname:
            subRot.append(float(df["_rlnTomoSubtomogramRot"][i]))
            subTilt.append(float(df["_rlnTomoSubtomogramTilt"][i]))
            subPsi.append(float(df["_rlnTomoSubtomogramPsi"][i]))

    angle_list = [[ri, ti, pi] for ri, ti, pi in zip(subRot, subTilt, subPsi)]

    return angle_list


def extract_angles_by_tomoname(df, tomoname):
    subRot, subTilt, subPsi = [], [], []
    for i, name in enumerate(df["_rlnTomoName"]):
        if name == tomoname:
            subRot.append(float(df["_rlnAngleRot"][i]))
            subTilt.append(float(df["_rlnAngleTilt"][i]))
            subPsi.append(float(df["_rlnAnglePsi"][i]))

    angle_list = [[ri, ti, pi] for ri, ti, pi in zip(subRot, subTilt, subPsi)]

    return angle_list
