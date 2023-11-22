# Create a spherical mask with soft edges for particles

import sys
import numpy as np
from emda2.core import iotools, plotter
import emda2.emda_methods2 as em


def create_soft_edged_mask(boxsize, diameter=None):
    # Create soft-edged-kernel. r1 is the radius of kernel in pixels
    from math import sqrt, cos

    mask = np.zeros(shape=(boxsize, boxsize, boxsize), dtype="float")
    if diameter is None:
        diameter = boxsize

    kernel = np.zeros(shape=(diameter, diameter, diameter), dtype="float")

    kx = ky = kz = diameter
    center = diameter // 2
    # print("center: ", center)
    r1 = diameter // 2
    r0 = r1 - 2
    # print("r1: ", r1, "r0: ", r0)
    for i in range(kx):
        for j in range(ky):
            for k in range(kz):
                dist = sqrt(
                    (i - center) ** 2 + (j - center) ** 2 + (k - center) ** 2
                )
                if dist < r1:
                    if dist < r0:
                        kernel[i, j, k] = 1
                    else:
                        kernel[i, j, k] = (
                            (1 + cos(np.pi * (dist - r0) / (r1 - r0)))
                        ) / 2.0
    # place kernel inside the box
    dx = (boxsize - diameter) // 2
    mask[dx : dx + diameter, dx : dx + diameter, dx : dx + diameter] = kernel
    return mask


if __name__ == "__main__":
    mapname = sys.argv[1]

    """ path = "/Users/Rangana/Tomo/eman2_projects/pipeline_components/extract-average-ribos/coordinate_test/"
    p1 = iotools.Map(
        path + "data/1_data.mrc"
    )  # Relion extracted particle - TS_01
    p1.read()
    p2 = iotools.Map(
        path + "emda_extracted_particles/ribos_0.mrc"
    )  # Manually extracted particle - TS_01.mrc
    p2.read()
    assert p1.workarr.shape[0] == p2.workarr.shape[0] """

    p1 = iotools.Map(mapname)  # Relion extracted particle
    p1.read()

    # create a spherical mask based on the particle's bozxsize
    boxsize = p1.workarr.shape[0]
    mask = create_soft_edged_mask(boxsize, diameter=int(boxsize * 0.5))

    m1 = iotools.Map("sphere_mask.mrc")
    m1.cell = p1.workcell
    m1.arr = mask
    m1.write()

    exit()

    # calculate FSC between two particles
    nbin, res_arr, bin_idx, sgrid = em.get_binidx(
        cell=p1.workcell, arr=p1.workarr
    )
    binfsc = em.fsc(
        f1=np.fft.fftshift(np.fft.fftn(p1.workarr * mask)),
        f2=np.fft.fftshift(np.fft.fftn(p2.workarr * mask)),
        bin_idx=bin_idx,
        nbin=nbin,
    )
    for i in range(len(binfsc)):
        print(i, res_arr[i], binfsc[i])

    plotter.plot_nlines(res_arr, [binfsc])
