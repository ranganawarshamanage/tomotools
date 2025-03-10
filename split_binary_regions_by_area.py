from skimage import measure
from more_itertools import sort_together
# from emda2.core import iotools
import argparse
# from scipy.ndimage import binary_erosion
import numpy as np
import mrcfile as mrc


parser = argparse.ArgumentParser(description="")
parser.add_argument("--segmented_tomogram", type=str, required=True, help="Specify segmented_tomogram.mrc")
parser.add_argument("--max_components", 
                    type=int, required=False, default=50,
                    help="maximum number of components (default=10 or less)")
args = parser.parse_args()


class Map:
    def __init__(self, name):
        self.name = name
        self.arr = None
        self.cell = None
        self.origin = None
        self.axorder = None
        self.paddedarr = None
        self.newcell = None

    def read(self):
        try:
            file = mrc.open(self.name)
            axorder = (
                file.header.mapc - 1,
                file.header.mapr - 1,
                file.header.maps - 1,
            )
            self.axorder = axorder
            # print('order: ', order)
            # axes_order = "".join(["XYZ"[i] for i in axorder])
            # print("Axes order: ", axes_order)
            self.arr = np.moveaxis(
                a=np.asarray(file.data, dtype="float"),
                source=(2, 1, 0),
                destination=axorder,
            )
            unit_cell = np.zeros(6, dtype="float")
            cell = file.header.cella[["x", "y", "z"]]
            unit_cell[:3] = cell.astype(
                [("x", "<f4"), ("y", "<f4"), ("z", "<f4")]
            ).view(("<f4", 3))
            unit_cell[3:] = float(90)
            self.cell = unit_cell
            self.origin = [
                1 * file.header.nxstart,
                1 * file.header.nystart,
                1 * file.header.nzstart,
            ]
            file.close()
            print(self.name, self.arr.shape, self.cell[:3])
            self.preprocess_read()
        except FileNotFoundError as e:
            print(e)

    def write(self):
        if self.name == "":
            filename = "new.mrc"
        else:
            filename = self.name
        if self.axorder is None:
            self.axorder = (0, 1, 2)
        if self.origin is None:
            self.origin = [0.0, 0.0, 0.0]
        self.arr = np.moveaxis(
            a=self.arr, source=self.axorder, destination=(2, 1, 0)
        )
        file = mrc.new(
            name=filename,
            data=np.float32(self.arr),
            compression=None,
            overwrite=True,
        )
        file.header.cella.x = self.cell[0]
        file.header.cella.y = self.cell[1]
        file.header.cella.z = self.cell[2]
        file.header.nxstart = self.origin[0]
        file.header.nystart = self.origin[1]
        file.header.nzstart = self.origin[2]
        file.close()

    def preprocess_read(self):
        arr = self.arr
        pixsize = self.cell[0] / self.arr.shape[0]
        tdim = [nd + 1 if nd % 2 != 0 else nd for nd in self.arr.shape]
        self.workcell = np.asarray(
            [pixsize * dim for dim in tdim] + [90.0 for _ in range(3)], "float"
        )
        self.workarr = padimage(arr, tdim)


def padimage(arr, tdim):
    if len(tdim) == 3:
        tnx, tny, tnz = tdim
    elif len(tdim) < 3:
        tnx = tny = tnz = tdim[0]
    else:
        raise SystemExit("More than 3 dimensions given. Cannot handle")
    nx, ny, nz = arr.shape
    assert tnx >= nx
    assert tny >= ny
    assert tnz >= nz
    dz = (tnz - nz) // 2 + (tnz - nz) % 2
    dy = (tny - ny) // 2 + (tny - ny) % 2
    dx = (tnx - nx) // 2 + (tnx - nx) % 2
    image = np.zeros((tnx, tny, tnz), arr.dtype)
    xstart, ystart, zstart = [0 if px == 1 else px for px in [dx, dy, dz]]
    xend, yend, zend = xstart + nx, ystart + ny, zstart + nz
    image[xstart:xend, ystart:yend, zstart:zend] = arr
    return image


def split_connected_regions(binary_array, cell, max_components, fid):
    pixel_size = cell[0] / binary_array.shape[0]
    blobs_labels, nlabels = measure.label(
        binary_array, 
        background=0, 
        connectivity=binary_array.ndim, 
        return_num=True
    )

    if nlabels < max_components:
        max_components = nlabels
        print("maximum number of components = ", max_components)

    regionprops = measure.regionprops(blobs_labels)

    blob_number = []
    blob_area = []
    blob_bbx = []
    blob_inertia = []
    for i in range(nlabels):
        blob_number.append(i+1)
        blob_area.append(regionprops[i].area)
        blob_bbx.append(regionprops[i].bbox)
        blob_inertia.append(regionprops[i].inertia_tensor_eigvals)

    (sblob_area, 
     sblob_number, 
     sblob_bbx, 
     sblob_inertia) = sort_together(
         [blob_area, blob_number, blob_bbx, blob_inertia], 
         reverse=True
         )

    
    for i in range(max_components):
        blob = binary_array * (blobs_labels == sblob_number[i])
        bbx = sblob_bbx[i]
        x =  bbx[3] - bbx[0]
        y =  bbx[4] - bbx[1]
        z =  bbx[5] - bbx[2]
        print(f"xcomponent_{i}.mrc: area: {sblob_area[i]}, bbox: {x, y, z}, inertiaT: {sblob_inertia[i]}")
        fid.write(f"xcomponent_{i}.mrc: area: {sblob_area[i]}, bbox: {x, y, z}, inertiaT: {sblob_inertia[i]}\n")
        fid.write(f"xcomponent_{i}.mrc x0y0z0x1y1z1: {bbx}\n")

        m2 = Map(f"xcomponent_{i}.mrc")
        m2.arr = blob[bbx[0]:bbx[3], bbx[1]:bbx[4], bbx[2]:bbx[5]]
        m2.cell = [pixel_size*sh for sh in [x, y, z]]
        m2.write()  
       

if __name__ == "__main__":
    # Read binary mask
    binary_mask_map = args.segmented_tomogram
    n_components = args.max_components

    fid = open("component_metadata.txt", "w")

    fid.write(f"Input tomogram: {binary_mask_map}\n")
    fid.write(f"Requested Nr. of componenets: {n_components}\n")

    m1 = Map(binary_mask_map)
    m1.read()

    # Make the tomogram binary
    amax = np.amax(m1.workarr)
    amin = np.amin(m1.workarr)
    if abs(amax - float(1)) > 1e-5:
        print(f"Max value of input is {amax}")
        print("Check input")

    # Ines data
    binary_data = (m1.workarr == -127).astype(int)
    m1.workarr = binary_data


    fid.write(f"Cell of the input tomogram: {m1.workcell}\n")
    fid.write(f"Shape of the input tomogram: {m1.workarr.shape}\n")

    split_connected_regions(
        binary_array=m1.workarr,
        cell=m1.workcell,
        max_components=n_components,
        fid=fid,
    )


