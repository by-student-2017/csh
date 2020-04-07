import h5py
import matplotlib.pyplot as plt

f = h5py.File("output.h5","r")

data = f["4DSTEM_simulation/data/realslices/annular_detector_depth0000/realslice"].value
plt.imshow(data)
plt.show()

