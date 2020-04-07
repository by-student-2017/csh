from PIL import Image
import numpy as np
import struct
#import matplotlib.pyplot as plt

fid = open("case.img","rb")

# int data
header_size = struct.unpack("i",fid.read(4))
print "header size: ", header_size[0]

paramSize = struct.unpack("i",fid.read(4))
print "parameter size: ", paramSize[0]

commentSize = struct.unpack("i",fid.read(4))
print "comment size: ", commentSize[0]

Nx = struct.unpack("i",fid.read(4))
print "Nx: number of pixels in x-direction, ", Nx[0]

Ny = struct.unpack("i",fid.read(4))
print "Ny: number of pixels in y-direction, ", Ny[0]

complexFlag = struct.unpack("i",fid.read(4))
print "Flag: complex(1) or real(0): ", complexFlag[0]

dataSize = struct.unpack("i",fid.read(4))
print "Data size (byte unit): ", dataSize[0]

version = struct.unpack("i",fid.read(4))
print "version: ",version[0]

# double data

t = struct.unpack("d",fid.read(8))
print "sample thickness or defocus: ", t[0]

dx = struct.unpack("d",fid.read(8))
print "dx: pixel size along x-direction, ", dx[0], " [Angstrom]"

dy = struct.unpack("d",fid.read(8))
print "dy: pixel size along y-direction, ", dy[0], " [Angstrom]"

# other data

for i in range(0,paramSize[0]):
  auxilliary_data = struct.unpack("d",fid.read(8))
  print auxilliary_data[0]

comment = []
for i in xrange(0,commentSize[0]):
  comment.append(struct.unpack("c",fid.read(1))[0])
print "comment"
print comment

# image data

fod = open("case.raw","wb")

actual_data = [[0 for i in range(Nx[0])] for j in range(Ny[0])]
for i in xrange(0,Nx[0]):
  for j in xrange(0,Ny[0]):
    actual_data[i][j] = struct.unpack("f",fid.read(dataSize[0]))[0]
    #
    #data = fid.read(dataSize[0])
    #fod.write(data)

data = np.array(actual_data,dtype="float32")
#print data.dtype
#print data.ndim
#print data.shape

#-----matplotlib-----
#plt.imshow(data)
#plt.show()

#-----PIL-----
new_data = (255 - ((data - np.min(data)) / (np.max(data) - np.min(data))) * 255).astype(np.uint8)
pil_img = Image.fromarray(new_data)
pil_img.show()
#new_pil_img = pil_img.convert("L")
#new_pil_img.show()
#new_pil_img.save("case.png")

fid.close()
fod.close()
