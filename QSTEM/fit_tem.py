from PIL import Image
import numpy as np
import struct
import random
import commands

# ----- input parameters ----- start ----- 
mixing_para = 0.05
max_move_atom_init = 1.0 # Angstrom unit
max_move_atom_low_limit = 0.01 # Angstrom unit
max_move_atom_high_limit = 1.5 # Angstrom unit

exchange_ratio = 0.001 # exhange element A for element B

dif_value_limit = 0.7 # exp(0.7) = 2.01
e_atom_limit = 1.1 # evaluate move atom

cycle_limit = 1.0e-8
reset_value = cycle_limit*1000

e_ran_low_limit = 0.995 # => not adopt new position
dif_value_factor = ((256/2)/2)**2

file_name = "tmp.cfg"
reference_image_file = "reference.png"
im = np.array(Image.open(reference_image_file).convert('L'),'f')
# ----- input parameters ----- end ----- 
commands.getoutput("cp "+file_name+" tmp.cfg")

# ----- auto reading ----- start ----- 
La = commands.getoutput("awk 'NR==3{printf \"%12.6f \", $3}' "+file_name)
print "lattice constant a: ", La, " Angstrom"
move_percent = max_move_atom_init / float(La) # [%/100] unit
move_percent_low_limit = max_move_atom_low_limit / float(La) # [%/100] unit
move_percent_high_limit = max_move_atom_high_limit / float(La) # [%/100] unit
# -----
total_char_line = commands.getoutput("grep -c [A-Za-z] "+file_name)
total_num_char_line = int(total_char_line) - 13
print "total number of elements: ", total_num_char_line
# -----
total_cfg_line = commands.getoutput("awk 'END{print NR}' "+file_name)
total_num_cfg_line = int(total_cfg_line) - 13 - total_num_char_line*2
print "total number of atoms: ", total_num_cfg_line
total_number_of_atoms = total_num_cfg_line
# -----
num_char_line = []
number_of_atoms_temp = []
for i in range(total_num_char_line):
  char_line = commands.getoutput("grep -n [A-Za-z] "+file_name+" | sed -e 's/:.*//g' | sed -ne '"+str(14+i)+"p'  | sed -r -e ':loop;N;b loop;s/\\n/ /g' -e 's/ +/ /g'")
  num_char_line.append(int(char_line))
  number_of_atoms_temp.append(int(char_line)-13-(i+1)*2)
number_of_atoms_temp.append(total_number_of_atoms)
print "boundary: ", number_of_atoms_temp
# ----- auto reading ----- end ----- 

def dif2sum(new_data,reference_data):
  return np.sum((new_data - reference_data)**2)

old_value = float((256-1)^2*2 + (256-1)^2*(256*256-2))
new_value = float((256-128)^2*256*256)
dif_value = (new_value-old_value)/old_value * dif_value_factor
n_ran = 16
e_atom = np.ones(total_number_of_atoms)
start_flag = 1
while (abs(dif_value) > cycle_limit):
  commands.getoutput("cp tmp.cfg case.cfg")
  #
  n_ran = n_ran - 15
  j = 0
  for i in range(1,total_num_char_line):
    if (n_ran > (number_of_atoms_temp[i]+2*i)):
      j = i
  n_ran = n_ran - 2*j
  #
  e_atom[n_ran] = e_atom[n_ran]*(1.0-mixing_para) + e_atom[n_ran]*np.exp(-dif_value)*mixing_para
  if (e_atom[n_ran] >= e_atom_limit):
    e_atom[n_ran] = e_atom_limit
    print "----- reset e_atom", e_atom[n_ran]
  print "----- move atom number and e_atom parameter value: ", n_ran, e_atom[n_ran]
  e_ran = random.uniform(e_ran_low_limit,1.0)
  n_ran = random.randrange(1,total_number_of_atoms)
  e_atom_ran = random.uniform(0.0,1.0)
  while (e_atom[n_ran] < e_atom_ran):
    n_ran = random.randrange(1,total_number_of_atoms)
  print " "
  print "----- new move atom number:", n_ran
  move_percent = move_percent*(1.0-mixing_para) + move_percent*e_atom[n_ran]*mixing_para
  if (move_percent <= move_percent_low_limit):
    move_percent == move_percent_low_limit
    print "reset move_percent (low limit)", move_percent_low_limit
  if (move_percent >= move_percent_high_limit):
    move_percent == move_percent_high_limit
    print "reset move_percent (high limit)", move_percent_high_limit
  print "----- moving limit (x,y,z), (1/La) unit: ", move_percent
  x_ran = random.uniform(-move_percent,move_percent)
  y_ran = random.uniform(-move_percent,move_percent)
  z_ran = random.uniform(-move_percent,move_percent)
  #
  n_ran = n_ran + 15
  j = 0
  for i in range(1,total_num_char_line):
    if (n_ran > (15+number_of_atoms_temp[i])):
      j = i
  n_ran = n_ran + 2*j
  #
  ex_ran = random.uniform(0.0,1.0)
  if (exchange_ratio >= ex_ran):
    ex_area_ran1 = 0
    ex_area_ran2 = 0
    while (ex_area_ran1 == ex_area_ran2):
      ex_area_ran1 = random.randrange(0,total_num_char_line)
      ex_area_ran2 = random.randrange(0,total_num_char_line)
    print "exchange area 1: ", ex_area_ran1
    print "exchange area 2: ", ex_area_ran2
    ex_at1_ran = random.randrange(number_of_atoms_temp[ex_area_ran1]+1,number_of_atoms_temp[ex_area_ran1+1]) + 15 + 2*ex_area_ran1
    ex_at2_ran = random.randrange(number_of_atoms_temp[ex_area_ran2]+1,number_of_atoms_temp[ex_area_ran2+1]) + 15 + 2*ex_area_ran2
    print "----- exchange (x y z) in ",ex_at1_ran," line for ",ex_at2_ran
    print "before"
    char_ex_at1 =  commands.getoutput("sed -n "+str(ex_at1_ran)+"p case.cfg")
    print char_ex_at1
    char_ex_at2 =  commands.getoutput("sed -n "+str(ex_at2_ran)+"p case.cfg")
    print char_ex_at2
    print "after"
    new_char_ex_at1 = char_ex_at2[0:38] + " " + char_ex_at1[39:len(char_ex_at1)]
    new_char_ex_at2 = char_ex_at1[0:38] + " " + char_ex_at2[39:]
    print new_char_ex_at1
    print new_char_ex_at2
    replace1 = commands.getoutput("sed -i -e 's/"+str(char_ex_at1)+"/"+str(new_char_ex_at1)+"/g' case.cfg")
    replace2 = commands.getoutput("sed -i -e 's/"+str(char_ex_at2)+"/"+str(new_char_ex_at2)+"/g' case.cfg")
  #
  print "----- replace parameter in cfg file -----"
  print str(n_ran), "line"
  char_old =  commands.getoutput("sed -n "+str(n_ran)+"p case.cfg")
  print char_old
  #
  new_x_char = commands.getoutput("awk 'NR=="+str(n_ran)+"{printf \"%12.6f \", ($1+"+str(x_ran)+"-int($1+"+str(x_ran)+"))}' case.cfg")
  if (float(new_x_char) < 0):
    new_x = float(new_x_char) + 1.0
  else:
    new_x = float(new_x_char)
  new_y_char = commands.getoutput("awk 'NR=="+str(n_ran)+"{printf \"%12.6f \", ($2+"+str(y_ran)+"-int($2+"+str(y_ran)+"))}' case.cfg")
  if (float(new_y_char) < 0):
    new_y = float(new_y_char) + 1.0
  else:
    new_y = float(new_y_char)
  new_z_char = commands.getoutput("awk 'NR=="+str(n_ran)+"{printf \"%12.6f \", ($3+"+str(z_ran)+"-int($3+"+str(z_ran)+"))}' case.cfg")
  if (float(new_z_char) < 0):
    new_z = float(new_z_char) + 1.0
  else:
    new_z = float(new_z_char)
  #
  char_new = commands.getoutput("awk 'NR=="+str(n_ran)+"{printf \"%12.6f %12.6f %12.6f %12.6f %12.6f \\n\", "+str(new_x)+", "+str(new_y)+", "+str(new_z)+", $4, $5}' case.cfg")
  print char_new
  #
  replace = commands.getoutput("sed -i -e 's/"+str(char_old)+"/"+str(char_new)+"/g' case.cfg")
  #print replace
  print "----- qstem calculation -----"
  run = commands.getoutput("$HOME/QSTEM/bin/stem3 tem.qsc")
  #print run
  copy = commands.getoutput("cp -f ./case/case_Proj.img case.img")
  #print copy
  # ----------------------------------------------------
  print "----- comparing calc image and exp -----"
  fid = open("case.img","rb")
  # int data
  header_size = struct.unpack("i",fid.read(4))
  #print "header size: ", header_size[0]
  #
  paramSize = struct.unpack("i",fid.read(4))
  #print "parameter size: ", paramSize[0]
  #
  commentSize = struct.unpack("i",fid.read(4))
  #print "comment size: ", commentSize[0]
  #
  Nx = struct.unpack("i",fid.read(4))
  #print "Nx: number of pixels in x-direction, ", Nx[0]
  #
  Ny = struct.unpack("i",fid.read(4))
  #print "Ny: number of pixels in y-direction, ", Ny[0]
  #
  complexFlag = struct.unpack("i",fid.read(4))
  #print "Flag: complex(1) or real(0): ", complexFlag[0]
  #
  dataSize = struct.unpack("i",fid.read(4))
  #print "Data size (byte unit): ", dataSize[0]
  #
  version = struct.unpack("i",fid.read(4))
  #print "version: ",version[0]
  #
  # double data
  t = struct.unpack("d",fid.read(8))
  #print "sample thickness or defocus: ", t[0]
  #
  dx = struct.unpack("d",fid.read(8))
  #print "dx: pixel size along x-direction, ", dx[0], " [Angstrom]"
  #
  dy = struct.unpack("d",fid.read(8))
  #print "dy: pixel size along y-direction, ", dy[0], " [Angstrom]"
  #
  # other data
  for i in range(0,paramSize[0]):
    auxilliary_data = struct.unpack("d",fid.read(8))
    print auxilliary_data[0]
  #
  comment = []
  for i in xrange(0,commentSize[0]):
    comment.append(struct.unpack("c",fid.read(1))[0])
  #print "comment"
  #print comment
  #
  # image data
  #fod = open("case.raw","wb")
  actual_data = [[0 for i in range(Nx[0])] for j in range(Ny[0])]
  for i in xrange(0,Nx[0]):
    for j in xrange(0,Ny[0]):
      actual_data[i][j] = struct.unpack("f",fid.read(dataSize[0]))[0]
      #
      #data = fid.read(dataSize[0])
      #fod.write(data)
  #
  data = np.array(actual_data,dtype="float32")
  #print data.dtype
  #print data.ndim
  #print data.shape
  #
  #-----matplotlib-----
  #plt.imshow(data)
  #plt.show()
  #
  #-----PIL-----
  new_data = (255 - ((data - np.min(data)) / (np.max(data) - np.min(data))) * 255).astype(np.uint8)
  pil_img = Image.fromarray(new_data)
  #pil_img.show()
  #new_pil_img = pil_img.convert("L")
  #new_pil_img.show()
  #new_pil_img.save("case.png")
  #
  fid.close()
  #fod.close()
  # ----------------------------------------------------
  #
  new_value = dif2sum(new_data,im)
  print "----- new value: ", new_value
  print "----- old value: ", old_value
  dif_value = (new_value-old_value)/old_value * dif_value_factor
  print "----- dif value: ", dif_value
  if (dif_value == 0.0 or dif_value >= dif_value_limit):
    dif_value = reset_value 
    print "----- if value = 0, reset value and exp(-reset): ", dif_value, np.exp(-dif_value)
    reset_flag = 1
  else:
    print "----- value and exp(-value): ", dif_value, np.exp(-dif_value)
    reset_flag = 0
  print "----- e_ran value:", e_ran
  if (np.exp(-dif_value) >= e_ran and reset_flag == 0):
    print "----- adopt new position -----"
    commands.getoutput("cp case.cfg tmp.cfg")
    old_value = new_value
  else: 
    print "----- not adopt new position -----"
  if (start_flag == 1):
    old_value = new_value
    start_flag = 0 
