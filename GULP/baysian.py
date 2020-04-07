from bayes_opt import BayesianOptimization
import numpy as np
import commands
import sys

file_tmp = 'mg2nih4_tmp.gin'
file_inp = 'mg2nih4_b.gin'
file_out = 'mg2nih4_b.got'

gulp_adress = "$HOME/gulp-5.1/Src/gulp"
timeout_s   = "timeout -s 9 5s" # kill command after 5 s

output_file_name = "baysian_files"
commands.getoutput("rm -f -r "+output_file_name)
commands.getoutput("mkdir "+output_file_name)

select_limit_y = 1.0e-12

keyword1 = "Final"
keyword2 = "cell"
target = [7.8827, 7.8827, 6.5006, 68.6808, 111.3192, 131.9938, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001]
y_str = [0,0,0,0,0,0,0,0,0,0,0,0] # initial data

# fitting parameters
x0  =  1.75896
x1  =  2.96068
X2  = -1.09386
x3  = 1080.0
x4  =  918.0
x5  = 5421.0
x6  =  0.25064
x7  =  0.24428
x8  =  0.24345
x9  =  5.04422
x10 =  5.36111
x11 = 38.7798
x12 = 22.2918
pbounds = {
   'x0': (1.6,1.9),
   'x1': (2.5,3.0),
   'x2': (-0.9,-1.2),
   'x3': (1000,1100),
   'x4': (800,1100),
   'x5': (5400,5500),
   'x6': (0.24,0.255),
   'x7': (0.24,0.25),
   'x8': (0.24,0.25),
   'x9': (3,10),
  'x10': (3,10),
  'x11': (30,45),
  'x12': (20,45),
} # boundary

count = 0

def descripter(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12):

  print "------------------------"
  global count
  count += 1
  #print count

  global no_flag
  no_flag = 0

  fi = open(file_tmp,'r')
  text = fi.read().replace('XXX0',str(x0))
  text = text.replace('XXX1',str(x1))
  text = text.replace('XXX2',str(x2))
  xc3 = (-(x0*3*2/3+x1) - x2*4) / 4
  text = text.replace('XXX3',str(xc3))
  text = text.replace('XXX4',str(x3))
  text = text.replace('XXX5',str(x4))
  text = text.replace('XXX6',str(x5))
  text = text.replace('XXX7',str(x6))
  text = text.replace('XXX8',str(x7))
  text = text.replace('XXX9',str(x8))
  text = text.replace('XX10',str(x9))
  text = text.replace('XX11',str(x10))
  text = text.replace('XX12',str(x11))
  #text = text.replace('XX12','0.0')
  text = text.replace('XX13',str(x12))
  fi.close

  with open(file_inp,'w') as f:
    print >> f, text

  commands.getoutput(timeout_s+" "+gulp_adress+" < "+file_inp+" > "+file_out)
  #commands.getoutput("cp -p "+file_inp+" ./"+output_file_name+"/"+file_inp+"."+str(count))
  #commands.getoutput("cp -p "+file_out+" ./"+output_file_name+"/"+file_out+"."+str(count))
  #commands.getoutput("cp -p "+file_inp+" ./"+output_file_name+file_inp+".`date \"+%Y%m%d_%H%M%S\"`")
  #commands.getoutput("cp -p "+file_out+" ./"+output_file_name+file_out+".`date \"+%Y%m%d_%H%M%S\"`")

  get_line = "awk '{if($1==\""+keyword1+"\" && $2==\""+keyword2+"\"){printf \"%i \",NR }else{}}' "+file_out
  line = commands.getoutput(get_line)
  #print line
  if (line == ""):
    #y_str[0] = 99.9
    #y_str[1] = 99.9
    #y_str[2] = 99.9
    #y_str[3] = 99.9
    #y_str[4] = 99.9
    #y_str[5] = 99.9
    #y_str[6] = 99.9
    #y_str[7] = 99.9
    #y_str[8] = 99.9
    #y_str[9] = 99.9
    #y_str[10] = 99.9
    #y_str[11] = 99.9
    y = 0.0
    print "time out or bad condition range"
  else:
    line_list = line.rstrip().split(" ")
    #print line_list, len(line_list)
    a_line = int(line_list[0]) + 3
    b_line = int(line_list[0]) + 4
    c_line = int(line_list[0]) + 5
    alpha_line = int(line_list[0]) + 6
    beta_line  = int(line_list[0]) + 7
    gamma_line = int(line_list[0]) + 8

    get_data = "awk '{if($1==\"a\" && NR=="+str(a_line)+"){printf($2)}else{}}' "+file_out
    y_str[0] = commands.getoutput(get_data)
    if (y_str[0] == "************"):
      y_str[0] = "99999.999999"
      no_flag = 1
    get_data = "awk '{if($1==\"b\" && NR=="+str(b_line)+"){printf($2)}else{}}' "+file_out
    y_str[1] = commands.getoutput(get_data)
    if (y_str[1] == "************"):
      y_str[1] = "99999.999999"
      no_flag = 1
    get_data = "awk '{if($1==\"c\" && NR=="+str(c_line)+"){printf($2)}else{}}' "+file_out
    y_str[2] = commands.getoutput(get_data)
    if (y_str[2] == "************"):
      y_str[2] = "99999.999999"
      no_flag = 1
    #
    get_data = "awk '{if($1==\"alpha\" && NR=="+str(alpha_line)+"){printf($2)}else{}}' "+file_out
    y_str[3] = commands.getoutput(get_data)
    if (y_str[3] == "************"):
      y_str[3] = "99999.999999"
      no_flag = 1
    get_data = "awk '{if($1==\"beta\" && NR=="+str(beta_line)+"){printf($2)}else{}}' "+file_out
    y_str[4] = commands.getoutput(get_data)
    if (y_str[4] == "************"):
      y_str[4] = "99999.999999"
      no_flag = 1
    get_data = "awk '{if($1==\"gamma\" && NR=="+str(gamma_line)+"){printf($2)}else{}}' "+file_out
    y_str[5] = commands.getoutput(get_data)
    if (y_str[5] == "************"):
      y_str[5] = "99999.999999"
      no_flag = 1
    #
    get_data = "awk '{if($1==\"a\" && NR=="+str(a_line)+"){printf($5)}else{}}' "+file_out
    y_str[6] = commands.getoutput(get_data)
    if (y_str[6] == "************"):
      y_str[6] = "99999.999999"
      no_flag = 1
    get_data = "awk '{if($1==\"b\" && NR=="+str(b_line)+"){printf($5)}else{}}' "+file_out
    y_str[7] = commands.getoutput(get_data)
    if (y_str[7] == "************"):
      y_str[7] = "99999.999999"
      no_flag = 1
    get_data = "awk '{if($1==\"c\" && NR=="+str(c_line)+"){printf($5)}else{}}' "+file_out
    y_str[8] = commands.getoutput(get_data)
    if (y_str[8] == "************"):
      y_str[8] = "99999.999999"
      no_flag = 1
    #
    get_data = "awk '{if($1==\"alpha\" && NR=="+str(alpha_line)+"){printf($5)}else{}}' "+file_out
    y_str[9] = commands.getoutput(get_data)
    if (y_str[9] == "************"):
      y_str[9] = "99999.999999"
      no_flag = 1
    get_data = "awk '{if($1==\"beta\" && NR=="+str(beta_line)+"){printf($5)}else{}}' "+file_out
    y_str[10] = commands.getoutput(get_data)
    if (y_str[10] == "************"):
      y_str[10] = "99999.999999"
      no_flag = 1
    get_data = "awk '{if($1==\"gamma\" && NR=="+str(gamma_line)+"){printf($5)}else{}}' "+file_out
    y_str[11] = commands.getoutput(get_data)
    if (y_str[11] == "************"):
      y_str[11] = "99999.999999"
      no_flag = 1

    y = 1 / ( 500*((float(y_str[0])-target[0])/target[0])**2 + 500*((float(y_str[1])-target[1])/target[1])**2 \
    + 500*((float(y_str[2])-target[2])/target[2])**2 + 50*((float(y_str[3])-target[3])/target[3])**2 \
    + 50*((float(y_str[4])-target[4])/target[4])**2 + 50*((float(y_str[5])-target[5])/target[5])**2 \
    + ((float(y_str[6])-target[6])/target[6])**2 + ((float(y_str[7])-target[7])/target[7])**2 \
    + ((float(y_str[8])-target[8])/target[7])**2 + ((float(y_str[9])-target[9])/target[9])**2 \
    + ((float(y_str[10])-target[10])/target[8])**2 + ((float(y_str[11])-target[11])/target[11])**2 )

    if (y > select_limit_y and no_flag == 0):
      commands.getoutput("cp -p "+file_inp+" ./"+output_file_name+"/"+file_inp+"."+str(count))
      commands.getoutput("cp -p "+file_out+" ./"+output_file_name+"/"+file_out+"."+str(count))

  print y, y_str

  return y

optimizer = BayesianOptimization(f=descripter, pbounds=pbounds)
optimizer.maximize(init_points=1, n_iter=2000, acq="ucb")
#acq = ucb, ei, poi, (default: ubc)
