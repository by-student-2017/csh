from scipy import optimize
import numpy as np
import commands
import sys

file_tmp = 'mg2nih4_tmp.gin'
file_inp = 'mg2nih4_nm.gin'
file_out = 'mg2nih4_nm.got'

gulp_adress = "$HOME/gulp-5.1/Src/gulp"
timeout_s   = "timeout -s 9 5s" # kill command after 5 s

output_file_name = "nm_files"
commands.getoutput("rm -f -r "+output_file_name)
commands.getoutput("mkdir "+output_file_name)
count = 0

select_limit_y = 1000
#limit = 0.001

keyword1 = "Final"
keyword2 = "cell"
target = [7.8827, 7.8827, 6.5006, 68.6808, 111.3192, 131.9938, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001]
y_str = [0,0,0,0,0,0,0,0,0,0,0,0] # initial data

# fitting parameters
x0 = [1.7583, 2.963, -1.044, 1079.8, 918.2, 5420.3, 0.25037, 0.24428, 0.24337, 5.04, 5.36, 38.78, 22.9] # initial data
b1 = np.array([[1.5,2.2],[2.5,3.2],[-0.85,-1.43],[1000,1200],[800,1400],[5400,5500],[0.2,0.3],[0.2,0.3],[0.2,0.3],[0.2,0.3],[0,50],[0,50],[0,50],[0,50]]) # boundary

def f(x):

  print "------------------------"
  global count
  count += 1
  print count

  global no_flag
  no_flag = 0

  fi = open(file_tmp,'r')
  text = fi.read().replace('XXX0',str(x[0]).replace("[","").replace("]",""))
  text = text.replace('XXX1',str(x[1]).replace("[","").replace("]",""))
  text = text.replace('XXX2',str(x[2]).replace("[","").replace("]",""))
  xc3 = (-(x[0]*3*2/3+x[1]) - x[2]*4) / 4
  text = text.replace('XXX3',str(xc3).replace("[","").replace("]",""))
  text = text.replace('XXX4',str(x[3]).replace("[","").replace("]",""))
  text = text.replace('XXX5',str(x[4]).replace("[","").replace("]",""))
  text = text.replace('XXX6',str(x[5]).replace("[","").replace("]",""))
  text = text.replace('XXX7',str(x[6]).replace("[","").replace("]",""))
  text = text.replace('XXX8',str(x[7]).replace("[","").replace("]",""))
  text = text.replace('XXX9',str(x[8]).replace("[","").replace("]",""))
  text = text.replace('XX10',str(x[9]).replace("[","").replace("]",""))
  text = text.replace('XX11',str(x[10]).replace("[","").replace("]",""))
  text = text.replace('XX12',str(x[11]).replace("[","").replace("]",""))
  #text = text.replace('XX12','0.0')
  text = text.replace('XX13',str(x[12]).replace("[","").replace("]",""))
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
    y = 99999.9
    print "time out or bad condition range"
  else:
    commands.getoutput("cp -p "+file_inp+" ./"+output_file_name+"/"+file_inp+"."+str(count))
    commands.getoutput("cp -p "+file_out+" ./"+output_file_name+"/"+file_out+"."+str(count))
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
    print y_str
    if (y_str[11] == "************"):
      y_str[11] = "99999.999999"
      no_flag = 1

    y = 500*((float(y_str[0])-target[0])/target[0])**2 + 500*((float(y_str[1])-target[1])/target[1])**2 \
    + 500*((float(y_str[2])-target[2])/target[2])**2 + 50*((float(y_str[3])-target[3])/target[3])**2 \
    + 50*((float(y_str[4])-target[4])/target[4])**2 + 50*((float(y_str[5])-target[5])/target[5])**2 \
    + ((float(y_str[6])-target[6])/target[6])**2 + ((float(y_str[7])-target[7])/target[7])**2 \
    + ((float(y_str[8])-target[8])/target[7])**2 + ((float(y_str[9])-target[9])/target[9])**2 \
    + ((float(y_str[10])-target[10])/target[8])**2 + ((float(y_str[11])-target[11])/target[11])**2

    if (y < select_limit_y and no_flag == 0):
      commands.getoutput("cp - p "+file_inp+" ./"+output_file_name+"/"+file_inp+"."+str(count))
      commands.getoutput("cp - p "+file_out+" ./"+output_file_name+"/"+file_out+"."+str(count))

  print y, x
  print "------------------------"

  #if (y < limit):
  #  sys.exit()

  return y

res = optimize.minimize(f,x0,method='Nelder-Mead',bounds=b1)
#res = optimize.minimize(f,x0,method='TNC',bounds=b1)
#res = optimize.minimize(f,x0,method='Powell',bounds=b1)
#res = optimize.minimize(f,x0,method='BFGS',bounds=b1)
