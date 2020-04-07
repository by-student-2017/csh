from scipy import optimize
import numpy as np
import commands
import sys

file_tmp = 'mgh2_tmp.gin'
file_inp = 'mgh2.gin'
file_out = 'mgh2.got'

gulp_adress = "$HOME/gulp-5.1/Src/gulp"

limit = 0.2143

keyword = "Frequency"
target = [180.65, 238.18, 299.64, 299.64, 488.61, 709.34, 738.37, 878.45, 988.57, 988.57, 1121.2, 1206.34, 1309.89, 1396.78, 1472.68]
y_str = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] # initial data

# fitting parameters
x0 = [1.4, -1.431, 1024.0, 5585.0, 0.2442, 0.2353, 10.3, 0.0, 600.0, 999999, 1.4] # initial data
b1 = np.array([[0,2.2],[-1,1],[950,1800],[5500,6000],[0.2,0.3],[0.2,0.3],[0,50],[0,50],[0,900],[9999,9999999999],[0,2.2]]) # boundary

def f(x):

  fi = open(file_tmp,'r')
  text = fi.read().replace('XXX0',str(x[0]).replace("[","").replace("]",""))
  xc1 = x[10] - x[0]
  text = text.replace('XXX1',str(xc1).replace("[","").replace("]",""))
  text = text.replace('XXX2',str(x[1]).replace("[","").replace("]",""))
  xc2 = -x[10]/2 - x[1]
  text = text.replace('XXX3',str(xc2).replace("[","").replace("]",""))
  text = text.replace('XXX4',str(x[2]).replace("[","").replace("]",""))
  text = text.replace('XXX5',str(x[3]).replace("[","").replace("]",""))
  text = text.replace('XXX6',str(x[4]).replace("[","").replace("]",""))
  text = text.replace('XXX7',str(x[5]).replace("[","").replace("]",""))
  text = text.replace('XXX8',str(x[6]).replace("[","").replace("]",""))
  text = text.replace('XXX9',str(x[7]).replace("[","").replace("]",""))
  text = text.replace('XX10',str(x[8]).replace("[","").replace("]",""))
  text = text.replace('XX11',str(x[9]).replace("[","").replace("]",""))
  fi.close

  with open(file_inp,'w') as f:
    print >> f, text

  commands.getoutput(gulp_adress+" < "+file_inp+" > "+file_out)

  get_line = "awk '{if($1==\""+keyword+"\"){printf \"%i \",NR }else{}}' "+file_out
  line = commands.getoutput(get_line)
  line_list = line.rstrip().split(" ")
  #print line_list, len(line_list)

  get_data = "awk '{if($1==\""+keyword+"\" && NR=="+line_list[0]+"){printf($5)}else{}}' "+file_out
  y_str[0] = commands.getoutput(get_data)
  get_data = "awk '{if($1==\""+keyword+"\" && NR=="+line_list[0]+"){printf($6)}else{}}' "+file_out
  y_str[1] = commands.getoutput(get_data)
  get_data = "awk '{if($1==\""+keyword+"\" && NR=="+line_list[0]+"){printf($7)}else{}}' "+file_out
  y_str[2] = commands.getoutput(get_data)
  #
  get_data = "awk '{if($1==\""+keyword+"\" && NR=="+line_list[1]+"){printf($2)}else{}}' "+file_out
  y_str[3] = commands.getoutput(get_data)
  get_data = "awk '{if($1==\""+keyword+"\" && NR=="+line_list[1]+"){printf($3)}else{}}' "+file_out
  y_str[4] = commands.getoutput(get_data)
  get_data = "awk '{if($1==\""+keyword+"\" && NR=="+line_list[1]+"){printf($4)}else{}}' "+file_out
  y_str[5] = commands.getoutput(get_data)
  get_data = "awk '{if($1==\""+keyword+"\" && NR=="+line_list[1]+"){printf($5)}else{}}' "+file_out
  y_str[6] = commands.getoutput(get_data)
  get_data = "awk '{if($1==\""+keyword+"\" && NR=="+line_list[1]+"){printf($6)}else{}}' "+file_out
  y_str[7] = commands.getoutput(get_data)
  get_data = "awk '{if($1==\""+keyword+"\" && NR=="+line_list[1]+"){printf($7)}else{}}' "+file_out
  y_str[8] = commands.getoutput(get_data)
  #
  get_data = "awk '{if($1==\""+keyword+"\" && NR=="+line_list[2]+"){printf($2)}else{}}' "+file_out
  y_str[9] = commands.getoutput(get_data)
  get_data = "awk '{if($1==\""+keyword+"\" && NR=="+line_list[2]+"){printf($3)}else{}}' "+file_out
  y_str[10] = commands.getoutput(get_data)
  get_data = "awk '{if($1==\""+keyword+"\" && NR=="+line_list[2]+"){printf($4)}else{}}' "+file_out
  y_str[11] = commands.getoutput(get_data)
  get_data = "awk '{if($1==\""+keyword+"\" && NR=="+line_list[2]+"){printf($5)}else{}}' "+file_out
  y_str[12] = commands.getoutput(get_data)
  get_data = "awk '{if($1==\""+keyword+"\" && NR=="+line_list[2]+"){printf($6)}else{}}' "+file_out
  y_str[13] = commands.getoutput(get_data)
  get_data = "awk '{if($1==\""+keyword+"\" && NR=="+line_list[2]+"){printf($7)}else{}}' "+file_out
  y_str[14] = commands.getoutput(get_data)
  print y_str

  y = ((float(y_str[0])-target[0])/target[0])**2 + ((float(y_str[1])-target[1])/target[1])**2 \
    + ((float(y_str[2])-target[2])/target[2])**2 + ((float(y_str[3])-target[3])/target[3])**2 \
    + ((float(y_str[4])-target[4])/target[4])**2 + ((float(y_str[5])-target[5])/target[5])**2 \
    + ((float(y_str[6])-target[6])/target[6])**2 + ((float(y_str[7])-target[7])/target[7])**2 \
    + ((float(y_str[8])-target[8])/target[7])**2 + ((float(y_str[9])-target[9])/target[9])**2 \
    + ((float(y_str[10])-target[10])/target[8])**2 + ((float(y_str[11])-target[11])/target[11])**2 \
    + ((float(y_str[12])-target[12])/target[12])**2 + ((float(y_str[13])-target[13])/target[13])**2 \
    + ((float(y_str[14])-target[14])/target[14])**2
  print y, x
  print "------------------------"

  if (y < limit):
    sys.exit()

  return y

#res = optimize.minimize(f,x0,method='Nelder-Mead',bounds=b1)
#res = optimize.minimize(f,x0,method='TNC',bounds=b1)
#res = optimize.minimize(f,x0,method='Powell',bounds=b1)
res = optimize.minimize(f,x0,method='BFGS',bounds=b1)
