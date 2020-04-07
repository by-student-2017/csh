from scipy import optimize
import numpy as np
import commands

file_tmp = 'mgh2_tmp.gin'
file_inp = 'mgh2.gin'
file_out = 'mgh2.got'

gulp_adress = "$HOME/gulp-5.1/Src/gulp"

keyword = "Frequency"
target1 = 180.65
target2 = 238.18
x0 = [1.5,-0.5]
b1 = np.array([[0,2.2],[-1,1]])

def f(x):
  
  fi = open(file_tmp,'r')
  text = fi.read().replace('XXX0',str(x[0]).replace("[","").replace("]",""))
  text = text.replace('XXX1',str(x[1]).replace("[","").replace("]",""))
  x2 = -(x[0]/2 + x[1])
  text = text.replace('XXX2',str(x2).replace("[","").replace("]",""))
  fi.close
  
  with open(file_inp,'w') as f:
    print >> f, text
  
  commands.getoutput(gulp_adress+" < "+file_inp+" > "+file_out)
  
  get_data = "awk '{if($1==\""+keyword+"\" && NR<=420){printf($5)}else{}}' mgh2.got"
  y1_str = commands.getoutput(get_data) 
  get_data = "awk '{if($1==\""+keyword+"\" && NR<=420){printf($6)}else{}}' mgh2.got"
  y2_str = commands.getoutput(get_data) 
  
  y = ((float(y1_str)-target1)/target1)**2 + ((float(y2_str)-target2)/target2)**2
  print y, x
  
  return y

res = optimize.minimize(f,x0,method='Nelder-Mead',bounds=b1)
#res = optimize.fmin_bfgs(f,x0)
#res = optimize.fmin_cg(f,x0)

