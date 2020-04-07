import random
from deap import creator, base, tools, algorithms
import numpy
import numpy as np
import commands
import sys
#----------------------------------------------------------------------
file_tmp = 'mg2nih4_tmp.gin'
file_inp = 'mg2nih4_ga.gin'
file_out = 'mg2nih4_ga.got'

gulp_adress = "$HOME/gulp-5.1/Src/gulp"
timeout_s   = "timeout -s 9 5s"

output_file_name = "ga_files"
commands.getoutput("rm -f -r "+output_file_name)
commands.getoutput("mkdir "+output_file_name)

select_limit_y = 100000.0

keyword1 = "Final"
keyword2 = "cell"
target = [7.8827, 7.8827, 6.5006, 68.6808, 111.3192, 131.9938, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001, 0.000001]
y_str = [0,0,0,0,0,0,0,0,0,0,0,0] # initial data

# fitting parameters
x0 = [1.73, 2.97, -1.0, 1155.4, 1015.4, 5420.3, 0.2425, 0.2473, 0.2371, 26.94, 4.489, 12.11, 37.3] # initial data
b1 = np.array([[1.6,1.8],[2.9,3],[-0.9,-1.1],[1000,1200],[900,1100],[5400,5500],[0.24,0.26],[0.24,0.26],[0.23,0.25],[0,50],[0,50],[0,50],[0,50]]) # boundary

count = 0
#----------------------------------------------------------------------
creator.create("FitnessMax", base.Fitness, weights=(-1.0,))
creator.create("Individual", numpy.ndarray, fitness=creator.FitnessMax)

toolbox = base.Toolbox()

n_gene = 13 # number of parameters
min_ind = numpy.ones(n_gene) * -1.0
max_ind = numpy.ones(n_gene) *  1.0
for i in range(n_gene):
  min_ind[i] = b1[i][0]
  max_ind[i] = b1[i][1]
  #print min_ind[i], max_ind[i]
#----------------------------------------------------------------------
def create_ind_uniform(min_ind, max_ind):
  ind = []
  for min, max in zip(min_ind, max_ind):
    ind.append(random.uniform(min, max))
  return ind
#----------------------------------------------------------------------
toolbox.register("create_ind", create_ind_uniform, min_ind, max_ind)
toolbox.register("individual", tools.initIterate, creator.Individual, toolbox.create_ind)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)
#----------------------------------------------------------------------
#def evalOneMax(individual):
#  return sum(individual),
#----------------------------------------------------------------------
def evalOneMax(individual):

  print "------------------------"
  global count
  count += 1
  print count

  global no_flag
  no_flag = 0

  fi = open(file_tmp,'r')
  text = fi.read().replace('XXX0',str(individual[0]).replace("[","").replace("]",""))
  text = text.replace('XXX1',str(individual[1]).replace("[","").replace("]",""))
  text = text.replace('XXX2',str(individual[2]).replace("[","").replace("]",""))
  xc3 = (-(individual[0]*3*2/3+individual[1]) - individual[2]*4) / 4
  text = text.replace('XXX3',str(xc3).replace("[","").replace("]",""))
  text = text.replace('XXX4',str(individual[3]).replace("[","").replace("]",""))
  text = text.replace('XXX5',str(individual[4]).replace("[","").replace("]",""))
  #text = text.replace('XXX6',str(individual[5]).replace("[","").replace("]",""))
  text = text.replace('XXX6','5420.3')
  text = text.replace('XXX7',str(individual[6]).replace("[","").replace("]",""))
  text = text.replace('XXX8',str(individual[7]).replace("[","").replace("]",""))
  text = text.replace('XXX9',str(individual[8]).replace("[","").replace("]",""))
  text = text.replace('XX10',str(individual[9]).replace("[","").replace("]",""))
  text = text.replace('XX11',str(individual[10]).replace("[","").replace("]",""))
  text = text.replace('XX12',str(individual[11]).replace("[","").replace("]",""))
  #text = text.replace('XX12','0.0')
  text = text.replace('XX13',str(individual[12]).replace("[","").replace("]",""))
  #text = text.replace('XX13','845.8')
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
    y = 99999999.9
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
      commands.getoutput("cp -p "+file_inp+" ./"+output_file_name+"/"+file_inp+"."+str(count))
      commands.getoutput("cp -p "+file_out+" ./"+output_file_name+"/"+file_out+"."+str(count))
    
  print y, individual
  print "------------------------"

  return y,
#----------------------------------------------------------------------
def cxTwoPointCopy(ind1, ind2):
  size = len(ind1)
  cxpoint1 = random.randint(1, size)
  cxpoint2 = random.randint(1, size-1)
  if (cxpoint2 >= cxpoint1):
    cxpoint2 += 1
  else:
    cxpoint1, cxpoint2 = cxpoint2, cxpoint1

  ind1[cxpoint1:cxpoint2], ind2[cxpoint2:cxpoint2] = ind2[cxpoint1:cxpoint2].copy(), ind1[cxpoint1:cxpoint2].copy()

  return ind1, ind2
#----------------------------------------------------------------------
def mutUniformDbl(individual, min_ind, max_ind, indpb):
  size = len(individual)
  for i, min, max in zip(xrange(size), min_ind, max_ind):
    if (random.random() < indpb):
      individual[i] = random.uniform(min, max)
  return indivisual,
#----------------------------------------------------------------------
toolbox.register("evaluate", evalOneMax)
toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)
toolbox.register("select", tools.selTournament, tournsize=3)
#----------------------------------------------------------------------
def main():
  random.seed(64)
  pop = toolbox.population(n=300)
  hof = tools.HallOfFame(1, similar=numpy.array_equal)
  stats = tools.Statistics(lambda ind: ind.fitness.values)
  stats.register("avg", numpy.mean)
  stats.register("std", numpy.std)
  stats.register("min", numpy.min)
  stats.register("max", numpy.max)
  algorithms.eaSimple(pop, toolbox, cxpb=0.5, mutpb=0.2, ngen=50000, stats=stats, halloffame=hof)
  return pop, stats, hof
#----------------------------------------------------------------------
if (__name__ == "__main__"):
  main()
#----------------------------------------------------------------------

