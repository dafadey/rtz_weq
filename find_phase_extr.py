#!/usr/bin/python3
import threading
import os
import math
import numpy

#Kerr = '0.0377' 
Kerr = '0.00754'
#Kerr = '0.017'
airDensity = '0.0122'
ampW2 = '0.00005'
default_prefix = 'NEW'
default_prefix_final = 'NEW_FINAL'
prefix = ''
skips = 'skipPics=1 skip1D=1 skipSources=1 dynamic_save_interval=0 saveW3envelope=0'
params=['Kerr='+Kerr,'airDensity='+airDensity] 
draw_interval = 'draw_interval=0'
#amps=['0.00005','0.000075','0.0001','0.00013','0.00017','0.0002','0.00025','0.0003']
#amps=['0.00015','0.000185','0.000225']

amps=['0.00005','0.000075','0.0001','0.00013', '0.00015','0.00017', '0.000185', '0.0002', '0.000225', '0.00025', '0.0003']

#amps=['0.0002','0.0004','0.0006','0.0008','0.001','0.0012','0.0014','0.0016', '0.0018', '0.002']

workers_count = 8

def cmd(amp, phase):
  cmd2run = './harm prefix="'+prefix+'" ampW0='+amp+' ampW2='+ampW2+' phase2w='+str(phase)
  for i in params:
    cmd2run = cmd2run + ' ' + i
  cmd2run = cmd2run +' '+draw_interval+' '+skips
  return cmd2run	

def datapath(amp, phase):
  path = prefix+'_harm_ampW0='+str(amp)+'_ampW2='+ampW2+'_phase2w='+str(phase)+'_'
  for i in params:
    path = path + i + '_'
  return path

w3FilePostfix = '/pics/res/field_energy0030969.dat'

def getW3(path):
  w3=0
  with open(path+w3FilePostfix , 'r') as file:
    line_id=0
    for line in file:
      word_id=0
      for word in line.split():
        if line_id == 4 and word_id == 2:
          w3 = word
        word_id+=1
      line_id+=1
  return float(w3)

class task:
  def __init__(self, cmd, path, amp, phase, w3, minmax='none'):
    self.cmd = cmd
    self.path = path
    self.amp = amp
    self.phase = phase
    self.w3 = w3
    self.minmax = minmax
  def __str__(self):
    return 'cmd='+self.cmd+', path='+self.path+', amp='+str(self.amp)+' phase='+str(self.phase)+' w3='+str(self.w3)

def run(arg):
    #print('running: ' + arg)
    os.system(arg)
    
def run_pool(pool):
  threads = []
  n = min(workers_count, len(pool));
  for i in range(0, n):
      t = threading.Thread(name='thread_'+str(i), target=run, args=[pool[i].cmd+' GPU2use='+str(i)+' 2>/dev/null 1>/dev/null'])
      threads.append(t)

  for i in range(0, n):
      threads[i].start()

  for i in range(0, n):
      threads[i].join()

#runs all tasks whose pathes do not exist
def run_tasks(tasks):
  pool=[]
  for t in tasks:
    if os.path.exists(t.path) == False or (os.path.exists(t.path) == True and os.path.exists(t.path+w3FilePostfix) == False):
      pool.append(t)
    if len(pool) == workers_count:
      run_pool(pool)
      pool=[]
  #last chunk
  if len(pool) != 0:
    run_pool(pool)
    pool=[]

  for t in tasks:
    t.w3 = getW3(t.path)

#analysis

def balance(x, y):
  return min(x/y, y/x)

critical_balance = 0.3

def phasemix(frac0, phase0, frac1, phase1): # phase1 > phase0
  if(abs(phase1 - phase0) > .5 * math.pi):
    res=.0
    if phase1 > phase0 :
      res = (phase0+math.pi) * frac0 + phase1 * frac1
    else:    
      res = phase0 * frac0 + (phase1+math.pi) * frac1
    return res if res < math.pi else res - math.pi
  else:
    return phase0*frac0+phase1*frac1

def gentasks(tasks):  
  phasew3 = []
  newtasks = []
  for i in range(0, len(tasks)):
    phasew3.append([tasks[i].phase, tasks[i].w3])
    if i+1 == len(tasks) or tasks[i+1].amp != tasks[i].amp: #add new tasks
      amp = tasks[i].amp
      maxw3=0
      minw3=1e17
      maxw3i=0
      minw3i=0
      for j in range(0, len(phasew3)):
        if phasew3[j][1] > maxw3:
          maxw3 = phasew3[j][1]
          maxw3i = j
        if phasew3[j][1] < minw3:
          minw3 = phasew3[j][1]
          minw3i = j
      n = len(phasew3)
      rminphasew3 = phasew3[(minw3i+1)%n]
      lminphasew3 = phasew3[minw3i-1]
      lmindiff = lminphasew3[1] - phasew3[minw3i][1]
      rmindiff = rminphasew3[1] - phasew3[minw3i][1]
      rmaxphasew3 = phasew3[(maxw3i+1)%n]
      lmaxphasew3 = phasew3[maxw3i-1]
      lmaxdiff = - lmaxphasew3[1] + phasew3[maxw3i][1]
      rmaxdiff = - rmaxphasew3[1] + phasew3[maxw3i][1]

      newphases = []
      
      if lmindiff < .3 * rmindiff:
        newphases.append(phasemix(.75,lminphasew3[0],.25,phasew3[minw3i][0]))
        newphases.append(phasemix(.25,lminphasew3[0],.75,phasew3[minw3i][0]))
      elif rmindiff < .3 * lmindiff:
        newphases.append(phasemix(.75,phasew3[minw3i][0],.25,rminphasew3[0]))
        newphases.append(phasemix(.25,phasew3[minw3i][0],.75,rminphasew3[0]))
      else:
        newphases.append(phasemix(.5,lminphasew3[0],.5,phasew3[minw3i][0]))
        newphases.append(phasemix(.5,rminphasew3[0],.5,phasew3[minw3i][0]))

      if lmaxdiff < .3 * rmaxdiff:
        newphases.append(phasemix(.75,lmaxphasew3[0],.25,phasew3[maxw3i][0]))
        newphases.append(phasemix(.25,lmaxphasew3[0],.75,phasew3[maxw3i][0]))
      elif rmaxdiff < .3 * lmaxdiff:
        newphases.append(phasemix(.75,phasew3[maxw3i][0],.25,rmaxphasew3[0]))
        newphases.append(phasemix(.25,phasew3[maxw3i][0],.75,rmaxphasew3[0]))
      else:
        newphases.append(phasemix(.5,lmaxphasew3[0],.5,phasew3[maxw3i][0]))
        newphases.append(phasemix(.5,rmaxphasew3[0],.5,phasew3[maxw3i][0]))
      
      for phase in newphases:
        newtasks.append(task(cmd(amp, phase), datapath(amp,phase), amp, phase, 0))
      
      phasew3 = []

  return newtasks
  
  
def inject(newtasks, tasks):
  for t in newtasks:
    amp = t.amp
    phase = t.phase
    prevphase = .0
    inside = True if amp == tasks[0].amp else False 
    found = False 
    for i in range(0, len(tasks)):
      
      if tasks[i].amp == amp:
        inside = True
        found = True
      if i+1 != len(tasks) and tasks[i].amp != amp:
        inside = False
        if found == True:
          #print('inserting '+str(amp)+' ' + str(phase) + ' at ' +str(i))
          tasks.insert(i, t)
          break
        else:
          continue
      nextphase = tasks[i].phase
      if prevphase < phase and phase < nextphase :
        #print('inserting '+str(amp)+' ' + str(phase) + ' at ' +str(i))
        tasks.insert(i, t)
        inside = False
        found = False
        break
      prevphase = tasks[i].phase
    
    if inside == True and found == True:
      #print('inserting '+str(amp)+' ' + str(phase) + ' at ' +str(len(tasks)))
      tasks.append(t)

def getOptimal(p0, w0, p1, w1, p2, w2):
  mat=numpy.array([[p0**2, p0, 1],[p1**2, p1, 1],[p2**2, p2, 1]])
  rhs=numpy.array([w0, w1, w2])
  [a,b,c] = numpy.linalg.solve(mat, rhs)
  return -b / (2*a)

def wrapPhases(p0,p1,p2):
  full = math.pi
  half = full/2
  correct = True if abs(p0 - p1) > half or abs(p0 - p2) > half or abs(p2 - p1) > half else False
  _p0 = p0
  _p1 = p1
  _p2 = p2
  if correct == True:
    _p0 = _p0 if _p0 > half else _p0+full
    _p1 = _p1 if _p1 > half else _p1+full
    _p2 = _p2 if _p2 > half else _p2+full
  return [_p0, _p1, _p2]

def genOptimalTasks(tasks):
  phasew3 = []
  newtasks = []
  for i in range(0, len(tasks)):
    phasew3.append([tasks[i].phase, tasks[i].w3])
    if i+1 == len(tasks) or tasks[i+1].amp != tasks[i].amp: #add new tasks
      amp = tasks[i].amp
      maxw3=0
      minw3=1e17
      maxw3i=0
      minw3i=0
      for j in range(0, len(phasew3)):
        if phasew3[j][1] > maxw3:
          maxw3 = phasew3[j][1]
          maxw3i = j
        if phasew3[j][1] < minw3:
          minw3 = phasew3[j][1]
          minw3i = j
      n = len(phasew3)
      rminphasew3 = phasew3[(minw3i+1)%n]
      lminphasew3 = phasew3[minw3i-1]
           
      rmaxphasew3 = phasew3[(maxw3i+1)%n]
      lmaxphasew3 = phasew3[maxw3i-1]

      [minp0,minp1,minp2] = wrapPhases(lminphasew3[0], phasew3[minw3i][0], rminphasew3[0])
      [maxp0,maxp1,maxp2] = wrapPhases(lmaxphasew3[0], phasew3[maxw3i][0], rmaxphasew3[0])
      
      #print('maximum by 3 points: '+str([maxp0,maxp1,maxp2]))
      #print('minimum by 3 points: '+str([minp0,minp1,minp2]))
      #print('--------------------')
      
      minp = getOptimal(minp0, lminphasew3[1], minp1, phasew3[minw3i][1], minp2, rminphasew3[1])
      maxp = getOptimal(maxp0, lmaxphasew3[1], maxp1, phasew3[maxw3i][1], maxp2, rmaxphasew3[1])
      
      minp = minp if minp < math.pi else minp-math.pi
      maxp = maxp if maxp < math.pi else maxp-math.pi


      newtasks.append(task(cmd(amp, minp), datapath(amp, minp), amp, minp, 0, 'min'))
      newtasks.append(task(cmd(amp, maxp), datapath(amp, maxp), amp, maxp, 0, 'max'))
      
      phasew3 = []

  return newtasks

#do the job

prefix='phases/' + default_prefix

phases_0 = []
for i in range(0, 6):
  phases_0.append(math.pi * i / 6)
  
tasks = []

for amp in amps:
  for phase in phases_0:
    tasks.append(task(cmd(amp, phase), datapath(amp, phase), amp, phase, 0))

def printtasks(tasks, inp):
  s=''
  amp=.0
  for i in range(0,len(tasks)):
    t = tasks[i]
    s+='('+str(t.phase)+','+str(t.w3)+'),'
    if i+1 == len(tasks) or tasks[i+1].amp != t.amp:
      print(inp+'.push(new pair[]{'+s+'});')
      s=''

run_tasks(tasks)

#for t in tasks:
#  print(str(t.amp)+' '+str(t.phase)+' '+str(t.w3))
#print('-----------------------')

print('pair points[][];')
print('pair newpoints[][];')
print('pair newnewpoints[][];')
print('pair finalpoints[][];')

printtasks(tasks,"points")

newtasks = gentasks(tasks)

#for t in newtasks:
#  print(str(t.amp)+' '+str(t.phase)+' '+str(t.w3))
#print('-----------------------')
run_tasks(newtasks)

inject(newtasks, tasks)

#for t in tasks:
#  print(str(t.amp)+' '+str(t.phase)+' '+str(t.w3))
#print('-----------------------')
printtasks(tasks,'newpoints')

newtasks = gentasks(tasks)

run_tasks(newtasks)

inject(newtasks, tasks)

printtasks(tasks,'newnewpoints')

skips='skipPics=0 skip1D=0 skipSources=0 dynamic_save_interval=1 saveW3envelope=1 saveFieldEnvelope=1'
prefix='phases/' + default_prefix_final
finaltasks = genOptimalTasks(tasks)
run_tasks(finaltasks)

printtasks(finaltasks,'finalpoints')

for t in finaltasks:
  prefix='phases/' + default_prefix_final
  inpath = datapath(t.amp, t.phase)
  prefix=default_prefix
  topath = datapath(t.amp, t.minmax)
  cmd = 'ln -s ' + inpath + ' ' + topath
  print(cmd)
  run(cmd)

#print("all done")
