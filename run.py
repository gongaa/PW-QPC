import os
import subprocess
import numpy as np

N = input('Enter N:')
K = input('Enter K:')
l = input('Enter list size:')
pz_min = input('Enter a range of pz, pz_min: ')
pz_max = input('Enter a range of pz, pz_max (not inclusive): ')
seed = input('Enter random seed:')
n = input('Enter number of samples: ')
while True:
    con = input('Enter construction method (PW, HPW, RM): ')
    if con in ['PW', 'HPW', 'RM']:
        break
exact = input("Enter exact interval (0 to not use exact): ")
euler = input("Running on Euler? 'y' or 'n': ") == 'y'
if euler:
    runtime = input("Enter runtime: ")

path = 'logs'
if not os.path.exists(path):    
    try:
        os.mkdir(path)
    except OSError as error:
        print(error)   

def run_exp(s):
    dir = 'N'+N+'_K'+K+'_l'+l+'_s'+s+'_n'+n+'_exact'+exact+'_'+con
    try:
        os.mkdir(path+'/'+dir)
    except OSError as error:
        print(error)  

    print("Your results will be saved under the directory logs/" + dir)

    for pz in np.arange(float(pz_min), float(pz_max), 0.01):
        cmd = "./build/apps/program -N "+N+" -K "+K+" -l "+l+" -seed "+s+" -n "+n+" -exact "+exact+ " -con "+con+" -pz "+str(round(pz,2))
        dest = "logs/" + dir + "/pz" + str(int(round(pz,2)*100)) + ".log"
        if euler:
            process = subprocess.Popen(['bsub', '-W', runtime+':00', cmd+' &> '+dest])
        else:
            process = subprocess.Popen(cmd.split(), stderr=open(dest, 'w'))

for s in seed.split(','):
    run_exp(s)