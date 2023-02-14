import os
import subprocess
import numpy as np

N = input('Enter N: ')
l = input('Enter list size: ')
pz_min = input('Enter a range of pz, pz_min: ')
pz_max = input('Enter a range of pz, pz_max (not inclusive): ')
seed = input('Enter random seed: ')
n = input('Enter number of samples: ')
while True:
    con = input('Enter construction method (PW, HPW, RM): ')
    if con in ['PW', 'HPW', 'RM']:
        break
euler = input("Running on Euler? 'y' or 'n': ") == 'y'
if euler:
    runtime = input("Enter runtime: ")

path = 'logs_codeword_thresh'
if not os.path.exists(path):    
    try:
        os.mkdir(path)
    except OSError as error:
        print(error)   

def run_exp(s, N):
    K = str(int(int(N)/2)+1)
    dir = 'N'+N+'_K'+K+'_l'+l+'_s'+s+'_n'+n+'_'+con
    try:
        os.mkdir(path+'/'+dir)
    except OSError as error:
        print(error)  

    print("Your results will be saved under the directory " + path + "/" + dir)

    for pz in np.arange(float(pz_min), float(pz_max), 0.001):
        cmd = "./build/apps/program -N "+N+" -Kz "+K+" -Kx "+K+" -l "+l+" -seed "+s+" -n "+n+" -con "+con+" -pz "+str(round(pz,3))
        dest = path + "/" + dir + "/pz" + str(int(round(pz,3)*1000)) + ".log"
        if euler:
            process = subprocess.Popen(['bsub', '-W', runtime+':00', cmd+' &> '+dest])
        else:
            process = subprocess.Popen(cmd.split(), stderr=open(dest, 'w'))

for s in seed.split(','):
    for NN in N.split(','):
        run_exp(s, NN)