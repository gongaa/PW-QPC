import os
import subprocess
import numpy as np
import sys

N = input('Enter N:')
while True:
    con = input('Enter construction method (PW, Q1): ')
    if con in ['PW', 'Q1']:
        break
l = input('Enter list size: ')
pz_min = input('Enter a range of pz, pz_min: ')
pz_max = input('Enter a range of pz, pz_max (not inclusive): ')
seed = input('Enter random seed (e.g. 0,24,31): ')
n = input('Enter number of samples: ')

euler = input("Running on Euler? 'y' or 'n': ") == 'y'
if euler:
    runtime = input("Enter runtime: ")

path = 'logs_Q1_PW'
if not os.path.exists(path):    
    try:
        os.mkdir(path)
    except OSError as error:
        print(error)   

Q1_dict = {'8': 4, '16': 7, '32': 8, '64': 23, '128': 16, '256': 91, 
           '512': 32, '1024': 363, '2048': 96, '4096': 1451}

if con == 'PW':
    Kz = int(int(N)/2)
    Kx = str(int(N) + 1 - Kz)
    Kz = str(Kz)
else:
    if N not in Q1_dict.keys():
        print("N not valid, abort")
        sys.exit()
    Kz = Q1_dict[N]
    Kx = str(int(N) + 1 - Kz)
    Kz = str(Kz)


def run_exp(s, Kz, Kx):
    dir = 'N'+N+'_Kz'+Kz+'_Kx'+Kx+'_l'+l+'_s'+s+'_n'+n+'_'+con
    try:
        os.mkdir(path+'/'+dir)
    except OSError as error:
        print(error)  

    print("Your results will be saved under the directory " + path + "/" + dir)

    for pz in np.arange(float(pz_min), float(pz_max), 0.01):
        cmd = "./build/apps/program -N "+N+" -Kz "+Kz+" -Kx "+Kx+" -l "+l+" -seed "+s+" -n "+n+" -con "+con+" -pz "+str(round(pz,2))
        dest = path + "/" + dir + "/pz" + str(int(round(pz,2)*100)) + ".log"
        if euler:
            process = subprocess.Popen(['bsub', '-W', runtime+':00', cmd+' &> '+dest])
        else:
            process = subprocess.Popen(cmd.split(), stderr=open(dest, 'w'))

for s in seed.split(','):
    run_exp(s, Kz, Kx)
    run_exp(s, Kx, Kz)