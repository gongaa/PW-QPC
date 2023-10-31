import os
import subprocess
import numpy as np
import sys

N = input('Enter N: ')
K = input('Enter K: ')
l = input('Enter list size: ')
px_min = input('Enter a range of px, px_min: ')
px_max = input('Enter a range of px, px_max (not inclusive): ')
delta_px = 0.01 # increasement
factor = 100 # for file saving
digits_to_keep = 2
# delta_px = 0.001 # increasement
# factor = 1000 # for file saving
# digits_to_keep = 3
seed = input('Enter random seeds: ')
n = input('Enter number of samples: ')
while True:
    con = input('Enter construction method (PW, HPW, RM, Q1): ')
    if con in ['PW', 'HPW', 'RM', 'Q1']:
        break

# the optimal info position at odd log(N) is subject to change
Q1_dict = {'8': 4, '16': 7, '32': 8, '64': 23, '128': 16, '256': 91, 
           '512': 32, '1024': 363, '2048': 96, '4096': 1451}

if con == 'Q1':
    if N not in Q1_dict.keys():
        print("Q1: N not valid, abort")
        sys.exit()
    Kz = Q1_dict[N]
    Kx = str(int(N) + 1 - Kz)
    Kz = str(Kz)

depolarize = input("Depolarizing channel? 'y' or 'n': ") == 'y'
euler = input("Running on Euler? 'y' or 'n': ") == 'y'
if euler:
    runtime = input("Enter runtime: ")

path = 'logs_depolarize_high_rate'
if not os.path.exists(path):    
    try:
        os.mkdir(path)
    except OSError as error:
        print(error)   

def run_exp(s, N, Kx, Kz, type='x', depolarize=0):
    if con != 'Q1':
        dir = 'N'+N+'_K'+Kx+'_l'+l+'_s'+s+'_n'+n+'_'+con
    else:
        dir = 'N'+N+'_Kx'+Kx+'_Kz'+Kz+'_l'+l+'_s'+s+'_n'+n+'_'+con
    if not os.path.exists(path+'/'+dir):
        os.mkdir(path+'/'+dir)

    print("Your results will be saved under the directory " + path + "/" + dir)

    for px in np.arange(float(px_min), float(px_max), delta_px):
        cmd = "./build/apps/program -N "+N+" -Kz "+Kz+" -Kx "+Kx+" -l "+l+" -seed "+s+" -n "+n+" -con "+con+" -px "+str(round(px,digits_to_keep))+" -dep "+depolarize+" -beta "+str(1.1127756)
        # cmd = "./build/apps/program -N "+N+" -Kz "+K+" -Kx "+K+" -l "+l+" -seed "+s+" -n "+n+" -con "+con+" -px "+str(round(px,digits_to_keep))+ " -version 1"
        dest = path + "/" + dir + "/p" + type + str(int(round(px,digits_to_keep)*factor)) + ".log"
        if euler:
            # process = subprocess.Popen(['bsub', '-W', runtime+':00', cmd+' &> '+dest])
            process = subprocess.Popen(['sbatch', '--time', runtime+':00:00', '--wrap', cmd+' &> '+dest])
        else:
            process = subprocess.Popen(cmd.split(), stderr=open(dest, 'w'))

for s in seed.split(','):
    if con != 'Q1':
        run_exp(s, N, K, K, 'x', str(int(depolarize)))
    else:
        run_exp(s, N, Kx, Kz, 'x') # in the Z basis, detect X-type noise
        run_exp(s, N, Kz, Kx, 'z') # in the X basis, detect Z-type noise