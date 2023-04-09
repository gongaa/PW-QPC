import os
import subprocess
import numpy as np

N = input('Enter N: ')
K = input('Enter K: ')
l = input('Enter list size: ')
px_min = input('Enter a range of pz, pz_min: ')
px_max = input('Enter a range of pz, pz_max (not inclusive): ')
delta_px = 0.01 # increasement
factor = 100 # for file saving
digits_to_keep = 2
seed = input('Enter random seeds: ')
n = input('Enter number of samples: ')
while True:
    con = input('Enter construction method (PW, HPW, RM): ')
    if con in ['PW', 'HPW', 'RM']:
        break
euler = input("Running on Euler? 'y' or 'n': ") == 'y'
if euler:
    runtime = input("Enter runtime: ")

path = 'logs'
if not os.path.exists(path):    
    try:
        os.mkdir(path)
    except OSError as error:
        print(error)   

def run_exp(s, N, K):
    dir = 'N'+N+'_K'+K+'_l'+l+'_s'+s+'_n'+n+'_'+con
    if not os.path.exists(path+'/'+dir):
        os.mkdir(path+'/'+dir)

    print("Your results will be saved under the directory " + path + "/" + dir)

    for pz in np.arange(float(px_min), float(px_max), delta_px):
        cmd = "./build/apps/program -N "+N+" -Kz "+K+" -Kx "+K+" -l "+l+" -seed "+s+" -n "+n+" -con "+con+" -pz "+str(round(pz,digits_to_keep))
        # cmd = "./build/apps/program -N "+N+" -Kz "+K+" -Kx "+K+" -l "+l+" -seed "+s+" -n "+n+" -con "+con+" -pz "+str(round(pz,digits_to_keep))+ " -version 1"
        dest = path + "/" + dir + "/pz" + str(int(round(pz,digits_to_keep)*factor)) + ".log"
        if euler:
            process = subprocess.Popen(['bsub', '-W', runtime+':00', cmd+' &> '+dest])
        else:
            process = subprocess.Popen(cmd.split(), stderr=open(dest, 'w'))

for s in seed.split(','):
    run_exp(s, N, K)