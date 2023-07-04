import numpy as np
# This script can facilitate the optimization of beta.
# For a given N and K, first adjust beta to meet your desired distance.
# Let n = log(N), n_half = n // 2.
# Distance can be at most 2 ** n_half, this distance can always be achieved
# for K <= {n \choose n_half}, using Reed-Muller (beta = 1).
# For why the following determines the distance of the code, 
# see the comments in src/Simulation/Simulation_polar_codeword.cpp.
beta = 2**(0.25) - 0.02
print("beta =", beta)
N = 1024
K = 42  # K = Kx + Kz - N
Kx = (N+K) // 2
Kz = Kx # symmetric code
print("N = {}, Kx = {}, Kz = {}".format(N, Kx, Kz))
n = int(np.log2(N))
pow = [beta**i for i in range(0, n)]
def cal_pw(i):
    b = bin(i)
    j = 0
    weight = 0.0
    for c in reversed(b[2:]):
        if int(c) == 1:
            weight += pow[j]
        j += 1
    return weight

pw_arr = [cal_pw(i) for i in range(N)]
pw_rank_arr = np.argsort(pw_arr)
# the distance is at least dx or dz
# log_dz = min([bin(int(c)).count("1") for c in pw_rank_arr[N-Kz:]])   # log(dz) of Cz = [N, Kz, dz] (X-type codeword)
# log_dx = n - max([bin(int(c)).count("1") for c in pw_rank_arr[:Kx]]) # log(dx) of Cx = [N, Kx, dx] (Z-type codeword)
# but the quantum code distance is actually the minimum weight of all the coset code (coset leader is some logical operator, add the stabilizers),
# this is not a linear code, and this minimum weight can be larger than dx and dz
# I will give a proof for the distance written in the table of the paper in the new version.
a = pw_rank_arr[N-Kz:Kx]
X_distance = min([bin(int(c)).count("1") for c in a])
Z_distance = n - max([bin(int(c)).count("1") for c in a])
print("X distance:", 2 ** X_distance)
print("Z distance:", 2 ** Z_distance)
print("distance:", 2 ** min(X_distance, Z_distance)) # d = min(dx, dz) for this CSS-constructed [[N, K, d]] QPC