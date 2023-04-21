# List Decoder for PW-QPC 

This repo contains the source code of the paper [Improved Logical Error Rate via List Decoding of Quantum Polar Codes](https://arxiv.org/pdf/2304.04743.pdf), where the Polarization Weight family of Quantum Polar Code is considered. 

This repo also contains some codes the classical error correction community could use. See the directory layout below for where the files are located at. Also see the commented lines in `main.cpp`.

My (classical) polar code implementation is a simplified version of [aff3ct](https://github.com/aff3ct/aff3ct). I also reimplemented the Dumer's list decoder [ECCLab](https://github.com/kshabunov/ecclab) for the Reed-Muller code. Alternatively, people can use the polar code with the Reed-Muller construction. They give very similar FER performance.

## Build and Run

Clone this repo, then build and run as:
```
make
mkdir logs
./build/apps/program -N 128 -Kx 65 -Kz 65 -l 16 -px 0.1 -n 10000 -con PW -seed 42 &> logs/Polar_N128_Kx65_Kz65_l16_pz10_n10000_conPW_seed42.log
```
Alternatively, you can run the executable by interacting with the run_SCL.py script. 

## Arguments
The independent bit- and phase-flip error model is used, hence only independent X and Z decoding is implemented. 
By default, it simulates an ($N,K_x,K_z$) PW-QPC in the Z basis, i.e., detect and correct bit-flip (X-type) errors. To simulate this code in the X basis (correct phase-flip errors), it suffices to simulate an ($N,K_z,K_x$) PW-QPC in the Z basis.

The possible options are
* `-N` - The length of the code, use a power of two.
* `-Kx`,`-Kz` - The number of information bits in X and Z basis. $K_x+K_z$ should be larger than $N$.
* `-l` - The list size. 
* `-n` - The number of samples. 
* `-px` - The probability of bit-flip error happening on each bit (independently).
* `-con` - The construction method. Can choose from `PW` `HPW` `RM` `Q1` `BEC`. Note that the quantum polar code constructed from `BEC` is not guaranteed to be CSS.
* `-seed` - The random seed used to generate the noise.
* `-version` - Choose from `0` codeword decoding (default). `1` syndrome decoding.
* `-interval` - The printing interval for showing intermediate results. Default value is $1000$.
* `-beta` - The $\beta$ used in the `PW` construction. Default value is $2^{1/4}$.

## Expected Runtime
For $N=1024,l=128$, $10^4$ samples take 4 hours on [Euler](https://scicomp.ethz.ch/wiki/FAQ). While for $N=1024,l=16$, $10^5$ samples take 2 hours.

## Directory Layout
    .
    ├── run_SCL.py                               # run ./build/apps/program 
    ├── distance.py                              # determine the distance of the code
    └── src                   
        ├── Util
        │   └── Frozen_bits_generator.cpp        # construction methods
        ├── Decoder
        │   ├── Decoder_polar.cpp                # polar code SCL decoder
        │   └── Decoder_RM_SCL.cpp               # Dumer's list decoder
        └── Simulation         
            ├── Simulation_polar_SCL.cpp         # classical polar code decoded using SCL
            ├── Simulation_RM_SCL.cpp            # classical RM code decoded using Dumer's list decoder
            ├── Simulation_polar_codeword.cpp    # QPC codeword decoder (a lot of comments here)
            └── Simulation_polar_syndrome.cpp    # QPC syndrome decoder

