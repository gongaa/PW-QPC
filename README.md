# QRM-List-Decoder

This repo contains the Reed-Muller code and the Dumer's list decoder in LLR form. My implementation is based on both [ECCLab](https://github.com/kshabunov/ecclab) and [aff3ct](https://github.com/aff3ct/aff3ct).
Compared to [ECCLab](https://github.com/kshabunov/ecclab), all the codes are in C++ style and all the calculations are based on log-likelihood ratio. Furthermore, the tree metric construction and (optimized) copy follows from [aff3ct](https://github.com/aff3ct/aff3ct) but on an unbalaneced tree. Alternatively, people can use the polar code with the Reed-Muller construction. They give very similar FER performance.

This repo also contains the polar code and the SCL decoder (taken from [aff3ct](https://github.com/aff3ct/aff3ct) but deleted the multi-kernel part).
The code organization is a much simplified version of [aff3ct](https://github.com/aff3ct/aff3ct), with the modulator demodulator modules and some design patterns peeled off. 
Various polar code construction methods are included, like the Arikan's BEC, the Gaussian Approximation, the Polarization Weight, the Higher Order Polarization Weight methods. I also wrote a Reed-Muller construction so that $K$ does not have to be $m\choose \leq r$.

The goal of this project is to explore for the quantum polar code, whether the codewords in the list can be used to approximate the Maximum Likelihood decoding. Our approach is to weigh each codeword differently according to the distance to the noise (that is compatible with the syndrome), and then choose the equivalence class that has the maximum probability.
Default decoding is syndrome decoding instead of codeword decoding.

## Build and Run

To be compiled on Euler, first clone this repo and then run the following before `make`:
```Shell
env2lmod
modoule load gcc/9.3.0
```
Then build and run as:
```
mkdir logs
make
./build/apps/program -N 1024 -K 513 -l 128 -pz 0.1 -n 10000 -con PW -exact 0 -seed 25
```
On Euler, instead of the last command, run the following
```
bsub -W 24:00 "./build/apps/program -N 1024 -K 513 -l 128 -pz 0.1 -n 10000 -con PW -exact 0 -seed 25 &> logs/Polar_N1024_K513_l128_pz10_n10000_conPW_exact0_seed25_syndrome.log"
```

### Arguments
We are using independent X and Z decoding. Currently the noise is generated and decoded unilaterally.

The possible options are
* `-N` - The length of the code, use a power of two so that it is CSS.
* `-K` - The number of information bits, use $\frac{N}{2}+1$ to make the degeneracy effect the most obvious.
* `-l` - The list size. For $N=1024, K=513, l=512, n=10000$ the code runs ~48h on Euler. You can ask for more memory on Euler using `-R "rusage[mem=2048]"`.
* `-n` - The number of samples. Ideally, the number of flips of the noise should follow a binomial distribution, but for $N=1024$, $10000$ samples are not enough to approximate the distribution well.
* `-pz` - The probability of Z error happening on each bit (independently).
* `-con` - The construction method. Can choose from `BEC` `RM` `PW` `HPW` `Q1`. Note that `BEC` suffers from the numerical problem so that the construction at $N\geq 256, K\geq 129$ is not CSS anymore. `HPW` gives the same construction as `PW` for $N=256, K=129$, but not for the other $K=\frac{N}{2}+1$.
* `-exact` - When `exact` is set to $0$, the number of flips of the noise is supposed to follow a binomial distribution. When `exact` is set to a non zero integer $t$, only noise that has flips $\lfloor N\cdot p_z - t \rfloor$ or $\lceil N\cdot p_z + t \rceil$ are generated.
* `-seed` - The random seed used to generate the noise.
* `-fast` - Choose from `0` syndrome decoding that works for high rate, not memory efficient. `1` syndrome decoding that is memory efficient but does not support high rate for now ($K-\frac{N}{2}$ too large will trigger segmentation fault). `2` codeword decoding, you should only use `Q1` construction together with this version.
* `-interval` - The print interval for showing intermediate results.


