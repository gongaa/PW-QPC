#ifndef SIMULATION_HPP_
#define SIMULATION_HPP_
#include <iostream>
#include <random>
#include <cassert>
#include <bits/stdc++.h>
#include "Test/Test_RM.hpp"
#include "Encoder/Encoder_RM.hpp"
#include "Decoder/Decoder_RM_SC.hpp"
#include "Decoder/Decoder_RM_SCL.hpp"
#include "Decoder/Decoder_RM_syndrome_SC.hpp"
#include "Channel/Channel.hpp"

int simulation_RM_SCL();

int simulation_RM_CSS(int m, int rx, int rz, int list_size);

int simulation_RM_degeneracy(int m, int rx, int rz, double p);

#endif