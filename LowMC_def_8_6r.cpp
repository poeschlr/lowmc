#include "LowMC_def_8_6r.h"
/*
 * LowMC matrices and constants
 ******************************/

/*
 * Linear layer matrices
 ********************/
std::string lgen_inst_lin_layer[6][8] = {
  { // Linear Layer - Round 0
    "11111010",
    "10010001",
    "01111101",
    "00100111",
    "00000010",
    "11100110",
    "00000011",
    "00101011"
  },
  { // Linear Layer - Round 1
    "01100010",
    "11011110",
    "01100000",
    "10001100",
    "10111000",
    "11001100",
    "00111011",
    "01010101"
  },
  { // Linear Layer - Round 2
    "11111010",
    "01000010",
    "01101010",
    "10100011",
    "00100010",
    "11100110",
    "10110100",
    "11011010"
  },
  { // Linear Layer - Round 3
    "11000001",
    "11100011",
    "01001010",
    "00110101",
    "00000100",
    "01000100",
    "10000100",
    "11001000"
  },
  { // Linear Layer - Round 4
    "01000101",
    "11100010",
    "00100011",
    "10001010",
    "11001001",
    "10000001",
    "01000111",
    "00111111"
  },
  { // Linear Layer - Round 5
    "11011100",
    "01110000",
    "10111011",
    "01111101",
    "11110000",
    "11100100",
    "10100111",
    "01101000"
  }
};

/*Round constants
 ********************/

std::string lgen_inst_round_constant[6] = {
  "10101000", // Round 0
  "11101010", // Round 1
  "01111100", // Round 2
  "11101100", // Round 3
  "11101100", // Round 4
  "01111111"  // Round 5
};

/*Round key matrices
 ********************/

std::string lgen_inst_key_matrix[7][8] = {
  { // Round key matrix - Round 0
    "01010010",
    "10010001",
    "00000010",
    "10101111",
    "01100111",
    "10101100",
    "00011000",
    "10010111"
  },
  { // Round key matrix - Round 1
    "10010011",
    "00001110",
    "10010100",
    "00000101",
    "11011101",
    "11100001",
    "11100010",
    "01011101"
  },
  { // Round key matrix - Round 2
    "00001110",
    "11111011",
    "00000101",
    "10110101",
    "10000110",
    "00001111",
    "11011110",
    "11101011"
  },
  { // Round key matrix - Round 3
    "01110101",
    "10110101",
    "01110110",
    "11010001",
    "01110000",
    "10000111",
    "01011111",
    "00101000"
  },
  { // Round key matrix - Round 4
    "10011111",
    "00101101",
    "11000011",
    "00110011",
    "10010011",
    "10110101",
    "10010010",
    "01110101"
  },
  { // Round key matrix - Round 5
    "00010101",
    "11100101",
    "00010001",
    "00011000",
    "11101010",
    "01111011",
    "11010110",
    "00111011"
  },
  { // Round key matrix - Round 6
    "11111000",
    "10000010",
    "01110000",
    "00001001",
    "10111010",
    "10100110",
    "00010111",
    "11110100"
  }
};

