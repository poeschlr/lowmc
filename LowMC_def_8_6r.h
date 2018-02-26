#ifndef _LOWMC_EXTERNAL_DEFINITION_H_
#define _LOWMC_EXTERNAL_DEFINITION_H_

#include <string>
#include <vector>

/*
 * LowMC matrices and constants
 ******************************/

const unsigned numofboxes = 1;    // Number of Sboxes
const unsigned blocksize = 8;   // Block size in bits
const unsigned keysize = 8; // Key size in bits
const unsigned rounds = 6; // Number of rounds

const unsigned identitysize = blocksize - 3*numofboxes;

/*
 * Linear layer matrices
 ********************/
extern std::string lgen_inst_lin_layer[6][8];/*Round constants
 ********************/

extern std::string lgen_inst_round_constant[6];/*Round key matrices
 ********************/

extern std::string lgen_inst_key_matrix[7][8];

#endif //_LOWMC_EXTERNAL_DEFINITION_H_