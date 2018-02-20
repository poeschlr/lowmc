#include "LowMC.h"
#include <iostream>
#include <string>
#include <bitset>
#include <cmath>
#include <stdlib.h>
#include <vector>
#include <algorithm>


//////////////////
//     MAIN     //
//////////////////

bool fullAdder(bool b1, bool b2, bool& carry){
	bool sum = (b1 ^ b2) ^ carry;
	carry = (b1 && b2) || (b1 && carry) || (b2 && carry);
	return sum;
}

template<unsigned int N>
void bitsetAdd(std::bitset<N>& x, const std::bitset<N>& y){
	bool carry = false;
	for (int i = 0; i < N; i++){
		x[i] = fullAdder(x[i], y[i], carry);
	}
}

int main (int argc, char *argv[]){
    // Example usage of the LowMC class
    // Instantiate a LowMC cipher instance called cipher using the key '1'.
    LowMC cipher(1);
//    block m = 0xFFD5;
//
//    std::cout << "Plaintext:" << std::endl;
//    std::cout << m << std::endl;
//    m = cipher.encrypt( m );
//    std::cout << "Ciphertext:" << std::endl;
//    std::cout << m << std::endl;
//    m = cipher.decrypt( m );
//    std::cout << "Encryption followed by decryption of plaintext:" << std::endl;
//    std::cout << m << std::endl;

    if(argc == 3){
    	std::string str_diff = std::string(argv[1]);
    	if(str_diff.length() != blocksize){
    		std::cout << "Wrong size of input diff" << std::endl;
    		return -1;
    	}
    	block in_diff = block(str_diff);
    	int num_rounds = 0;
    	try{
    		num_rounds = std::stoi(argv[2]);
    	}
		catch (const std::invalid_argument& ia) {
			  std::cerr << "Invalid argument for number of rows: " << ia.what() << std::endl;
			  return -1;
		}

		std::vector<block> diffset;
    	for(unsigned i=0; i < std::pow(2,blocksize); i++){
    		block i1 = block(i);
    		block i2 = i1^in_diff;
    		block r1 = cipher.encrypt(i1, num_rounds);
    		block r2 = cipher.encrypt(i2, num_rounds);
    		block odiff = r1 ^ r2;

    		if(std::find(diffset.begin(), diffset.end(), odiff) == diffset.end()) {
    		    diffset.push_back(odiff);
    		}
    	}
    	std::cout << in_diff << std::endl;
    	std::cout << diffset.size() << std::endl;
    	std::vector<block>::iterator it;
    	for(it = diffset.begin(); it != diffset.end(); it++)    {
			std::cout<< *it << std::endl;
    	}
    }
    //cipher.print_matrices();
   
    return 0;
}
