#!/usr/bin/env python

# Sage needs python 2
# Set SAGE_LOCAL to the same location as used if you rund the sage console
# (type SAGE_LOCAL into console to find it out.)
# example:
# export SAGE_LOCAL='/usr/lib64/sagemath/local'

from sage.all import *
import argparse
import yaml

sbox_size = 3

def getInputForGoodTrail(lowmc_inst):
    affine_matrixes = [matrix(GF(2), M) for M in lowmc_inst['linear_layers']]
    blocksize = lowmc_inst['settings']['blocksize']
    num_sboxes = lowmc_inst['settings']['num_sboxes']
    rounds_with_prop1 = ceil(blocksize/(num_sboxes*sbox_size))-1
    current_affine_trail = affine_matrixes[0]

    round = 0

    sb_activation = current_affine_trail[-sbox_size*num_sboxes:]

    while rank(sb_activation) < blocksize-1:
        round += 1
        current_affine_trail = affine_matrixes[round]*current_affine_trail
        new_equations = current_affine_trail[-sbox_size*num_sboxes:].rows()

        if round < rounds_with_prop1:
            sb_activation = matrix(GF(2), sb_activation.rows() + new_equations)
        else:
            if round == rounds_with_prop1:
                posibility_space = sb_activation.right_kernel()

            idx = sbox_size*num_sboxes
            while idx > 0 and rank(sb_activation) < blocksize:
                sb_activation = matrix(GF(2), sb_activation.rows()+[new_equations[-idx]])
                idx -= 1

    print(posibility_space)
    print('\n-----------------\n')

    print(sb_activation)


    print('\n-----------------\nInput diff with longest trail:')
    print(sb_activation.right_kernel())

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate the constands for a LowMC instance.')
    parser.add_argument('-d','--definition', type=str, nargs=1, default='matrices_and_constants.yaml')
    parser.add_argument('-s', '--num_sboxes', type=int, nargs='?', default=1)
    args = parser.parse_args()


    with open(args.definition[0], 'r') as config_stream:
        try:
            lowmc_instance = yaml.load(config_stream)
        except yaml.YAMLError as exc:
            print(exc)
            exit()

    lowmc_instance['settings']['num_sboxes'] = args.num_sboxes

    getInputForGoodTrail(lowmc_instance)
