#!/usr/bin/env python

# Sage needs python 2
# Set SAGE_LOCAL to the same location as used if you rund the sage console
# (type SAGE_LOCAL into console to find it out.)
# example:
# export SAGE_LOCAL='/usr/lib64/sagemath/local'

from __future__ import print_function
#import Exception

from ddt import *
from sage.all import *
import argparse
import yaml

SBOX_SIZE = 3

class NotEnoughDegreesOfFreedom(Exception):
    def __init__(required_degrees, degrees):
        self.required_degrees = required_degrees
        self.degrees = degrees

    def __str__():
        return "Posibility space not enough degrees of freedom. Would need {:d} has {:d}".format(self.required_degrees, self.degrees)

class LowMC():
    def __init__(self, lowmc_instance_description):
        self.blocksize = lowmc_instance_description['settings']['blocksize']
        self.num_sboxes = lowmc_instance_description['settings']['num_sboxes']
        self.rounds = lowmc_instance_description['settings']['rounds']
        self.rounds_with_prop1 = ceil(self.blocksize/(self.num_sboxes*SBOX_SIZE))-1
        print("Create matrixes from yaml data ", end="")
        self.affine_matrixes = [matrix(GF(2), M) for M in lowmc_instance_description['linear_layers']]
        self.inv_affine_matrixes = [M**-1 for M in self.affine_matrixes]
        self.key_matrixes = [matrix(GF(2), M) for M in lowmc_instance_description['roundkey_matrices']]
        self.round_constants = [vector(GF(2), v) for v in lowmc_instance_description['round_constants']]
        print("[   Done   ]")

def getInputForGoodTrail(lowmc_inst):
    affine_matrixes = lowmc_inst.affine_matrixes
    blocksize = lowmc_inst.blocksize
    num_sboxes = lowmc_inst.num_sboxes
    rounds_with_prop1 = lowmc_inst.rounds_with_prop1
    current_affine_trail = affine_matrixes[0]

    round = 0

    sb_activation = current_affine_trail[-SBOX_SIZE*num_sboxes:]
    posibility_space = []

    while rank(sb_activation) < blocksize-1:
        round += 1
        current_affine_trail = affine_matrixes[round]*current_affine_trail
        new_equations = current_affine_trail[-SBOX_SIZE*num_sboxes:].rows()

        if round < rounds_with_prop1:
            sb_activation = matrix(GF(2), sb_activation.rows() + new_equations)
        else:
            if round == rounds_with_prop1:
                posibility_space.append(sb_activation.right_kernel())
                old_rank = rank(sb_activation)

            idx = SBOX_SIZE*num_sboxes
            while idx > 0 and rank(sb_activation) < blocksize-1:
                sb_activation = matrix(GF(2), sb_activation.rows()+[new_equations[-idx]])
                idx -= 1
                if rank(sb_activation) > old_rank:
                    old_rank = rank(sb_activation)
                    posibility_space.append(sb_activation.right_kernel())

    return posibility_space

def int_to_gf2_vector(input, size):
    return vector(GF(2), [x for x in list('{0:0{width}b}'.format(input, width=size))])


def sort_pos_by_trail_len(posibility_space):
    sorted_space = [posibility_space[-1].basis()[0]]
    for i in range(1,len(posibility_space)):
        current_basis = posibility_space[-1-i].basis()
        for j in range(len(current_basis)):
            if current_basis[j] not in sorted_space:
                sorted_space.append(current_basis[j])

    return matrix(GF(2), list(reversed(sorted_space)))

def get_optimal_ddiff_of_len(posibility_space, size):
    max_len = 2**posibility_space[0].dimension() - 1
    if size > max_len:
        raise NotEnoughDegreesOfFreedom(max_len, size)
    sorted_space = sort_pos_by_trail_len(posibility_space)
    print(sorted_space)
    d_diff = []
    i=1
    while len(d_diff) < size:
        v = int_to_gf2_vector(i, sorted_space.nrows())
        r = v*sorted_space
        if r not in d_diff:
            d_diff.append(r)
        i += 1
        if i > 2**sorted_space.nrows()-1:
            raise NotEnoughDegreesOfFreedom(2**sorted_space.nrows()-1, i)

    return d_diff


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate the constands for a LowMC instance.')
    parser.add_argument('-d','--definition', type=str, nargs=1)
    parser.add_argument('-s', '--num_sboxes', type=int, nargs='?', default=1)
    args = parser.parse_args()

    with open(args.definition[0], 'r') as config_stream:
        try:
            lowmc_instance = yaml.load(config_stream)
        except yaml.YAMLError as exc:
            print(exc)
            exit()

    lowmc_instance['settings']['num_sboxes'] = args.num_sboxes

    lowmc = LowMC(lowmc_instance)

    posibility_space=getInputForGoodTrail(lowmc)
    print(posibility_space)
    print('--------------------------')

    #print(sort_pos_by_trail_len(posibility_space))
    print('\n'.join(['{}'.format(v) for v in get_optimal_ddiff_of_len(posibility_space, 3)]))
