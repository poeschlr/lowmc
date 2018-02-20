#!/usr/bin/env python

# Sage needs python 2
# Set SAGE_LOCAL to the same location as used if you rund the sage console
# (type SAGE_LOCAL into console to find it out.)
# example:
# export SAGE_LOCAL='/usr/lib64/sagemath/local'

from __future__ import print_function
from __future__ import division
#import Exception

from ddt import *
from sage.all import *
import argparse
import yaml

SBOX_SIZE = 3

class NotEnoughDegreesOfFreedom(Exception):
    def __init__(self, required_degrees, degrees):
        self.required_degrees = required_degrees
        self.degrees = degrees

    def __str__(self):
        return "Posibility space not enough degrees of freedom. Would need at least {:d} has {:d}".format(self.required_degrees, self.degrees)

def int_to_gf2_vector(input, size):
    return vector(GF(2), [x for x in list('{0:0{width}b}'.format(input, width=size))])

def init_diff_propagation(table):
    result = {}
    for i in range(len(table)):
        idx = int_to_gf2_vector(i, SBOX_SIZE).row()[0]
        result[idx] = [int_to_gf2_vector(o, SBOX_SIZE) for o in table[i]]
    return result

def print_list_of_vectors(list):
    for v in list:
        print(v)

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
        self.diff_propagation_forward = init_diff_propagation(possible_out_d)
        self.diff_propagation_backward = init_diff_propagation(possible_in_d)
        print("[   Done   ]")

    def getInputForGoodTrail(self):
        current_affine_trail = self.affine_matrixes[0]
        num_sbox_bits = SBOX_SIZE*self.num_sboxes
        round = 0

        sb_activation = current_affine_trail.matrix_from_rows_and_columns(
            range(self.blocksize-num_sbox_bits,self.blocksize),
            range(self.blocksize-num_sbox_bits)
            )
        self.posibility_space = []

        desired_rank = self.blocksize-num_sbox_bits-1
        while rank(sb_activation) < desired_rank:
            round += 1
            current_affine_trail = self.affine_matrixes[round]*current_affine_trail
            new_equations = []
            for eq in current_affine_trail[-num_sbox_bits:].rows():
                new_equations.append(eq[:-num_sbox_bits])

            if round < self.rounds_with_prop1-1:
                sb_activation = matrix(GF(2), sb_activation.rows() + new_equations)
            else:
                if round == self.rounds_with_prop1-1:
                    self.posibility_space.append(sb_activation.right_kernel())
                    old_rank = rank(sb_activation)
                    #print("Num guaranteed rounds={:d}, rank={:d}".format(round, old_rank))

                idx = num_sbox_bits
                while idx > 0 and rank(sb_activation) < desired_rank:
                    sb_activation = matrix(GF(2), sb_activation.rows()+[new_equations[-idx]])
                    idx -= 1
                    if rank(sb_activation) > old_rank:
                        old_rank = rank(sb_activation)
                        self.posibility_space.append(sb_activation.right_kernel())

        zero_bits = [0]*num_sbox_bits

        for i in range(len(self.posibility_space)):
            current_basis = self.posibility_space[i].basis()
            resulting_basis = []
            for j in range(len(current_basis)):
                resulting_basis.append(current_basis[j].list() + zero_bits)
            self.posibility_space[i] = resulting_basis
        return self.posibility_space


    def sort_pos_by_trail_len(self):
        sorted_space = [self.posibility_space[-1][0]]
        for i in range(1,len(self.posibility_space)):
            current_basis = self.posibility_space[-1-i]
            for j in range(len(current_basis)):
                if current_basis[j] not in sorted_space:
                    sorted_space.append(current_basis[j])

        self.sorted_space = matrix(GF(2), list(reversed(sorted_space)))

    def get_optimal_ddiff_of_len(self, size):
        max_len = 2**len(self.posibility_space[0]) - 1
        if size > max_len:
            raise NotEnoughDegreesOfFreedom(max_len, size)
        self.sort_pos_by_trail_len()
        #print("sorted")
        #print(self.sorted_space)
        #print("")
        d_diff = []
        i=0
        while len(d_diff) < size:
            i += 1
            if i > 2**self.sorted_space.nrows()-1:
                raise NotEnoughDegreesOfFreedom(i, 2**self.sorted_space.nrows()-1)

            v = int_to_gf2_vector(i, self.sorted_space.nrows())
            r = v*self.sorted_space
            if r not in d_diff:
                d_diff.append(r)

        return d_diff

    def set_bits_at_sbox_out(self, sbox_idx, replacement_source, source):
        result = copy(replacement_source)
        for bit_idx in range(SBOX_SIZE):
            target_start = self.blocksize - SBOX_SIZE*(sbox_idx + 1)
            result[target_start + bit_idx] = source[bit_idx]
        return result


    def _permutation_diff(self, in_diff, array, sbox_idx):
        if sbox_idx == len(array)-1:
            return [self.set_bits_at_sbox_out(sbox_idx, in_diff, o) for o in array[sbox_idx]]

        result = []

        lower_res = self._permutation_diff(in_diff, array, sbox_idx + 1)
        for o in array[sbox_idx]:
            for lr in lower_res:
                result.append(self.set_bits_at_sbox_out(sbox_idx, lr, o))
        return result

    def propagate_diff_forward(self, in_diff, round):
        od = []
        for sbox_idx in range(self.num_sboxes):
            id = in_diff[-SBOX_SIZE*(sbox_idx+1): None if sbox_idx == 0 else -SBOX_SIZE*sbox_idx].row()[0]
            x = self.diff_propagation_forward[id]
            od.append(x)

        return [self.affine_matrixes[round]*v for v in self._permutation_diff(in_diff, od, 0)]

    def _permutation_ddiff(self, diffs, ddiff_idx):
        if ddiff_idx == 0:
            return [[d] for d in diffs[0]]

        result = []

        lower_res = self._permutation_ddiff(diffs, ddiff_idx-1)

        for lr in lower_res:
            for d in diffs[ddiff_idx]:
                result.append(lr + [d])
        return result

    def propagate_ddiff_forward(self, in_ddiff, round):
        #print_list_of_vectors(in_ddiff)
        #print("------------------------------")
        out_diffs = []
        for in_diff in in_ddiff:
            #print(self.propagate_diff_forward(in_diff, round))
            out_diffs.append(self.propagate_diff_forward(in_diff, round))

        return self._permutation_ddiff(out_diffs,len(out_diffs)-1)
        #print_list_of_vectors(out_diffs)

    def propagate_ddiff_forward_till_round(self, in_ddiff, round):
        ddiffs_current_round = self.propagate_ddiff_forward(in_ddiff, 0)

        for r in range(1, round):
            ddiffs_after_round = []
            for ddiff in ddiffs_current_round:
                ddiffs_after_round += self.propagate_ddiff_forward(ddiff, r)
                #ToDo remove dublicates.

            ddiffs_current_round = ddiffs_after_round
            print("Num ddiffs after round {} is {}".format(r, len(ddiffs_current_round)))

        return ddiffs_after_round


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate the constands for a LowMC instance.')
    parser.add_argument('-d', '--definition', type=str, nargs=1)
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

    posibility_space=lowmc.getInputForGoodTrail()
    #print(posibility_space)
    #print('--------------------------')

    #print(sort_pos_by_trail_len(posibility_space))
    try:
        ddiff = lowmc.get_optimal_ddiff_of_len(2)
    except Exception as e:
        print(e)
        exit()
    #print('\n'.join(['{}'.format(v) for v in ddiff]))

    v=vector(ddiff[0])
    #print('--------------------------')
    x = lowmc.affine_matrixes[3]*v
    #print(x)
    #print('--------------------------')
    #print_list_of_vectors(lowmc.propagate_diff_forward(x,4))


    print_list_of_vectors(ddiff)
    print('\n\n--------------------------')

    print_list_of_vectors(lowmc.propagate_ddiff_forward_till_round(ddiff,7))
