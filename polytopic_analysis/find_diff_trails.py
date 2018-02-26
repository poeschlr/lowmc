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
SBox = [0x00, 0x01, 0x03, 0x06, 0x07, 0x04, 0x05, 0x02]
invSBox = [0x00, 0x01, 0x07, 0x02, 0x05, 0x06, 0x03, 0x04]
DDIFF_SIZE = 3
REFERENCE_PT = 0

class NotEnoughDegreesOfFreedom(Exception):
    def __init__(self, required_degrees, degrees):
        self.required_degrees = required_degrees
        self.degrees = degrees

    def __str__(self):
        return "Posibility space not enough degrees of freedom. Would need at least {:d} has {:d}".format(self.required_degrees, self.degrees)

def init_diff_propagation(table):
    result = {}
    for i in range(len(table)):
        idx = to_gf2_vector(i, SBOX_SIZE).row()[0]
        result[idx] = [to_gf2_vector(o, SBOX_SIZE) for o in table[i]]
    return result

verbosity_level = 0

try:
    basestring  # attempt to evaluate basestring
    def isstr(s):
        return isinstance(s, basestring)
except NameError:
    def isstr(s):
        return isinstance(s, str)

def to_gf2_vector(input, size):
    if type(input) is int:
        return vector(GF(2), [x for x in list('{0:0{width}b}'.format(input, width=size))])
    elif isstr(input):
        return to_gf2_vector(int(input, 0), size)
    else:
        raise TypeError('can only convert int or correctly formated string to gf2 vector, got {} with type {}'.format(input, type(input)))

def print_list(lst, required_verbosity = 0):
    if verbosity_level >= required_verbosity:
        ml = 0
        for l in lst:
            s = str(l)
            if len(s) > ml:
                ml = len(s)
            print(l, end='\n')
        print('-'*ml, end='\n\n')

def debug_output(msg, required_verbosity = 0, end='\n'):
    if verbosity_level >= required_verbosity:
        print(msg, end=end)

class LowMC():

    def __init__(self, lowmc_instance_description):
        self.blocksize = lowmc_instance_description['settings']['blocksize']
        self.keysize = lowmc_instance_description['settings']['keysize']
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

        self.set_key(lowmc_instance_description['settings']['key'])

    def set_key(self, key):
        debug_output("Setup roundkeys from key ", 1, end="")
        key = to_gf2_vector(key, self.keysize)
        self.round_keys = []
        for r in range(self.rounds + 1):
            self.round_keys.append(self.key_matrixes[r]*key)

        zv = VectorSpace(GF(2), self.blocksize).zero()
        self.reduced_round_keys = [copy(zv) for i in range(self.rounds + 1)]
        prev_rk_rest = copy(zv)

        sbox_bits = self.num_sboxes*SBOX_SIZE
        first_sbox_bit_idx = self.blocksize - sbox_bits

        for r in range(self.rounds):
            cr = self.rounds - r
            temp_rk = prev_rk_rest + self.round_keys[cr]
            temp_rk = self.inv_affine_matrixes[cr - 1]*temp_rk


            for i in range(self.blocksize):
                if i < first_sbox_bit_idx:
                    prev_rk_rest[i] = temp_rk[i]
                else:
                    self.reduced_round_keys[cr][i] = temp_rk[i]


        self.reduced_round_keys[0] = prev_rk_rest + self.round_keys[0]

        debug_output("[   Done   ]", 1)

    def substitution(self, input, inverse=False):
        sbox_bits = self.num_sboxes*SBOX_SIZE
        ct = copy(input)
        for i in range(sbox_bits):
            ct[self.blocksize - i - 1]=0
        for sbox_idx in range(self.num_sboxes):
            target_end = self.blocksize - SBOX_SIZE*(sbox_idx)
            sb_in = input[target_end - SBOX_SIZE : target_end]
            if inverse:
                sb_out = to_gf2_vector(invSBox[int(''.join(map(str, sb_in)), 2)]<<(sbox_idx*SBOX_SIZE), self.blocksize)
            else:
                sb_out = to_gf2_vector(SBox[int(''.join(map(str, sb_in)), 2)]<<(sbox_idx*SBOX_SIZE), self.blocksize)
            ct += sb_out

        return ct


    def encrypt(self, input, rounds=None):
        rd = rounds if rounds is not None else self.rounds
        ct = input + self.round_keys[0]
        for r in range(rd):
            debug_output('round {:d}'.format(r+1), 2)
            debug_output(ct, 2)
            ct = self.substitution(ct)
            debug_output(ct, 2)
            ct = self.affine_matrixes[r]*ct
            debug_output(ct, 2)
            ct += self.round_constants[r]
            debug_output(ct, 2)
            ct += self.round_keys[r+1]
            debug_output(ct, 2)

        return ct

    def encrypt_reduced(self, input, rounds=None):
        rd = rounds if rounds is not None else self.rounds
        ct = input + self.reduced_round_keys[0]
        for r in range(rd):
            debug_output('round {:d}'.format(r+1), 2)
            debug_output(ct, 2)
            ct = self.substitution(ct)
            debug_output(ct, 2)
            ct += self.reduced_round_keys[r+1]
            debug_output(ct, 2)
            ct = self.affine_matrixes[r]*ct
            debug_output(ct, 2)
            ct += self.round_constants[r]
            debug_output(ct, 2)

        return ct

    def decrypt(self, input, rounds=None):
        rd = rounds if rounds is not None else self.rounds
        pt = copy(input)
        for rn in range(rd):
            r = self.rounds-rn
            pt += self.round_keys[r]
            pt += self.round_constants[r-1]
            pt =  self.inv_affine_matrixes[r-1]*pt
            pt =  self.substitution(pt, inverse=True);
        return pt + self.round_keys[0]

    def decrypt_reduced(self, input, rounds=None):
        rd = rounds if rounds is not None else self.rounds
        pt = copy(input)
        for rn in range(rd):
            r = self.rounds-rn
            pt += self.round_constants[r-1]
            pt =  self.inv_affine_matrixes[r-1]*pt
            pt += self.reduced_round_keys[r]
            pt =  self.substitution(pt, inverse=True);
        return pt + self.reduced_round_keys[0]

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

            v = to_gf2_vector(i, self.sorted_space.nrows())
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


    def _permutation_diff(self, in_diff, resulting_diffs, sbox_idx):
        if sbox_idx == len(resulting_diffs)-1:
            return [self.set_bits_at_sbox_out(sbox_idx, in_diff, o) for o in resulting_diffs[sbox_idx]]

        result = []

        lower_res = self._permutation_diff(in_diff, resulting_diffs, sbox_idx + 1)
        for o in resulting_diffs[sbox_idx]:
            for lr in lower_res:
                result.append(self.set_bits_at_sbox_out(sbox_idx, lr, o))
        return result

    def propagate_diff(self, _in_diff, round, direction = 'F'):
        if direction != 'F':
            in_diff = self.inv_affine_matrixes[round]*_in_diff
        else:
            in_diff = _in_diff


        #print(type(in_diff))

        od = []
        for sbox_idx in range(self.num_sboxes):
            sbox_diff_bits = in_diff[-SBOX_SIZE*(sbox_idx+1): None if sbox_idx == 0 else -SBOX_SIZE*sbox_idx].row()[0]
            if direction == 'F':
                x = self.diff_propagation_forward[sbox_diff_bits]
            else:
                x = self.diff_propagation_backward[sbox_diff_bits]
            od.append(x)

        result = self._permutation_diff(in_diff, od, 0)
        if direction == 'F':
            return [self.affine_matrixes[round]*v for v in result]
        return result

    def _permutation_ddiff(self, diffs, ddiff_idx):
        if ddiff_idx == 0:
            return [[d] for d in diffs[0]]

        result = []

        lower_res = self._permutation_ddiff(diffs, ddiff_idx-1)

        for lr in lower_res:
            for d in diffs[ddiff_idx]:
                result.append(lr + [d])
        return result

    def propagate_ddiff(self, in_ddiff, round, direction = "F"):
        #print_list(in_ddiff)
        #print("------------------------------")
        out_diffs = []
        for in_diff in in_ddiff:
            #print(self.propagate_diff_forward(in_diff, round))
            out_diffs.append(self.propagate_diff(in_diff, round, direction))

        return self._permutation_ddiff(out_diffs, len(out_diffs)-1)
        #print_list(out_diffs)

    def propagate_ddiff_forward_till_round(self, in_ddiff, round):
        ddiffs_current_round = self.propagate_ddiff(in_ddiff, 0, 'F')

        ddiffs_after_round = ddiffs_current_round # make it work for trails of lenght 1
        for r in range(1, round):
            ddiffs_after_round = []
            for ddiff in ddiffs_current_round:
                ddiffs_after_round += self.propagate_ddiff(ddiff, r, 'F')
                #ToDo remove dublicates.
            print('ddiffs after round {}: {}'.format(r, len(ddiffs_after_round)))
            ddiffs_current_round = ddiffs_after_round
            #print("Num ddiffs after round {} is {}".format(r, len(ddiffs_current_round)))

        return ddiffs_after_round

    def propagate_ddiff_backward_from_to_round(self, in_ddiff, from_round, to_round):
        ddiffs_current_round = [in_ddiff]#self.propagate_ddiff(in_ddiff, from_round, 'B')

        ddiffs_after_round = ddiffs_current_round # make it work for trails of lenght 1
        for r in range(from_round - 1, to_round-1, -1):
            ddiffs_after_round = []
            for ddiff in ddiffs_current_round:
                ddiffs_after_round += self.propagate_ddiff(ddiff, r, 'B')
                #ToDo remove dublicates.

            ddiffs_current_round = ddiffs_after_round
            print('ddiffs before round {}: {}'.format(r, len(ddiffs_after_round)))

        return ddiffs_after_round

def check_collision(ddiff_1, ddiff_2):
    for d in ddiff_1:
        if d in ddiff_2:
            return True

    return False

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate the constands for a LowMC instance.')
    parser.add_argument('-d', '--definition', type=str, nargs=1)
    parser.add_argument('-s', '--num_sboxes', type=int, nargs='?', default=1)
    parser.add_argument('-k', '--key', type=str, nargs='?', default=1)
    parser.add_argument('-v', '--verbose', action='count')
    args = parser.parse_args()

    with open(args.definition[0], 'r') as config_stream:
        try:
            lowmc_instance = yaml.load(config_stream)
        except yaml.YAMLError as exc:
            print(exc)
            exit()

    lowmc_instance['settings']['num_sboxes'] = args.num_sboxes
    lowmc_instance['settings']['key'] = args.key


    lowmc = LowMC(lowmc_instance)

    lowmc.getInputForGoodTrail()
    try:
        ddiff_pt_side = lowmc.get_optimal_ddiff_of_len(DDIFF_SIZE)
    except Exception as e:
        print(e)
        exit()

    pt = [to_gf2_vector(REFERENCE_PT, lowmc.blocksize)]
    for d in ddiff_pt_side:
        pt.append(pt[0]+d)

    print_list(ddiff_pt_side)
    print_list(pt)

    ct = [lowmc.encrypt(p) for p in pt]
    print_list(ct)
    ct[0] += to_gf2_vector(7,lowmc.blocksize)

    ddiff_ct_side = [ct[0] + ct[i] for i in range(1,len(ct))]
    print(type(ddiff_pt_side[0]))
    print(type(ddiff_ct_side[0]))

    round_mid = ceil(lowmc.rounds/2)+1
    ddiff_m_f = lowmc.propagate_ddiff_forward_till_round(ddiff_pt_side,round_mid)
    print('forward probagation resulted in {} ddiffs'.format(len(ddiff_m_f)))

    ddiff_m_b = lowmc.propagate_ddiff_backward_from_to_round(ddiff_ct_side,lowmc.rounds,round_mid)
    print('backward probagation resulted in {} ddiffs'.format(len(ddiff_m_b)))

    print(check_collision(ddiff_m_b, ddiff_m_f))
    # x = vector(GF(2), [0,0,1,0,0,1,1,0])
    # print(lowmc.affine_matrixes[0])
    # print('\n------------------------------\n')
    # print(lowmc.affine_matrixes[0]*x)
