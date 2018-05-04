#!/usr/bin/env python

# Sage needs python 2
# Set SAGE_LOCAL to the same location as used if you rund the sage console
# (type SAGE_LOCAL into console to find it out.)
# example:
# export SAGE_LOCAL='/usr/lib64/sagemath/local'

from __future__ import print_function
from __future__ import division
#import Exception

import os
os.environ["SAGE_LOCAL"] = '/usr/lib64/sagemath/local'

import sys
sys.path.append("../")
from generate_matrices import instantiate_matrix, grain_ssg

from ddt import *
from sage.all import *
import argparse
import yaml
import traceback
import time

SBOX_SIZE = 3
SBox = mq.SBox(0x00, 0x01, 0x03, 0x06, 0x07, 0x04, 0x05, 0x02)
invSBox = mq.SBox(0x00, 0x01, 0x07, 0x02, 0x05, 0x06, 0x03, 0x04)
REFERENCE_PT = 0

num_ddiffs_after_round = []
num_ddiffs_before_round = []

def log(logline, logfile):
    print(logline)
    logfile.write(logline+'\n')


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
    elif type(input) is list and len(input) == size:
        return vector(GF(2), input)
    else:
        raise TypeError('can only convert int or correctly formated string to gf2 vector, got {} with type {}'.format(input, type(input)))

def gf2_to_int(input):
    p = len(input) - 1
    r = 0
    for x in input:
        r += int(x)*2**p
        p -= 1
    return int(r)

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

class DDiff(object):
    HASH_SIZE = 2**32
    def __init__(self, ddiff, other = None, new_history_item = None):
        if other is not None and type(other) is not DDiff:
            raise TypeError('init for DDiff expects other to be None or of type DDiff.')
        #ToDo check that ddiff is a list of gf(2) vectors.
        self.ddiff = ddiff

        self.history = []
        if other is not None:
            for other_trail in other.history:
                self.history.append(copy(other_trail))
        else:
            self.history.append([])

        if new_history_item is not None:
            for trail in self.history:
                trail.append(new_history_item)

    def merge(self, other):
        if self == other:
            self.history.extend(other.history)
        else:
            raise ValueError('Merging d-diffs (history) only allowed if they are equal.')

    def __eq__(self, other):
        if type(other) == DDiff:
            return self.ddiff == other.ddiff
        elif type(other) == list:
            return self.ddiff == other
        else:
            return False

    def __neq__(self,other):
        return not self.__eq__(other)

    def __str__(self):
        str_diffs = ''
        ml = 0
        for d in self.ddiff:
            s = str(d)
            if len(s) > ml:
                ml = len(s)
            str_diffs += s +'\n'

        str_hist = ''
        for h in self.history:
            s = str(h)
            if len(s) > ml:
                ml = len(s)
            str_hist += s +'\n'

        return '-'*ml + '\n' + str_hist + '-'*ml + '\n' + str_diffs + '-'*ml + '\n\n'

    def __getitem__(self, idx):
        return self.ddiff[idx]

    def __iter__(self):
        for d in self.ddiff:
            yield d

    def __hash__(self):
        dd = self.ddiff[0]
        for i in range(1,len(self.ddiff)):
            dd += self.ddiff[i]
        dd = gf2_to_int(dd)
        return dd % self.HASH_SIZE



class DDiffHashMap(object):
    def __init__(self):
        self.ddiffs = {}
        self.length = 0

    def append(self, ddiff):
        if type(ddiff) == list:
            ddiff_to_append = DDiff(ddiff)
        elif type(ddiff) == DDiff:
            ddiff_to_append = ddiff
        else:
            raise TypeError('DDiff hashmap only supports appending DDiffs or lists of diffs (= list of GF(2) vectors)')


        h = hash(ddiff_to_append)
        if h in self.ddiffs:
            if type(self.ddiffs[h]) is list:
                if ddiff_to_append not in self.ddiffs[h]:
                    self.length += 1
                    self.ddiffs[h].append(ddiff_to_append)
                else:
                    i = self.ddiffs[h].index(ddiff_to_append)
                    self.ddiffs[h][i].merge(ddiff_to_append)
            else:
                if self.ddiffs[h] != ddiff_to_append:
                    self.ddiffs[h] = [self.ddiffs[h], ddiff_to_append]
                    self.length += 1
                else:
                    self.ddiffs[h].merge(ddiff_to_append)
        else:
            self.ddiffs[h] = ddiff_to_append
            self.length += 1

    def __contains__(self, ddiff):
        h = hash(ddiff)
        if h in self.ddiffs:
            if type(self.ddiffs[h]) is list:
                return ddiff in self.ddiffs[h]
            return self.ddiffs[h] == ddiff
        return False

    def find(self, ddiff):
        h = hash(ddiff)
        if h in self.ddiffs:
            if type(self.ddiffs[h]) is list:
                for d in self.ddiffs[h]:
                    if d == ddiff:
                        return d
            elif self.ddiffs[h] == ddiff:
                return self.ddiffs[h]
        return None


    def __len__(self):
        return self.length

    def __str__(self):
        r = ''
        for ddiff in self:
            r += str(ddiff)
        return r

    def strDebug(self, h):
        return str(self.ddiffs[h])

    def __iter__(self):
        for h in self.ddiffs:
            if type(self.ddiffs[h]) is list:
                for ddiff in self.ddiffs[h]:
                    yield ddiff
            else:
                print('not a list')
                yield self.ddiffs[h]

class LowMC(object):
    def generate_random(self, generator_settings):
        self.blocksize = generator_settings.get('blocksize', 32)
        self.keysize = generator_settings.get('keysize', 32)
        self.num_sboxes = generator_settings.get('num_sboxes', 1)
        self.rounds = generator_settings.get('rounds', 24)
        self.rounds_with_prop1 = ceil(self.blocksize/(self.num_sboxes*SBOX_SIZE))-1

        print("Generate random matrixes ", end="")

        gen = grain_ssg()
        self.affine_matrixes = []
        self.inv_affine_matrixes = []
        for _ in range(self.rounds):
            M = matrix(GF(2), instantiate_matrix(self.blocksize, self.blocksize, gen))
            self.affine_matrixes.append(M)
            self.inv_affine_matrixes.append(M**-1)

        self.round_constants = []
        for _ in range(self.rounds):
            constant = vector(GF(2), [next(gen) for _ in range(self.blocksize)])
            self.round_constants.append(constant)

        self.key_matrixes = []
        for _ in range(self.rounds + 1):
            mat = matrix(GF(2),instantiate_matrix(self.blocksize, self.keysize, gen))
            self.key_matrixes.append(mat)

        print("[   Done   ]")

        self.diff_propagation_forward = init_diff_propagation(possible_out_d)
        self.diff_propagation_backward = init_diff_propagation(possible_in_d)

        self.possible_reduced_keys = self.init_possible_reduced_keys()

        self.set_key(generator_settings.get('key', [next(gen) for _ in range(self.keysize)]))

    def from_description_file(self, lowmc_instance_description):
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
        self.possible_reduced_keys = self.init_possible_reduced_keys()

        self.set_key(lowmc_instance_description['settings']['key'])


    def __init__(self, lowmc_instance_description=None, generator_settings=None):
        if lowmc_instance_description:
            self.from_description_file(lowmc_instance_description)
        elif generator_settings:
            self.generate_random(generator_settings)


    def init_possible_reduced_keys(self, sbox_idx=0):
        pk = []
        if sbox_idx < self.num_sboxes-1:
            pk_l = self.init_possible_reduced_keys(sbox_idx+1)

        for i in range(2**SBOX_SIZE):
            vi = to_gf2_vector(i<<(sbox_idx*SBOX_SIZE), self.blocksize)
            if sbox_idx == self.num_sboxes-1:
                pk.append(vi)
            else:
                for k in pk_l:
                    pk.append(k+vi)

        # print('idx: {}'.format(sbox_idx))
        # print_list(pk)
        return pk



    def set_key(self, key):
        debug_output("Setup roundkeys from key ", 1, end="")
        #print(key)
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
                #sb_out = to_gf2_vector(invSBox[int(''.join(map(str, sb_in)), 2)]<<(sbox_idx*SBOX_SIZE), self.blocksize)
                sb_out = vector(GF(2),
                    [0]*(self.blocksize-SBOX_SIZE*(sbox_idx+1)) +
                    invSBox(sb_in) +
                    [0]*(SBOX_SIZE*sbox_idx))
            else:
                sb_out = vector(GF(2),
                    [0]*(self.blocksize-SBOX_SIZE*(sbox_idx+1)) +
                    SBox(sb_in) +
                    [0]*(SBOX_SIZE*sbox_idx))
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


    def decrypt_round_reduced(self, input, round, round_key=None):
        pt = copy(input)
        pt += self.round_constants[round-1]
        pt =  self.inv_affine_matrixes[round-1]*pt
        if round_key is None:
            pt += self.reduced_round_keys[round]
        else:
            pt += round_key
        pt =  self.substitution(pt, inverse=True);
        return pt

    def decrypt_reduced(self, input, rounds=None):
        rd = rounds if rounds is not None else self.rounds
        pt = copy(input)
        for rn in range(rd):
            r = self.rounds-rn
            pt += self.round_constants[r-1]
            pt =  self.inv_affine_matrixes[r-1]*pt
            pt += self.reduced_round_keys[r]
            pt =  self.substitution(pt, inverse=True);
        if rd == self.rounds:
            return pt + self.reduced_round_keys[0]
        else:
            return pt

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
                old_rank = rank(sb_activation)
            else:
                if round == self.rounds_with_prop1-1:
                    self.posibility_space.append(sb_activation.right_kernel())
                    #print(self.posibility_space[-1])
                    old_rank = rank(sb_activation)
                    print("rank after rounds with propability 1: {}, expected: {}, round: {}".format(old_rank, 3*(round), round))
                    #print("Num guaranteed rounds={:d}, rank={:d}".format(round, old_rank))

                idx = num_sbox_bits
                while idx > 0 and rank(sb_activation) < desired_rank:
                    sb_activation = matrix(GF(2), sb_activation.rows()+[new_equations[-idx]])
                    idx -= 1
                    if rank(sb_activation) > old_rank:
                        old_rank = rank(sb_activation)
                        self.posibility_space.append(sb_activation.right_kernel())
            print("after round {} rank = {}".format(round,old_rank))

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

    def propagate_ddiff(self, in_ddiff, round, inverse = False):
        if inverse:
            in_ddiff_local = [self.inv_affine_matrixes[round]*d for d in in_ddiff]
        else:
            if type(in_ddiff) == DDiff:
                in_ddiff_local = in_ddiff.ddiff
            elif type(in_ddiff) == list:
                in_ddiff_local = in_ddiff
            else:
                raise TypeError('in_ddiff must be of type DDiff or list of diffs. Type is {}'.format(type(in_ddiff)))


        out_ddiffs = []

        for possible_anchor in self.possible_reduced_keys:
            out_ddiff = []

            anchor_after_sbox = self.substitution(possible_anchor, inverse)

            for in_diff in in_ddiff_local:

                out_diff = anchor_after_sbox + self.substitution(possible_anchor + in_diff, inverse)
                if not inverse:
                    out_diff = self.affine_matrixes[round] * out_diff
                out_ddiff.append(out_diff)

            out_ddiffs.append(DDiff(out_ddiff, in_ddiff, possible_anchor))

        return out_ddiffs

    def propagate_ddiff_inactive(self, in_ddiff, round, inverse = False):
        if inverse:
            return DDiff([self.inv_affine_matrixes[round]*d for d in in_ddiff])
        else:
            return DDiff([self.affine_matrixes[round]*d for d in in_ddiff])



    def propagate_ddiff_forward_till_round(self, in_ddiff, round):
        global num_ddiffs_after_round

        ddiff_current_round = in_ddiff

        ddiffs_after_round = [ddiff_current_round] # make it work for trails of lenght 1

        #print_list(ddiffs_after_round)

        for r in range(self.rounds_with_prop1):
            ddiff_current_round = self.propagate_ddiff_inactive(ddiff_current_round, r, inverse=False)
            debug_output('ddiff for inactive round {} calculated'.format(r), 1)
            num_ddiffs_after_round[r] = 1

        ddiffs_current_round = [ddiff_current_round]
        ddiffs_after_round = ddiffs_current_round


        for r in range(self.rounds_with_prop1, round):

            #print('round {}, last round? {}'.format(r, last_round))
            ddiffs_after_round = DDiffHashMap()

            for ddiff in ddiffs_current_round:
                for dd in self.propagate_ddiff(ddiff, r, inverse=False):
                    ddiffs_after_round.append(dd)
            debug_output('ddiffs after round {}: {}'.format(r, len(ddiffs_after_round)), 1)
            num_ddiffs_after_round[r] = len(ddiffs_after_round)
            ddiffs_current_round = ddiffs_after_round
            #print("Num ddiffs after round {} is {}".format(r, len(ddiffs_current_round)))

        return ddiffs_after_round

    def propagate_ddiff_backward_from_to_round(self, in_ddiff, from_round, to_round):
        global num_ddiffs_before_round
        ddiffs_current_round = [in_ddiff]#self.propagate_ddiff(in_ddiff, from_round, 'B')

        ddiffs_after_round = ddiffs_current_round # make it work for trails of lenght 1
        for r in range(from_round - 1, to_round-1, -1):
            ddiffs_after_round = DDiffHashMap()
            for ddiff in ddiffs_current_round:
                for dd in self.propagate_ddiff(ddiff, r, inverse=True):
                    ddiffs_after_round.append(dd)
                #ToDo remove dublicates.

            ddiffs_current_round = ddiffs_after_round
            debug_output('ddiffs before round {}: {}'.format(r, len(ddiffs_after_round)), 1)
            num_ddiffs_before_round[r] = len(ddiffs_after_round)
            #print('ddiffs before round {}: {}'.format(r, len(ddiffs_after_round)))

        return ddiffs_after_round

def check_collision(ddiff_forward, ddiff_backward):
    possible_trails = []
    for d_b in ddiff_backward:
        d_f = ddiff_forward.find(d_b)
        if d_f is not None:
            x = {'back': d_b.history,
                'forward': d_f.history}
            possible_trails.append(x)

    return possible_trails


def attack(lowmc, logfile, only_trail = False):
    global num_ddiffs_after_round
    num_ddiffs_after_round = [0]*lowmc.rounds
    global num_ddiffs_before_round
    num_ddiffs_before_round = [0]*lowmc.rounds

    t_start=time.time()
    lowmc.getInputForGoodTrail()
    try:
        ddiff_pt_side = lowmc.get_optimal_ddiff_of_len(args.ddiff_size)
    except NotEnoughDegreesOfFreedom as e:
        print(e)
        print_list(lowmc.posibility_space)
        exit()
    else:
        traceback.print_exc()

    if only_trail:
        print_list(lowmc.posibility_space)
        print_list(ddiff_pt_side)
        return

    t_init=time.time()
    pt = [to_gf2_vector(REFERENCE_PT, lowmc.blocksize)]
    for d in ddiff_pt_side:
        pt.append(pt[0]+d)

    #print_list(ddiff_pt_side)
    #print_list(pt)

    ct = [lowmc.encrypt(p) for p in pt]
    #print_list(ct)

    # ddiff_ct_side = [ct[0] + ct[i] for i in range(1,len(ct))]
    #print(time.strftime("%H:%M:%S", time.localtime()))
    t_enc=time.time()
    round_mid = lowmc.rounds_with_prop1 + ceil((lowmc.rounds-lowmc.rounds_with_prop1)/2)
    ddiff_m_f = lowmc.propagate_ddiff_forward_till_round(ddiff_pt_side,round_mid)
    #print('forward probagation resulted in {} ddiffs'.format(len(ddiff_m_f)))
    #print(num_ddiffs_after_round)
    #print(ddiff_m_f)
    #print(time.strftime("%H:%M:%S", time.localtime()))
    t_forward=time.time()

    ddiff_ct_side = DDiff([ct[0] + ct[i] for i in range(1,len(ct))])
    #print(ddiff_ct_side)
    ddiff_m_b = lowmc.propagate_ddiff_backward_from_to_round(ddiff_ct_side, lowmc.rounds, round_mid)
    #print(num_ddiffs_before_round)
    #print(time.strftime("%H:%M:%S", time.localtime()))
    t_backward=time.time()

    #exit()
    #print('trails:')
    possible_trails = check_collision(ddiff_backward=ddiff_m_b, ddiff_forward=ddiff_m_f)
    #print_list(possible_trails)
    #print(time.strftime("%H:%M:%S", time.localtime()))
    t_collision=time.time()

    rb = lowmc.blocksize - 3
    round_keys = []
    for i in range(lowmc.rounds + 1):
        round_keys.append([])

    #print_list(round_keys)
    num_trails_forward = 0;
    num_trails_backward = 0;

    for collision in possible_trails:
        num_trails_forward += len(collision['forward'])
        num_trails_backward += len(collision['back'])
        for t in collision['back']:
            i = lowmc.rounds-1
            c = ct[0]
            for a in t:
                c += lowmc.round_constants[i]
                c = lowmc.inv_affine_matrixes[i]*c
                ck = vector(GF(2), [0]*rb + list((c + a)[-3:]))
                if ck not in round_keys[i+1]:
                    round_keys[i+1].append(ck)
                c += ck
                c = lowmc.substitution(c, inverse=True)
                i -= 1
    #print(time.strftime("%H:%M:%S", time.localtime()))
    t_roundkeys=time.time()

    #############################################################
    #                     Print sumary                          #
    #############################################################
    #print_list(lowmc.reduced_round_keys)
    #print_list(round_keys)

    keys_ok = True
    keys_found = 0
    for i in range(lowmc.rounds + 1):
        if len(round_keys[i]) >= 1:
            if lowmc.reduced_round_keys[i] not in round_keys[i]:
                keys_ok = False
            else:
                keys_found += 1

    log('{bs:d}, {ks:d}, {nr:d}, {dds:d}, {dd_f:s}, {dd_b:s}, '\
            '{t_find:.3f}, {t_enc:.3f}, {t_forward:.3f}, '\
            '{t_backward:.3f}, {t_collision:.3f}, {t_roundkeys:.3f}, '\
            '{num_collisions:d}, {num_trails_forward:d}, {num_trails_backward:d}, '\
            '{key_ok:s}'.format(
        bs=lowmc.blocksize, ks=lowmc.keysize, nr=lowmc.rounds, dds=3,
        dd_f=str(num_ddiffs_after_round[:round_mid]).replace(',',';'),
        dd_b=str(num_ddiffs_before_round[round_mid:]).replace(',',';'),
        t_find=t_init - t_start, t_enc=t_enc - t_init,
        t_forward=t_forward - t_enc, t_backward=t_backward - t_forward,
        t_collision=t_collision - t_backward, t_roundkeys=t_roundkeys - t_collision,
        num_collisions=len(possible_trails),
        num_trails_forward=num_trails_forward,
        num_trails_backward=num_trails_backward,
        key_ok='{:d} roundkeys found -- '.format(keys_found) + ('[  OK  ]' if keys_ok else '[ ERROR ]')
    ), logfile)





if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate the constands for a LowMC instance.')
    parser.add_argument('-d', '--definition', type=str, nargs='?', default=None)
    parser.add_argument('-a', '--auto', type=str, nargs='?', default=None)
    parser.add_argument('-s', '--num_sboxes', type=int, nargs='?', default=1)
    parser.add_argument('-k', '--key', type=str, nargs='?', default=1)
    parser.add_argument('-v', '--verbose', action='count')
    parser.add_argument('-z', '--ddiff_size', type=int, nargs='?', default=3)
    parser.add_argument('--only_trail', action='store_true')
    args = parser.parse_args()

    if args.definition:
        with open(args.definition, 'r') as config_stream:
            try:
                lowmc_instance = yaml.load(config_stream)
            except yaml.YAMLError as exc:
                print(exc)
                exit()

        lowmc_instance['settings']['num_sboxes'] = args.num_sboxes
        lowmc_instance['settings']['key'] = args.key


        lowmc = LowMC(lowmc_instance_description=lowmc_instance)
        attack(lowmc, logfile, args.only_trail)

    elif args.auto:
        with open(args.auto, 'r') as config_stream:
            try:
                auto_def = yaml.load(config_stream)
            except yaml.YAMLError as exc:
                print(exc)
                exit()
        with open('log -- {}.csv'.format(time.strftime("%Y_%m_%d - %H_%M_%S", time.localtime())), 'w') as logfile:
            log('{bs:s}, {ks:s}, {nr:s}, {dds:s}, {dd_f:s}, {dd_b:s}, '\
                    '{t_find:s}, {t_enc:s}, {t_forward:s}, '\
                    '{t_backward:s}, {t_collision:s}, {t_roundkeys:s}, '\
                    '{num_collisions:s}, {num_trails_forward:s}, {num_trails_backward:s}, '\
                    '{key_ok:s}'.format(
                bs='blocksize', ks='keysize', nr='rounds', dds='d-diff size',
                dd_f='# d-diffs foward',
                dd_b='# d-diffs backward',
                t_find='t find trail', t_enc='t encrypt',
                t_forward='t forward', t_backward='t backward',
                t_collision='t find collision', t_roundkeys='t calc roundkeys',
                num_collisions='# collisions',
                num_trails_forward='# trails forward',
                num_trails_backward='# trails backward',
                key_ok='test result'
            ), logfile)
            for definition in auto_def:
                if 'setup' in definition:
                    lowmc = LowMC(generator_settings=definition['setup'])
                else:
                    continue

                for test in definition['tests']:
                    if 'rounds' in test:
                        lowmc.rounds = test['rounds']
                    repeat = test.get('repeat', 1)
                    key = test.get('key', 'random')
                    gen = grain_ssg()
                    for i in range(repeat):
                        if key == 'random':
                            k = [next(gen) for _ in range(lowmc.keysize)]
                            print('random key {}'.format(k))
                        elif type(key) is list:
                            k = key[i]
                        else:
                            print('Key must be random or list')
                            continue
                        #print('setup new key')
                        lowmc.set_key(k)

                        attack(lowmc, logfile, args.only_trail)
                        if args.only_trail:
                            break
                    if args.only_trail:
                        break



    else:
        print('No instance definition given')
        exit()
