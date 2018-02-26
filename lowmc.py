#!/usr/bin/env python

# Sage needs python 2
# Set SAGE_LOCAL to the same location as used if you rund the sage console
# (type SAGE_LOCAL into console to find it out.)
# example:
# export SAGE_LOCAL='/usr/lib64/sagemath/local'

from __future__ import print_function
from __future__ import division
#import Exception

from sage.all import *
import argparse
import yaml

SBOX_SIZE = 3
SBox = [0x00, 0x01, 0x03, 0x06, 0x07, 0x04, 0x05, 0x02]
invSBox = [0x00, 0x01, 0x07, 0x02, 0x05, 0x06, 0x03, 0x04]

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
    if verbosity_level > required_verbosity:
        for l in lst:
            print(l,end='\n\n')

def debug_output(msg, required_verbosity = 0, end='\n'):
    if verbosity_level >= required_verbosity:
        print(msg, end=end)

class LowMC():
    def __init__(self, lowmc_instance_description):
        self.blocksize = lowmc_instance_description['settings']['blocksize']
        self.num_sboxes = lowmc_instance_description['settings']['num_sboxes']
        self.keysize = lowmc_instance_description['settings']['keysize']
        self.rounds = lowmc_instance_description['settings']['rounds']

        debug_output("Create matrixes from yaml data ", 1, end="")
        self.affine_matrixes = [matrix(GF(2), M) for M in lowmc_instance_description['linear_layers']]
        self.inv_affine_matrixes = [M**-1 for M in self.affine_matrixes]
        self.key_matrixes = [matrix(GF(2), M) for M in lowmc_instance_description['roundkey_matrices']]
        self.round_constants = [vector(GF(2), v) for v in lowmc_instance_description['round_constants']]
        debug_output("[   Done   ]", 1)

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

def command_line():
    while True:
        print("lowmc> ", end='')
        try:
            line = sys.stdin.readline()
        except KeyboardInterrupt:
            print(" good bye ;)")
            break

        if not line:
            print(" good bye ;)")
            break

        command = line.rstrip().split(" ")

        if command[0] == 'exit':
            print(" good bye ;)")
            return

        if command[0] == 'print':
            if len(command) < 2:
                print('print command needs more input')
                continue

            if command[1] == 'keymatrix':
                if len(command) > 2:
                    print(lowmc.key_matrixes[int(command[2])])
                else:
                    print_list(lowmc.key_matrixes)
            elif command[1] == 'roundkey':
                if len(command) > 2:
                    print(lowmc.round_keys[int(command[2])])
                else:
                    print_list(lowmc.round_keys)
            elif command[1] == 'reduced_roundkey':
                if len(command) > 2:
                    print(lowmc.reduced_round_keys[int(command[2])])
                else:
                    print_list(lowmc.reduced_round_keys)
            else:
                print('unknown print command "{}"'.format(command[1]))

            continue

        if command[0] == 'sbox':
            if len(command) < 2:
                print('sbox command expects an input parameter')
                continue
            pt = int(command[1],0)
            if pt < 0 or pt > len(SBox) - 1:
                print('too many input bits for sbox command')
                continue

            print('0b{:0{width}b}'.format(SBox[pt], width=SBOX_SIZE))

        if command[0].startswith('enc'):
            if len(command) < 2:
                print('encryption command expects the plaintext as the first parameter')
                continue
            pt = command[1]
            try:
                block = to_gf2_vector(pt, lowmc.blocksize)
            except TypeError as e:
                print('Illegal input')
                debug_output(e, 2)
                continue

            rounds = int(command[2]) if len(command) > 2 else None

            ct = lowmc.encrypt(block, rounds) if 'std' in command[0] else lowmc.encrypt_reduced(block, rounds)
            print(ct)

        if command[0].startswith('dec'):
            if len(command) < 2:
                print('description command expects the plaintext as the first parameter')
                continue
            ct = command[1]
            try:
                block = to_gf2_vector(ct, lowmc.blocksize)
            except TypeError as e:
                print('Illegal input')
                debug_output(e, 2)
                continue

            rounds = int(command[2]) if len(command) > 2 else None

            pt = lowmc.decrypt(block, rounds) if 'std' in command[0] else lowmc.decrypt_reduced(block, rounds)
            #pt = lowmc.decrypt(block, rounds)
            print(pt)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate the constands for a LowMC instance.')
    parser.add_argument('-d', '--definition', type=str, nargs=1)
    parser.add_argument('-s', '--num_sboxes', type=int, nargs='?', default=1)
    parser.add_argument('-k', '--key', type=str, nargs=1)
    parser.add_argument('-v', '--verbose', action='count')
    parser.add_argument('-c', '--command_line', action='store_true')
    parser.add_argument('-p', '--plaintext', type=str, nargs='*')
    args = parser.parse_args()

    if args.verbose:
        verbosity_level = args.verbose

    with open(args.definition[0], 'r') as config_stream:
        try:
            lowmc_instance = yaml.load(config_stream)
        except yaml.YAMLError as exc:
            print(exc)
            exit()

    lowmc_instance['settings']['num_sboxes'] = args.num_sboxes
    lowmc_instance['settings']['key'] = args.key[0]

    lowmc = LowMC(lowmc_instance)

    if args.command_line:
        command_line()

    if args.plaintext:
        for pt in args.plaintext:
            p = to_gf2_vector(pt, lowmc.blocksize)
            c = lowmc.encrypt(p)
            print(c)
