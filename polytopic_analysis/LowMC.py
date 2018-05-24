from __future__ import print_function
from __future__ import division

import sys
sys.path.append("../")
from generate_matrices import instantiate_matrix, grain_ssg

from sage.all import *
from utils import *
from copy import copy

SBOX_SIZE = 3
SBox = mq.SBox(0x00, 0x01, 0x03, 0x06, 0x07, 0x04, 0x05, 0x02)
invSBox = mq.SBox(0x00, 0x01, 0x07, 0x02, 0x05, 0x06, 0x03, 0x04)

class LowMC(object):
    def generate_random(self, generator_settings):
        r""" Generate a random lowmc instance

        Generates a LowMC instance using the random generator options
        of the reference implementation.

        :param generator_settings: (``dict``) --
           dictionary that holds the parameters for the instance.
           Settings are:
           * *blocksize* (``int``) default: 32
           * *keysize* (``int``) default: 32
           * *num_sboxes* (``int``) number of s-boxes per round, default: 1
           * *rounds* (``int``) default: 24
           * *key* (``int``, ``string``), default: generate random key
        """

        self.blocksize = generator_settings.get('blocksize', 32)
        self.keysize = generator_settings.get('keysize', 32)
        self.num_sboxes = generator_settings.get('num_sboxes', 1)
        self.rounds = generator_settings.get('rounds', 24)

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

        self.set_key(generator_settings.get('key', [next(gen) for _ in range(self.keysize)]))

    def from_description_file(self, lowmc_instance_description):
        r""" Generate a lowmc instance from description file

        :param lowmc_instance_description: (``dict``) --
           dictionary that holds the parsed yaml file content:
           * *settings* sub dict holding instance settings
               * *blocksize* (``int``)
               * *keysize* (``int``)
               * *num_sboxes* (``int``) number of s-boxes per round
               * *rounds* (``int``)
               * *key* (``int``, ``string``)
           * *linear_layers* (list of matrixes)
           * *roundkey_matrices* (list of matrixes)
           * *round_constants* (list of vectors)
        """

        self.blocksize = lowmc_instance_description['settings']['blocksize']
        self.keysize = lowmc_instance_description['settings']['keysize']
        self.num_sboxes = lowmc_instance_description['settings']['num_sboxes']
        self.rounds = lowmc_instance_description['settings']['rounds']

        print("Create matrixes from yaml data ", end="")
        self.affine_matrixes = [matrix(GF(2), M) for M in lowmc_instance_description['linear_layers']]
        self.inv_affine_matrixes = [M**-1 for M in self.affine_matrixes]
        self.key_matrixes = [matrix(GF(2), M) for M in lowmc_instance_description['roundkey_matrices']]
        self.round_constants = [vector(GF(2), v) for v in lowmc_instance_description['round_constants']]
        print("[   Done   ]")

        self.set_key(lowmc_instance_description['settings']['key'])

    def __init__(self, lowmc_instance_description=None, generator_settings=None):
        if lowmc_instance_description:
            self.from_description_file(lowmc_instance_description)
        elif generator_settings:
            self.generate_random(generator_settings)

    def set_key(self, key):
        r""" Setup the key cycle for the given key

        This uses the key matrixes (setup at init stage) to
        generate both the full roundkeys and the reduced roundkeys.

        :param key: (``int``, ``string``)

        """
        debug_output("Setup roundkeys from key ", 1, end="")
        #print(key)
        key = to_gf2_vector(key, self.keysize)
        self.round_keys = []
        for r in range(self.rounds + 1):
            self.round_keys.append(self.key_matrixes[r]*key)

        # Calculate reduced round keys
        # Every reduced round key is influenced by the round key after
        # it and the linear layer of the next round.
        zv = VectorSpace(GF(2), self.blocksize).zero()
        # TODO: do we need this initial copy? (It does not look like we ever do a read access)
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
        r""" Substitution round (s-boxes)

        :param input: (``GF(2) vector``)
        :param inverse: (``bool``) default: False
        """
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
        r""" Encrypt the plaintext for given number of rounds

        :param input: (``GF(2) vector``)
        :param rounds: (``int``)
            If none, the input is encrypted over all rounds (This is the default)
            This assumes rounds <= number of rounds setup in the init step

        :returns: cyphertext (``GF(2) vector``)
        """
        if rounds is not None:
            rd = rounds
            if rounds > self.rounds:
                raise ValueError("Number of encryption rounds requested is too high.")
        else:
            rd =self.rounds

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
        r""" Encrypt the plaintext for given number of rounds using the reduced round keys

        :param input: (``GF(2) vector``)
        :param rounds: (``int``)
            If none, the input is encrypted over all rounds (This is the default)
            This assumes rounds <= number of rounds setup in the init step

        :returns: cyphertext (``GF(2) vector``)
        """
        if rounds is not None:
            rd = rounds
            if rounds > self.rounds:
                raise ValueError("Number of encryption rounds requested is too high.")
        else:
            rd =self.rounds

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
        r""" Decrypt the plaintext for given number of rounds

        :param input: (``GF(2) vector``)
        :param rounds: (``int``)
            If none, the input is encrypted over all rounds (This is the default)
            This assumes rounds <= number of rounds setup in the init step

        :returns: plaintext (``GF(2) vector``)
        """
        if rounds is not None:
            rd = rounds
            if rounds > self.rounds:
                raise ValueError("Number of decryption rounds requested is too high.")
        else:
            rd =self.rounds

        pt = copy(input)
        for rn in range(rd):
            r = self.rounds-rn
            pt += self.round_keys[r]
            pt += self.round_constants[r-1]
            pt =  self.inv_affine_matrixes[r-1]*pt
            pt =  self.substitution(pt, inverse=True);
        return pt + self.round_keys[0]

    def decrypt_round_reduced(self, input, round, round_key=None):
        r""" Decrypt single round using reduced round keys

        :param input: (``GF(2) vector``)
        :param round: (``int``)
            Determines the round constands (or more precisely the index of them)
            This assumes rounds <= number of rounds setup in the init step
        :param round_key: (``GF(2) vector``)
            If given, this round key will be used.
            If None (=default) the roundkey setup in the setupKey step is used

        :returns: decrypted (``GF(2) vector``)
        """
        if round > self.rounds-1:
            raise ValueError("Number of decryption rounds requested is too high.")

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
        r""" Decrypt the plaintext for given number of rounds using reduced round keys

        :param input: (``GF(2) vector``)
        :param rounds: (``int``)
            If none, the input is encrypted over all rounds (This is the default)
            This assumes rounds <= number of rounds setup in the init step

        :returns: plaintext (``GF(2) vector``)
        """
        if rounds is not None:
            rd = rounds
            if rounds > self.rounds:
                raise ValueError("Number of decryption rounds requested is too high.")
        else:
            rd =self.rounds

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
