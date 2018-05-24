from __future__ import print_function
from __future__ import division

from LowMC import *
from DDiff import *

# from ddt import *
from sage.all import *
from utils import *

class LowMCPoly(LowMC):
    def generate_random(self, generator_settings):
        LowMC.generate_random(self, generator_settings)

        self.max_ddiff_size = generator_settings.get('max_ddiff_size', 3)

        # self.diff_propagation_forward = LowMCPoly._init_diff_propagation(possible_out_d)
        # self.diff_propagation_backward = LowMCPoly._init_diff_propagation(possible_in_d)

        self.possible_reduced_keys = self.init_possible_reduced_keys()

    def from_description_file(self, lowmc_instance_description):
        LowMC.from_description_file(self, lowmc_instance_description)

        self.max_ddiff_size = lowmc_instance_description['settings'].get('max_ddiff_size', 3)

        # self.diff_propagation_forward = LowMCPoly._init_diff_propagation(possible_out_d)
        # self.diff_propagation_backward = LowMCPoly._init_diff_propagation(possible_in_d)

        self.possible_reduced_keys = self.init_possible_reduced_keys()

    def __init__(self, lowmc_instance_description=None, generator_settings=None):
        if lowmc_instance_description:
            self.from_description_file(lowmc_instance_description)
        elif generator_settings:
            self.generate_random(generator_settings)

        self.rounds_with_prop1 = floor((self.blocksize-log(self.max_ddiff_size+1,2))/(self.num_sboxes*SBOX_SIZE))


    def init_possible_reduced_keys(self, sbox_idx=0):
        r""" Recursive initialization for the set of all possible (reduced) round keys

        For every sbox there are 8 possible round keys (3 bits)
        To make the attack easier this function generates a helper datastructure holding all
        possible round keays depending on the number of sboxes.

        :param sbox_idx: (``int``) --
           parameter used for recursive calls. Do not use if calling from outside.
        """
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
        return pk

    # @staticmethod
    # def _init_diff_propagation(table):
    #     # TODO: Check if this would require knowledge of how many sboxes are in the cypher.
    #     result = {}
    #     for i in range(len(table)):
    #         idx = to_gf2_vector(i, SBOX_SIZE).row()[0]
    #         result[idx] = [to_gf2_vector(o, SBOX_SIZE) for o in table[i]]
    #     return result

    def propSpaceAfterInitRound(self):
        r""" Calculate all possible diffs that go through the first substitution layer

        Only needed for when the selected d-diff size is too large to get through more
        than one round. (Workaround for toy ciphers)
        """
        num_sbox_bits = SBOX_SIZE*self.num_sboxes

        zero_bits = [0]*num_sbox_bits

        remaining_bits = self.blocksize - num_sbox_bits

        prop_space = []
        for i in range(1, remaining_bits**2-1):
            prop_space.append([x for x in list('{0:0{width}b}'.format(i, width=remaining_bits))]+zero_bits)

        return prop_space

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
            old_rank = 0
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
        if self.rounds_with_prop1 == 1:
            self.posibility_space.insert(0, self.propSpaceAfterInitRound())
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
        self.num_ddiffs_after_round = [0]*self.rounds

        ddiff_current_round = in_ddiff

        ddiffs_after_round = [ddiff_current_round] # make it work for trails of lenght 1

        #print_list(ddiffs_after_round)

        for r in range(self.rounds_with_prop1):
            ddiff_current_round = self.propagate_ddiff_inactive(ddiff_current_round, r, inverse=False)
            debug_output('ddiff for inactive round {} calculated'.format(r), 1)
            self.num_ddiffs_after_round[r] = 1

        ddiffs_current_round = [ddiff_current_round]
        ddiffs_after_round = ddiffs_current_round


        for r in range(self.rounds_with_prop1, round):

            #print('round {}, last round? {}'.format(r, last_round))
            ddiffs_after_round = DDiffHashMap()

            for ddiff in ddiffs_current_round:
                for dd in self.propagate_ddiff(ddiff, r, inverse=False):
                    ddiffs_after_round.append(dd)
            debug_output('ddiffs after round {}: {}'.format(r, len(ddiffs_after_round)), 1)
            self.num_ddiffs_after_round[r] = len(ddiffs_after_round)
            ddiffs_current_round = ddiffs_after_round
            #print("Num ddiffs after round {} is {}".format(r, len(ddiffs_current_round)))

        return ddiffs_after_round

    def propagate_ddiff_backward_from_to_round(self, in_ddiff, from_round, to_round):
        self.num_ddiffs_before_round = [0]*self.rounds
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
            self.num_ddiffs_before_round[r] = len(ddiffs_after_round)
            #print('ddiffs before round {}: {}'.format(r, len(ddiffs_after_round)))

        return ddiffs_after_round
