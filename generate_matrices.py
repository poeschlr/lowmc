import argparse
import yaml


def save_datafile(file, linlayers, round_constants, roundkey_matrices):
    with open(file, 'w') as matfile:
        s = 'LowMC matrices and constants\n'\
            '============================\n'\
            'Block size: ' + str(blocksize) + '\n'\
            'Key size: ' + str(keysize) + '\n'\
            'Rounds: ' + str(rounds) + '\n\n'
        matfile.write(s)
        s = 'Linear layer matrices\n'\
            '---------------------'
        matfile.write(s)
        for r in range(rounds):
            s = '\nLinear layer ' + str(r + 1) + ':\n'
            for row in linlayers[r]:
                s += str(row) + '\n'
            matfile.write(s)

        s = '\nRound constants\n'\
              '---------------------'
        matfile.write(s)
        for r in range(rounds):
            s = '\nRound constant ' + str(r + 1) + ':\n'
            s += str(round_constants[r]) + '\n'
            matfile.write(s)

        s = '\nRound key matrices\n'\
              '---------------------'
        matfile.write(s)
        for r in range(rounds + 1):
            s = '\nRound key matrix ' + str(r) + ':\n'
            for row in roundkey_matrices[r]:
                s += str(row) + '\n'
            matfile.write(s)

def save_as_yaml(file, linlayers, round_constants, roundkey_matrices):
    data = {
        'settings':{
            'blocksize': blocksize,
            'keysize': keysize,
            'rounds': rounds
        },
        'linear_layers': linlayers,
        'round_constants': round_constants,
        'roundkey_matrices': roundkey_matrices
    }
    with open(file, 'w') as outfile:
        yaml.dump(data, outfile)

def save_as_cpp(file, linlayers, round_constants, roundkey_matrices):
    cpp_file = file.replace('.h', '.cpp')
    with open(cpp_file, 'w') as matfile:
        s = '#include "' + file + '"\n'\
            '/*\n'\
            ' * LowMC matrices and constants\n'\
            ' ******************************/\n\n'#\
            # 'const unsigned numofboxes = ' + str(num_sboxes) + ';    // Number of Sboxes\n'\
            # 'const unsigned blocksize = ' + str(blocksize) + ';   // Block size in bits\n'\
            # 'const unsigned keysize = ' + str(keysize) + '; // Key size in bits\n'\
            # 'const unsigned rounds = ' + str(rounds) + '; // Number of rounds\n\n'\
            # 'const unsigned identitysize = blocksize - 3*numofboxes;\n\n'
        matfile.write(s)
        s = '/*\n'\
            ' * Linear layer matrices\n'\
            ' ********************/\n'
        matfile.write(s)
        s = 'std::string lgen_inst_lin_layer[{num_rounds:d}][{blocksize:d}] = {{\n'.format(
                num_rounds=rounds, blocksize=blocksize
            )
        for r in range(rounds):
            s += '  {{ // Linear Layer - Round {:d}\n'.format(r)
            for row in linlayers[r]:
                s += '    "{:s}",\n'.format(''.join(map(str, row)))
            s = s[:-2] + '\n'
            s += '  }}{:s}\n'.format(',' if r<rounds-1 else '')
        s += '};\n\n'
        matfile.write(s)

        s = '/*Round constants\n'\
            ' ********************/\n\n'
        matfile.write(s)
        s = 'std::string lgen_inst_round_constant[{num_rounds:d}] = {{\n'.format(num_rounds=rounds)
        for r in range(rounds):
            s += '  "{:s}"{:s} // Round {:d}\n'.format(
                ''.join(map(str,round_constants[r])),
                ',' if r<rounds-1 else ' ', r
            )
        s += '};\n\n'
        matfile.write(s)

        s = '/*Round key matrices\n'\
            ' ********************/\n\n'
        matfile.write(s)

        s = 'std::string lgen_inst_key_matrix[{num_rounds:d}][{blocksize:d}] = {{\n'.format(
                num_rounds=rounds+1, blocksize=blocksize
            )
        for r in range(rounds + 1):
            s += '  {{ // Round key matrix - Round {:d}\n'.format(r)
            for row in roundkey_matrices[r]:
                s += '    "{:s}",\n'.format(''.join(map(str, row)))
            s = s[:-2] + '\n'
            s += '  }}{:s}\n'.format(',' if r<rounds else '')
        s += '};\n\n'
        matfile.write(s)

    with open(file, 'w') as matfile:
        s = '#ifndef _LOWMC_EXTERNAL_DEFINITION_H_\n'\
            '#define _LOWMC_EXTERNAL_DEFINITION_H_\n\n'\
            '#include <string>\n'\
            '#include <vector>\n\n'\
            '/*\n'\
            ' * LowMC matrices and constants\n'\
            ' ******************************/\n\n'\
            'const unsigned numofboxes = ' + str(num_sboxes) + ';    // Number of Sboxes\n'\
            'const unsigned blocksize = ' + str(blocksize) + ';   // Block size in bits\n'\
            'const unsigned keysize = ' + str(keysize) + '; // Key size in bits\n'\
            'const unsigned rounds = ' + str(rounds) + '; // Number of rounds\n\n'\
            'const unsigned identitysize = blocksize - 3*numofboxes;\n\n'
        matfile.write(s)
        s = '/*\n'\
            ' * Linear layer matrices\n'\
            ' ********************/\n'
        matfile.write(s)
        s = 'extern std::string lgen_inst_lin_layer[{num_rounds:d}][{blocksize:d}];'.format(
                num_rounds=rounds, blocksize=blocksize
            )
        matfile.write(s)

        s = '/*Round constants\n'\
            ' ********************/\n\n'
        matfile.write(s)
        s = 'extern std::string lgen_inst_round_constant[{num_rounds:d}];'.format(
                num_rounds=rounds, blocksize=blocksize
            )
        matfile.write(s)

        s = '/*Round key matrices\n'\
            ' ********************/\n\n'
        matfile.write(s)

        s = 'extern std::string lgen_inst_key_matrix[{num_rounds:d}][{blocksize:d}];'.format(
                num_rounds=rounds+1, blocksize=blocksize
            )
        matfile.write(s)

        s = '\n\n#endif //_LOWMC_EXTERNAL_DEFINITION_H_'
        matfile.write(s)

def main():
    ''' Use the global parameters `blocksize`, `keysize` and `rounds`
        to create the set of matrices and constants for the corresponding
        LowMC instance. Save those in a file named
        `matrices_and_constants.dat`.
    '''
    gen = grain_ssg()
    linlayers = []
    for _ in range(rounds):
        linlayers.append(instantiate_matrix(blocksize, blocksize, gen))

    round_constants = []
    for _ in range(rounds):
        constant = [next(gen) for _ in range(blocksize)]
        round_constants.append(constant)

    roundkey_matrices = []
    for _ in range(rounds + 1):
        mat = instantiate_matrix(blocksize, keysize, gen)
        roundkey_matrices.append(mat)


    for file in output_file:
        if file.lower().endswith('.yaml'):
            save_as_yaml(file, linlayers, round_constants, roundkey_matrices)
        elif file.lower().endswith('.h'):
            save_as_cpp(file, linlayers, round_constants, roundkey_matrices)
        else:
            save_datafile(file, linlayers, round_constants, roundkey_matrices)

def instantiate_matrix(n, m, gen):
    ''' Instantiate a matrix of maximal rank using bits from the
        generatator `gen`.
    '''
    while True:
        mat = []
        for _ in range(n):
            row = []
            for _ in range(m):
                row.append(next(gen))
            mat.append(row)
        if rank(mat) >= min(n, m):
            return mat

def rank(matrix):
    ''' Determine the rank of a binary matrix. '''
    # Copy matrix
    mat = [[x for x in row] for row in matrix]

    n = len(matrix)
    m = len(matrix[0])
    for c in range(m):
        if c > n - 1:
            return n
        r = c
        while mat[r][c] != 1:
            r += 1
            if r >= n:
                return c
        mat[c], mat[r] = mat[r], mat[c]
        for r in range(c + 1, n):
            if mat[r][c] == 1:
                for j in range(m):
                    mat[r][j] ^= mat[c][j]
    return m


def grain_ssg():
    ''' A generator for using the Grain LSFR in a self-shrinking generator. '''
    state = [1 for _ in range(80)]
    index = 0
    # Discard first 160 bits
    for _ in range(160):
        state[index] ^= state[(index + 13) % 80] ^ state[(index + 23) % 80]\
                        ^ state[(index + 38) % 80] ^ state[(index + 51) % 80]\
                        ^ state[(index + 62) % 80]
        index += 1
        index %= 80
    choice = False
    while True:
        state[index] ^= state[(index + 13) % 80] ^ state[(index + 23) % 80]\
                        ^ state[(index + 38) % 80] ^ state[(index + 51) % 80]\
                        ^ state[(index + 62) % 80]
        choice = state[index]
        index += 1
        index %= 80
        state[index] ^= state[(index + 13) % 80] ^ state[(index + 23) % 80]\
                        ^ state[(index + 38) % 80] ^ state[(index + 51) % 80]\
                        ^ state[(index + 62) % 80]
        if choice == 1:
            yield state[index]
        index += 1
        index %= 80

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate the constands for a LowMC instance.')
    parser.add_argument('-b', '--blocksize', type=int, nargs='?', default=256)
    parser.add_argument('-k', '--keysize', type=int, nargs='?', default=80)
    parser.add_argument('-r', '--rounds', type=int, nargs='?', default=12)
    parser.add_argument('-s', '--num_sboxes', type=int, nargs='?', default=1)
    parser.add_argument('-o','--output', type=str, nargs='+', default='matrices_and_constants.dat')
    args = parser.parse_args()

    blocksize = args.blocksize
    keysize = args.keysize
    rounds = args.rounds
    output_file = args.output
    num_sboxes = args.num_sboxes

    main()
