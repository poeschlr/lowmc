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

from utils import *
from sage.all import *
import argparse
import yaml
import traceback
import time

from LowMCPoly import *

REFERENCE_PT = 0

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

    rb = lowmc.blocksize - lowmc.num_sboxes*3
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
                ck = vector(GF(2), [0]*rb + list((c + a)[-lowmc.num_sboxes*3:]))
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

    write_log_ln('{bs:d}, {ks:d}, {nr:d}, {dds:d}, {dd_f:s}, {dd_b:s}, '\
            '{t_find:.3f}, {t_enc:.3f}, {t_forward:.3f}, '\
            '{t_backward:.3f}, {t_collision:.3f}, {t_roundkeys:.3f}, '\
            '{num_collisions:d}, {num_trails_forward:d}, {num_trails_backward:d}, '\
            '{key_ok:s}'.format(
        bs=lowmc.blocksize, ks=lowmc.keysize, nr=lowmc.rounds, dds=3,
        dd_f=str(lowmc.num_ddiffs_after_round[:round_mid]).replace(',',';'),
        dd_b=str(lowmc.num_ddiffs_before_round[round_mid:]).replace(',',';'),
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
        lowmc_instance['settings']['max_ddiff_size'] = args.ddiff_size


        lowmc = LowMCPoly(lowmc_instance_description=lowmc_instance)
        attack(lowmc, logfile, args.only_trail)

    elif args.auto:
        with open(args.auto, 'r') as config_stream:
            try:
                auto_def = yaml.load(config_stream)
            except yaml.YAMLError as exc:
                print(exc)
                exit()
        with open('log -- {}.csv'.format(time.strftime("%Y_%m_%d - %H_%M_%S", time.localtime())), 'w') as logfile:
            write_log_ln('D-Diff size, {:d}'.format(args.ddiff_size), logfile)
            write_log_ln('{bs:s}, {ks:s}, {nr:s}, {dds:s}, {dd_f:s}, {dd_b:s}, '\
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
                    lowmc = LowMCPoly(generator_settings=definition['setup'])
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
