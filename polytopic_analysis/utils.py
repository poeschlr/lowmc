from __future__ import print_function
from __future__ import division

from sage.all import *


class NotEnoughDegreesOfFreedom(Exception):
    def __init__(self, required_degrees, degrees):
        self.required_degrees = required_degrees
        self.degrees = degrees

    def __str__(self):
        return "Posibility space not enough degrees of freedom. Would need at least {:d} has {:d}".format(self.required_degrees, self.degrees)

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

class OutputHandler():
    def __init__(self, verbosity_level, logfile):
        self.verbosity_level = verbosity_level
        self.logfile = logfile

    def print_list(lst, required_verbosity = 0):
        if self.verbosity_level >= required_verbosity:
            ml = 0
            for l in lst:
                s = str(l)
                if len(s) > ml:
                    ml = len(s)
                print(l, end='\n')
            print('-'*ml, end='\n\n')

    def debug_output(msg, required_verbosity = 0, end='\n'):
        if self.verbosity_level >= required_verbosity:
            print(msg, end=end)

    def write_log_ln(logline):
        print(logline)
        self.logfile.write(logline+'\n')
