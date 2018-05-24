from __future__ import print_function
from __future__ import division

from copy import copy
from utils import *

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
