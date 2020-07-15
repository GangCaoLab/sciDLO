import typing as t
import random
from itertools import islice

import editdistance.bycython
edit_dist = editdistance.bycython.eval
from Bio.Seq import reverse_complement as rc

class KeepDist(object):

    def __init__(self,
                 d: float,
                 dist_func: t.Callable = edit_dist
                 ):
        self.set = []
        self.d = d
        self.dist_func = dist_func

    def add(self, e):
        self.set.append(e)

    def filter(self, up_stream: t.Iterable) -> t.Iterable:
        """Let elements in stream keep distance with existing elements."""
        for e in up_stream:
            for e2 in self.set:
                if self.dist_func(e, e2) <= self.d:
                    break
            else:
                self.add(e)
                yield e

def barcode_gen(code_len, code_book='ATCG'):
    codes = set()
    while 1:
        code = []
        for _ in range(code_len):
            code.append(random.choice(code_book))
        code = "".join(code)
        if code not in codes:
            codes.add(code)
            yield code

def max_rep(seq):
    res = 1
    cnt = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i-1]:
            cnt += 1
        else:
            if cnt > res:
                res = cnt
            cnt = 1
    if cnt > res:
        res = cnt
    return res


def filter_tendem(gen):
    for code in gen:
        if max_rep(code) > 2:
            continue
        else:
            yield code

if __name__ == "__main__":
    codes = barcode_gen(8)
    codes = filter_tendem(codes)
    acc = KeepDist(3)
    codes = acc.filter(codes)
    codes = islice(codes, 96)
    #print("index\tbarcode\tlinker_F\tlinker_R")
    for i, code in enumerate(codes):
        a1 = "GTCGGA" + code + "G"
        a2 = rc(a1)[::-1]
        a1 = "TA" + a1
        a2 = a2 + "GATC"
        a2 = a2[::-1]
        #items = [str(i), code, a1, a2]
        #print("\t".join(items))
        if i < 10:
            idx = '0' + str(i)
        else:
            idx = str(i)
        print(f"MseI-linker-{idx}-F\t{a1}")
        print(f"MseI-linker-{idx}-R\t{a2}")
