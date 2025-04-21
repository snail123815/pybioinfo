import re
import time
from copy import deepcopy
from functools import reduce

from Bio import Phylo


def makeList(ele1, ele2):
    if not isinstance(ele1, list):
        ele1 = [ele1]
    if not isinstance(ele2, list):
        ele2 = [ele2]
    return ele1 + ele2


def flattenList(lst):
    # print(flattenList([[[2, 3], [42, 'lskdj', 3]], [23, 8888, 88]]))
    if all(not isinstance(ele, list) for ele in lst):
        return lst
    return flattenList(reduce(makeList, lst, []))


def getTimeStr():
    return time.strftime("%z, %a, %d %b %Y, %H:%M:%S", time.localtime())


def timeDiffStr(a):
    d = abs(time.time() - a)
    h = int(d // 3600)
    return str(h).zfill(2) + time.strftime(":%M:%S", time.gmtime(d))


def safe_name(name: str) -> str:
    return re.sub(r"[ _:,();{}+*'\"[\]\/\t\n]+", "_", name)
