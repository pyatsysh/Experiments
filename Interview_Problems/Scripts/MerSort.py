"""
Merge-sort recursive algorithm. From MIT algorithm book.


We will import deque. Without it, we would need to sort in descending order and reverse lists,
because popping from front of list is not optimal in python.

2Do: I used deques here. Maybe, it is possible to implement this with stock python, using lists?
"""
from collections import deque
import itertools  # to slice deck
import math # to use floor()


def merge(l1, l2):
    """ Merges two sorted deques l1 and l2 into one"""
    l = deque([])
    while len(l1) and len(l2):
        if l1[0] < l2[0]:
            l.append(l1.popleft())
        elif l2[0] < l1[0]:
            l.append(l2.popleft())
        else:
            l.extend((l1.popleft(), l2.popleft()))
    l.extend(l1)
    l.extend(l2)
    return l


def merge_sort(l):
    if len(l) == 1: return l
    ind = math.floor(len(l)/2)
    l1 = merge_sort(deque(itertools.islice(l, 0, ind)))
    l2 = merge_sort(deque(itertools.islice(l, ind, len(l))))
    return merge(l1, l2)


l = deque([5, 3, 2, 1, 4, 5, 5])
print(merge_sort(l))

