""" Here I just test out quirks and language idioms """

# comprehension quarks
# nested_list = [[0, 1], [2, 3], [4, 5], [6, 7], [8, 9]]
# result = [x for sublist in nested_list for x in sublist if x % 3 == 0]
# b = [a for z in [1, 2] for x in [z] for a in [x]]

import re

str = 'purple alic.e-b@google.com monkey dishwasher'
match = re.search(r'[\w.-]+@[\w.-]+', str)
if match:
    print(match.group())  # 'alice-b@google.com'


str = 'purple alice@google.com, blah monkey bob@abc.com blah dishwasher'
tuples = re.findall(r'([\w\.-]+)@([\w\.-]+)', str)
print(tuples)  ## [('alice', 'google.com'), ('bob', 'abc.com')]
for tuple in tuples:
    print(tuple[0], end = " ")  ## username
    print(tuple[1])  ## host

