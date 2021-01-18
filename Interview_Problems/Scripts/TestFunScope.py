""" Understanding assignments and scoping in python.

1. Assigning to var is binding name to object.
"""
# 1 -----------------
#
# Numbers are value types, which means the actual values are immutable.
# If you do x=3; x += 2, you aren't turning the number 3 into the number 5;
# you're just making x point to 5 instead of 3. The 3 is still out there unchanged,
# and any variables pointing to it will still see 3 as their value.
# ------------------

# (1) this is the equivalent of x = y = 2
x = 2  # 'x' is bound to 2
y = x  # 'y' is bound to 2, NOT to x, because x is not an object. The object s 2
y += 1  # 'y' is bound to 3, i.e.,
# bound different object int(3) to name 'y'. This has no effect on x
print(x, y)

# (2) this is equivalent to i = j = [1, 2, 3]
i = [1, 2, 3]
j = i  # both, i and j, are bound to the same list object
i[0] = 5
print(j)

# No contradiction between (1) and (2):
i = [1, 2, 3]
j = i
i = 5 # now bind i to int(5)
print(j)