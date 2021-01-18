""" Insertion sort (from alg book) 

Using for loop is ugly. Here a while loop is more appropriate

The main difference is that if in a while loop I manually increment/decrement counter,
after the loop finishes, the return counter value is changed in preparation for the next
iteration

In a for-loop, it will be the last used counter

TODO: implement using a while loop
"""

import random 
l = [1, 2, 3, 4, 5, 5, 4, 0, 0]
print(l)
random.shuffle(l)
print(l)

# if you need to break out of a for-loop, consider using a while loop instead!
for i in range(1, len(l)):
    a = l[i]
    for j in range(i-1, -1, -1):  # from i-1 down to 0 in reverse order
        if l[j] > a:  # if l[j] is bigger than l[i], move l[j] up
            l[j+1] = l[j]
        else:
            j = j+1  # THIS IS UGLY! a while loop over j would avoid this.
                     # Without this, in line 21 would have l[j+1], which would not sort [2, 1]
            break
    l[j] = a

print(l) 


