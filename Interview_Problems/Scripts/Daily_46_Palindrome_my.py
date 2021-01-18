# This is Daily Coding problem #46 
# ''' find largest palindrome in a string'''


def q(s):
    r = [[]]  # unique part of pol
    t = [list(s[0])]  # tester: temporary stack of chars (string)

    ind = []  # indexa of candidates for deletion

    ismid = [0]  # is middle a char
    lmax = 0
    rmax = ''



    def update_max(j):

        # update lmax and rmax with r[j] if it has bigger length

        nonlocal lmax, rmax

        if j > len(r)-1 or len(r[j]) == 0: return

        l = 2 * len(r[j]) - ismid[j]  # check for need to update winner
        if l >= lmax: # equal because rmax needs updating in that case
            lmax = l
            rmax = ''
            # transform list r[j] into polyndrome string
            for x in r[j]: rmax += x
            rmax = rmax[-1:-len(rmax) + ismid[j] - 1:-1] + rmax
        return None

    for i in range(1, len(s)):  # run down the string

        # UPDATE ALL EXISTING CANDIDATES, INCLUDING THE ONE JUST CREATED

        for j in range(0, len(t)):  # run down t

            #
            if len(t[j]) == 0:
                update_max(j)
                ind.append(j)
                continue

            if t[j][-1] != s[i]:
                update_max(j)  # check
                ind.append(j)  # earmark for del
            else:
                r[j].append(t[j].pop())  # pop from top of t[j] to bottom of r[j]
                if len(t[j]) == 0:
                    update_max(j)  # check
                    ind.append(j)  # earmark for del

        # DELETE EARMARKED CANDIDATES
        if len(ind) > 0:
            for index in sorted(ind, reverse=True):  # must inv-sort ind to not mess up loop index
                del t[index]
                del r[index]
                del ismid[index]
            ind = []  # clear list of candidates for deletion

        # CREATE NEW CANDIDATE IF MATCH FOUND

        # add element to t: a list, containing chars from s, up to i-1 or to i-2 (inclusive)
        if s[i] == s[i - 1]:  # no mid char
            t.append(list(s[:i - 1]))
            r.append(list(s[i]))
            ismid.append(0)

        if i > 1 and s[i] == s[i - 2]:  # there is mid char (avoid looking at first character)
            t.append(list(s[:i - 2]))
            r.append(list(s[i - 1:i + 1]))
            ismid.append(1)

    # from unprocessed candidates, process only the longest, which can only be 1-st or 2-nd
    if len(r) > 0:
        update_max(0)
        update_max(1)

    return lmax, rmax, r, t, ismid


#     0123456789
s1 = 'qabbaabbaw'
s2 = 'bananas'
s3 = 'qqqqqqqq'
s4 = 'qqqqabbaqqqqabba'

print(q(s1))
print(q(s2))
print(q(s4))
