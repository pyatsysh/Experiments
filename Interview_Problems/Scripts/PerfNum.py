""" Daily coding problem 70. Perfect numbers.

Perfect number is such that its digits add up to 10. Given n, find n-th perfect number.

Lets obtain a list of perfect nums up to n. It is easier to test.
We will first brute-force it to get a test case. Then we will smart it.
"""

# Note that for perf number p, p-10 is divisible by 9.


def sum_digits(n):  # should be re-written using %-operator. Also, implement recursively
    s = 0
    p = n
    while p > 0:
        p = p // 10
        s += n-10*p
        n = p
    return s


def perf_num_brute(n): # brute-force approach
    p = 10
    counter = 0
    while counter < n:
        p += 9
        if sum_digits(p) == 10:
            counter += 1
            print(p)


def perf_num_smart(n): # smart way
    if n <= 9:
        return 19+(n-1)*9
    a, i, s = 1, 9, 1
    while i < n:
        s = sum_digits(a)
        i = i+10-s+1
        a += 1
    # roll back last change
    a -= 1
    i -= 10-s+1
    k = n-i-1
    return a*100+(10-s+9*k)


def main():
    perf_num_brute(30)

    print(perf_num_smart(30))

main()
