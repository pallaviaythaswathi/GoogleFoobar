def solution(n):
    n = int(n)
    operations = 0
  
    while n > 1:
        #if n is even, divide by 2 using bit manipulation
        if n % 2 == 0:
            n = n >> 1
        else:
            #create as many 0 in the Least Significant Bit as possible.
            n = (n - 1) if (n == 3 or n % 4 == 1) else (n + 1)

        operations += 1
    return operations