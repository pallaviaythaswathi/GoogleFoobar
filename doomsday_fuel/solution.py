#The problem is an Absorbing Markov Chain
#The solution uses Matrix Inverse using Gaussian Elimination

m = [[0, 1, 0, 0, 0, 1], [4, 0, 0, 3, 2, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]]
print(solution(m))

from fractions import Fraction

def lcm(a, b):
    return long(a*b/gcd(a,b))

def simplify(a, b):
    g = gcd(a, b)
    return Fraction(long(a/g), long(b/g))


def gcd(a, b):
    def gcd1(a, b):
        if b == 0:
            return a
        return gcd1(b, a%b)
    return gcd1(abs(a), abs(b))

def transform(m):
    sums = list(map(sum, m))
    bool_indices = list(map(lambda x: x == 0, sums))
    indices = set([i for i, x in enumerate(bool_indices) if x])
    m2 = []
    for i in range(len(m)):
        m2.append(list(map(lambda x: Fraction(0, 1) if(sums[i] == 0) else simplify(x, sums[i]), m[i])))
    transform_matrix = []
    zero_matrix = []
    for i in range(len(m2)):
        if i not in indices:
            transform_matrix.append(m2[i])
        else:
            zeros_matrix.append(m2[i])
    transform_matrix.extend(zeros_matrix)
    tmat = []
    for i in range(len(transform_matrix)):
        tmat.append([])
        extend_mat = []
        for j in range(len(transform_matrix)):
            if j not in indices:
                tmat[i].append(transform_matrix[i][j])
            else:
                extend_mat.append(transform_matrix[i][j])
        tmat[i].extend(extend_mat)
    return [tmat, len(zeros_matrix)]

def copy_mat(m):
    copy_matrix = []
    for i in range(len(m)):
        copy_matrix.append([])
        for j in range(len(m[i])):
            copy_matrix[i].append(Fraction(m[i][j].numerator, m[i][j].denominator))
    return copy_matrix

def gauss_elmination(m, values):
    matrix = copy_mat(m)
    for i in range(len(matrix)):
        index = -1
        for j in range(i, len(matrix)):
            if matrix[j][i].numerator != 0:
                index = j
                break
        if index == -1:
            raise ValueError('Gauss elimination failed!')
        matrix[i], matrix[index] = matrix[index], matrix[j]
        values[i], values[index] = values[index], values[i]
        for j in range(i+1, len(matrix)):
            if matrix[j][i].numerator == 0:
                continue
            ratio = -matrix[j][i]/matrix[i][i]
            for k in range(i, len(matrix)):
                matrix[j][k] += ratio * matrix[i][k]
            values[j] += ratio * values[i]
    res = [0 for i in range(len(matrix))]
    for i in range(len(matrix)):
        index = len(matrix) -1 -i
        end = len(matrix) - 1
        while end > index:
            values[index] -= matrix[index][end] * res[end]
            end -= 1
        res[index] = values[index]/matrix[index][index]
    return res

def transpose(m):
    transpose_matrix = []
    for i in range(len(m)):
        for j in range(len(m)):
            if i == 0:
                transpose_matrix.append([])
            transpose_matrix[j].append(m[i][j])
    return transpose_matrix

def inverse(m):
    tmatrix = transpose(m)
    mat_inverse = []
    for i in range(len(tmatrix)):
        values = [Fraction(int(i==j), 1) for j in range(len(m))]
        mat_inv.append(gauss_elmination(tmatrix, values))
    return mat_inverse

def mat_mult(m1, m2):
    res = []
    for i in range(len(m1)):
        res.append([])
        for j in range(len(m2[0])):
            res[i].append(Fraction(0, 1))
            for k in range(len(m1[0])):
                res[i][j] += m1[i][k] * m2[k][j]
    return res

def splitQR(m, lengthR):
    lengthQ = len(m) - lengthR
    Q = []
    R = []
    for i in range(lengthQ):
        Q.append([int(i==j)-m[i][j] for j in range(lengthQ)])
        R.append(m[i][lengthQ:])
    return [Q, R]

def solution(m):
    res = transform(m)
    if res[1] == len(m):
        return [1, 1]
    Q, R = splitQR(*res)
    inv = inverse(Q)
    res = mat_mult(inv, R)
    row = res[0]
    l = 1
    for item in row:
        l = lcm(l, item.denominator)
    res = list(map(lambda x: long(x.numerator*l/x.denominator), row))
    res.append(l)
    return res
