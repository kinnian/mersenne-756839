import os

def V(B):
    L = len(B)
    x = sum(2**(8*i)*B[i] for i in range(L))
    return x

def to_string(x):
    B = "{0:b}".format(x)
    return B

def random_mod(m):
    v = m
    while v >= m:
        B = os.urandom(3)
        v = sum(2**(8*i)*B[i] for i in range(3))

    return v

def generate_h_sparse_string(B, h):
    for i in range(h):
        B[i] = 1
    for i in range(h-n):
        B[h+i] = 0
    i = h-1
    while i >= 0:
        j = random_mod(n-i)
        a = B[i]
        B[i] = B[i+j]
        B[i+j] = a
        i = i-1

    return B

def key_pair(h, K, P):
    A_f = []
    A_g = []
    A_f = generate_h_sparse_string(A_f, h)
    A_g = generate_h_sparse_string(A_g, h)
    A_r = os.urandom(K)

    f = V(A_f)
    g = V(A_g)
    R = V(A_r)
    R = R % P
    T = f*R + g
    T = T % P
    PK = to_string(R) + to_string(T)
    SK = to_string(f)

    return PK, SK

def kem_enc(h, PK, K, P, S, rho):
    SS = os.urandom(32)
    A_a = []
    A_b1 = []
    A_b2 = []
    A_a = random_mod(A_a, h)
    A_b1 = random_mod(A_b1, h)
    A_b2 = random_mod(A_b2, h)
    a = V(A_a)
    b1 = V(A_b1)
    b2 = V(A_b2)
    R = V(PK[0:K-1])
    T = V(PK[K:2*K-1])
    C1 = (a*R + b1) % P
    C2 = (a*T + b2) % P
    M = []
    for i in range(255):
        if S[i] == 0:
            for j in range (i*rho/8, (i+1)*rho/8 -1):
                M[j] = 0
        else:
            for j in range (i*rho/8, (i+1)*rho/8 -1):
                M[j] = 255

    CT = to_string(C1) + to_string(operator.xor(C2[0:32*rho - 1], M))
    return CT


