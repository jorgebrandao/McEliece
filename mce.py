#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#######################################################################
################# Convolutional McEliece Cryptosystem #################
#######################################################################

#########               IMPORT LIST
import numpy as np
from numpy.linalg import inv
from numpy.linalg import det
from sympy import symbols, Matrix, zeros, eye, expand


#########               FUNCTIONS
# str2numpy and numpy2str
def str2bytes(s):
    return s.encode('utf-8')

def bytes2str(b):
    return b.decode('utf-8')

def bytes2numpy(b):
    return np.array([int(e) for e in ''.join([format(v, '08b') for v in b])])

def numpy2bytes(n):
    return bytes(bytearray([int(e, 2) for e in split_str(''.join([str(v) for v in n]))]))

def str2numpy(s):
    return bytes2numpy(str2bytes(s))

def numpy2str(n):
    return bytes2str(numpy2bytes(n))

def split_str(s):
    return [s[i:i+8] for i in range(0, len(s), 8)]



##              MATRIX S(z) GENERATOR
def SMatrixZ(r, grau):# 
    z = symbols('z')
    SMatD = eye(r)*z
    for j in range(1, r):
        for i in range(j):
            for k in range (2 , grau+1):
                SMatD[i,j] = SMatD[i,j]+np.random.choice([z**k, 0, 0])
    SMat = SMatD
    while SMat == SMatD:
        SMat = Matrix((np.random.permutation(SMatD.T)).T)
    return SMat



##              MATRIX T(z) GENERATOR
def TDiagZ(G):
    z = symbols('z')
    c = G.shape[1]
    diagonal = [1, 1/z, z]
    diagTz = zeros(c)
    for i in range(c):
        diagTz[i,i] = np.random.choice(diagonal)
    return diagTz
def TMatrixZDiag(G):
    z = symbols('z')
    c = G.shape[1]
    ii1_z = [0, 1/z]
    ii1 = [0, 1]
    iiz = [0, z]
    diagTz = TDiagZ(G)
    for j in range(1, c):
        if diagTz[j,j] == 1/z:
            for i in range(j):
                diagTz[i,j] =  np.random.choice(ii1_z)
        if diagTz[j,j] == 1:
            for i in range(j):
                diagTz[i,j] =  np.random.choice(ii1)
        if diagTz[j,j] == z:
            for i in range(j):
                diagTz[i,j] =  np.random.choice(iiz)
    return diagTz
def remove(mx):
    r = mx.shape[0]
    c = mx.shape[1]
    for i in range(r-1):
        for j in range(i, c-1):
            for a in range(j+1, r):
                if mx[i,j] == mx[i,a]:
                    mx[i,a] = 0
    return mx
def TMatrixZ(G):
    TzDiag = TMatrixZDiag(G)
    TzDiagRem = remove(TzDiag)
    while TzDiagRem.det() != 1:
        TzDiag = TMatrixZDiag(G)
        TzDiagRem = remove(TzDiag)

    TMat = TzDiagRem

    while TMat == TzDiagRem:
        TMat = Matrix((np.random.permutation(TzDiagRem.T)).T)
    return TMat




#     Split a matrix in the diferent coefs (mod2)
def splitMatrix(G_, begin, end):
    z = symbols('z')

    expandedG_ = expand(G_) #expand poly of G_ (it's faster than simplify(G_))

    row = G_.shape[0]
    col = G_.shape[1]
    l = row*col #=len(G_)

    nmat = len(list(range(begin, end + 1)))

    G__ = [zeros(row, col) for j in range(nmat)]

    for n in range(nmat):
        for i in range (l):
            G__[n][i] = (expandedG_[i].coeff(z, end -n)) % 2#.coeff(x, n) gives the coefficient of x**n
    return G__


# multiply a message for G__ (public key)
def multp4MAT(u, MAT):
    r = MAT[0].shape[0]
    c = MAT[0].shape[1]

    resto = len(u) % r
    if resto != 0:
        u = np.append(u, np.zeros((r-resto), dtype=int))

    a = len(u)
    nmat = len(MAT)
    nzeros = np.zeros((nmat-1)*r, dtype=int)
    u = np.append(nzeros, u)
    u = np.append(u, nzeros)

    MAT = np.array(MAT).reshape(r*nmat, c)

    tfinal = int((nmat-1)*c+a/r*c) #len final message
    
    uMAT = np.zeros(tfinal, dtype=int)

    for i in range(int(tfinal/c)):
        uMAT[i*c:i*c+c] = np.dot(u[i*r: i*r + nmat*r], MAT)
    return uMAT % 2


def errorprob(x):# returns 1 if p<=x and 0 if p>x
    err = np.random.rand(1)
    if err <= x:
        return 1
    elif err > x:
        return 0

#Add a error vector with error probability peso/intervalomin
def addErro(ci, peso, colT):
    v = int(len(ci))

    proberro = peso/(3*colT)#3 because nmatofT is 3

    verro = np.zeros(v+2*colT, dtype=int)#2  =  nmatofT - 1

    for j in range(0, v, colT):
        for i in range(colT):
            soma = sum(verro[j : j+2*colT + i])
            if soma < peso:
                verro[j + 2*colT + i] = errorprob(proberro)
    cierro = verro[2*colT:] + ci   

    return cierro % 2


#Calculates syndrom and sums error padron
#notes: this function only works for 1 error. For t>1 we need to input the error padron matrix manually
def sumsind(cT, H):
    r = H.shape[0]
    c = H.shape[1]

    a = int(len(cT)/r)

    poserro = np.eye(r, dtype=int)

    sind = np.zeros(a*c, dtype=int)

    for i in range(a):
        sind[i*c:i*c+c] = np.dot(cT[i*r:i*r+r], H) % 2

        for j in range(r):
            if np.array_equal(sind[i*c:i*c+c], H[j]):
                cT[i*r:i*r+r] = cT[i*r:i*r+r] + poserro[j]
    return cT % 2

#multiply for final S
def multp4Sdecrypt(cTe, S, menorgrau, maiorgrau, T):
    SS0 = splitMatrix(S, 1, 1)[0]
    SS = splitMatrix(S, menorgrau+1, maiorgrau)
    c = T.shape[1]
    r = SS[0].shape[0]
    nmat = int(len(SS) + 1)
    IS0 = np.array(SS0.inv())
    a = int(len(cTe)/c)

    ctefinal = np.zeros((a-nmat-3)*r, dtype=int)#3 = nmat T
    original = np.zeros((a-nmat-3)*r, dtype=int)#3 = nmat T
    for i in range(a-nmat-3):#(3 = nmat de T)
        ctefinal[i*r : i*r+r] = cTe[(i+2)*c : (i+2)*c+r]#2 = nmatT-1

    nzeros = np.zeros((nmat-1)*r, dtype=int)
    ctefinal = np.append(nzeros, ctefinal)
    ctefinal = np.append(ctefinal, nzeros)

    SS = np.array(SS).reshape(r*(nmat-1), r)


    for i in range(a-nmat-3):
        ctefinal[(nmat-1+i)*r:(nmat+i)*r] = np.dot((ctefinal[(nmat-1+i)*r:(nmat+i)*r] - np.dot(ctefinal[i*r: (nmat-1+i)*r], SS)), IS0)
        original[i*r: (i+1)*r] = ctefinal[(nmat-1+i)*r:(nmat+i)*r]

    return original % 2

