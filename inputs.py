#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#######################################################################
################# Convolutional McEliece Cryptosystem #################
#######################################################################


from mce import*

ordemS = 20

############## Private Key

G = Matrix([
[1, 0, 0, 1, 1, 0],
[0, 1, 0, 1, 0, 1],
[0, 0, 1, 0, 1, 1]
])

H = np.array([
    [1,1,0],
    [1,0,1],
    [0,1,1],
    [1,0,0],
    [0,1,0],
    [0,0,1]])

r = G.shape[0]

S = SMatrixZ(r, ordemS) # (rows, maiorgrau)

T = TMatrixZ(G)

P = T.inv()

############## Public Key
G_ = S*G*P
t = 1


############### Encrypt
def encrypt(s, G_=G_, t=t):
    G__ = splitMatrix(G_, 0, ordemS + 1)#splitMatrix(Matriz, menorgrau, maiorgrau)
    colT = G_.shape[1]
    peso = t
    u = Matrix([str2numpy(s)])
    ci = multp4MAT(u, G__)
    c = addErro(ci, peso, colT)
    return c

############### Decrypt
def decrypt(c, S=S, T=T, H=H, ordemS=ordemS):
    TT = splitMatrix(T, -1, 1)
    cT = multp4MAT(c, TT)
    cTe = sumsind(cT, H)
    original = multp4Sdecrypt(cTe, S, 1, ordemS, T) #1 is the low order of S
    return numpy2str(original)
