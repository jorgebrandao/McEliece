########################################################################
################### Poly McEliece Cryptosystem - RS ####################
########################################################################

#########               IMPORT LIST

from sympy import *
import numpy as np
from numpy.linalg import inv


########################################################################
#########               GERAR CHAVE
def keygen(size, g):
    a = symbols('a')
    key = np.zeros(size, dtype=object)

    grau = Poly(g).total_degree()

    keydegree = 2**grau

    for i in range(size):
        lista = list(range(keydegree))
        exp = np.random.choice(lista)
        if exp == 0:
            key[i] = 0
        else:
            key[i] = rem(a**(exp), g, domain=GF(2))
    return key


#########               GERAR MATRIZ GERADORA
def matgen(k,n,g):#Matriz geradora G
    a = symbols('a')
    mat = np.zeros(k*n, dtype=object)


    for i in range(n):
        mat[i] = 1


    j=n
    for c in range(k-1):
        for i in range(n):
            mat[i+j] = prem(mat[i+j-n]*a**(i+1), g, domain=GF(2))
            #mat[i+j] = mat[i+j-n]*a**(i+1)
        j=j+n

    return mat.reshape(k,n)

########################################################################
##              GERAR MATRIZ S(z)
def SMatrixZ(r, grau):# (ainda é preciso confirmar (por observação) se existem Zs na matriz. Não coloquei isto a verificar automaticamente porque para matrizes grandes a probabilidade de não ter nenhum Z é muito baixa e não compensa (em termos de tempo) fazer a verificação)
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
########################################################################
##              GERAR MATRIZ T(z)
def TDiagZ(n):#criar diagonal para matriz T(z)
    z = symbols('z')
    c = n
    diagonal = [1, 1/z, z]
    diagTz = zeros(c)
    for i in range(c):
        diagTz[i,i] = np.random.choice(diagonal)
    return diagTz
def TMatrixZDiag(n):# pega numa matriz com diagonal 1/z, 1 ou z e coloca para cima, na mesma coluna, apenas coisas do mesmo tipo do elemento (i,i) dessa coluna ou zero (de forma aleatória)
    z = symbols('z')
    c = n
    ii1_z = [0, 1/z]
    ii1 = [0, 1]
    iiz = [0, z]
    diagTz = TDiagZ(n)
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
def remove(mx):#como restrição da matriz T temos que em cada linha não pode existir dois elmentos iguais logo, esta função vai deixar o primeiro e meter os seguintes iguais a zero
    r = mx.shape[0]
    c = mx.shape[1]
    for i in range(r-1):
        for j in range(i, c-1):
            for a in range(j+1, r):
                if mx[i,j] == mx[i,a]:
                    mx[i,a] = 0
    return mx
def TMatrixZ(n): #matriz T(Z) final. Faz a permutação das colunas de TMatrixZDiag(G)
    TzDiag = TMatrixZDiag(n)
    TzDiagRem = remove(TzDiag)
    while TzDiagRem.det() != 1:
        TzDiag = TMatrixZDiag(n)
        TzDiagRem = remove(TzDiag)

    TMat = TzDiagRem

    while TMat == TzDiagRem:
        TMat = Matrix((np.random.permutation(TzDiagRem.T)).T)
    return TMat
########################################################################
#     Separar uma matriz nas várias componentes mod2 (por exemplo o G_)
#     Neste caso, vai colocar a matriz de menor grau (menor expoente) como G__[0], etc...
def splitMatrix(G_, begin, end):#criou a matriz G_ final mod 2 separada nos varios expoentes
    z = symbols('z')

    expandedG_ = expand(G_) #expand poly of G_ (it's faster than simplify(G_))

    row = G_.shape[0]
    col = G_.shape[1]
    l = row*col #=len(G_)

    nmat = len(list(range(begin, end + 1)))#numero de matrizes de zeros que terei que colocar na matriz de matrizes

    G__ = [zeros(row, col) for j in range(nmat)]# cria uma matriz de matrizes

    for n in range(nmat):
        for i in range (l):
            G__[n][i] = rem((expandedG_[i].coeff(z, end -n)), g, domain=GF(2))#.coeff(x, n) gives the coefficient of x**n
    return G__
########################################################################
# multiplica uma mensagem por G__ (public key)
def multp4MAT(u, MAT, g):
    r = MAT[0].shape[0]
    c = MAT[0].shape[1]

    resto = len(u) % r
    if resto != 0: #acrescenta zeros ao vetor u (menssagem) de forma a que o resto da divisão de u por r seja 0
        u = np.append(u, np.zeros((r-resto), dtype=object))

    b = len(u)
    nmat = len(MAT)#numero de matrizes após o split
    nzeros = np.zeros((nmat-1)*r, dtype=object)
    u = np.append(nzeros, u)
    u = np.append(u, nzeros)

    MAT = np.array(MAT).reshape(r*nmat, c)

    tfinal = int((nmat-1)*c+b/r*c) #tamanho da mensagem final encryptada
    
    uMAT = np.zeros(tfinal, dtype=object)#vetor de zeros

    for i in range(int(tfinal/c)):
        uMAT[i*c:i*c+c] = np.dot(u[i*r: i*r + nmat*r], MAT)
    for i in range(tfinal):
        uMAT[i] = rem(uMAT[i], g, domain=GF(2))
    return uMAT
########################################################################
### Cria vetor de erros
def errorprob(x):# retorna um 1 com p<=x e 0 com p>x
    err = np.random.rand(1)
    if err <= x:
        return 1
    elif err > x:
        return 0

def tvetor(ci, t, n):#adiciona um erro com determinado peso t e garante que não há mais do que p erros por qqer intervalomin
    v = int(len(ci))
    proberro = t/(3*n)#uso 3 pois é o numero de mat de T

    verro = np.zeros(v+2*n, dtype=object)#uso 2 pois é o numero de mat de T - 1

    for j in range(0, v, n):
        for i in range(n):
            soma = sum(verro[j : j+2*n + i])
            if soma < t:
                verro[j + 2*n + i] = errorprob(proberro)
    return verro[2*n:]

def errvetor(ci, t, n):
    a = symbols('a')
    vetorT = tvetor(ci, t, n)

    v = int(len(ci))

    Fn = keygen(v, g)#vetor em Fn com alphas

    for i in range(v):
        vetorT[i] = vetorT[i]*Fn[i]

    return vetorT

def sumerr(ci, t, n, g):
    vetorT = errvetor(ci, t, n)
    vetorT = vetorT + ci
    for i in range(int(len(ci))):
        vetorT[i] = rem(vetorT[i], g, domain=GF(2))
    return  vetorT

########################################################################
#multiplica pelo S final
def multp4Sdecrypt(cTe, S, n, grauS, menorgrau=1):
    SS0 = splitMatrix(S, 1, 1)[0]
    SS = splitMatrix(S, menorgrau+1, grauS)
    r = SS[0].shape[0]
    nmat = int(len(SS) + 1)
    IS0 = np.array(SS0.inv())
    a = int(len(cTe)/n)

    ctefinal = np.zeros((a-nmat-3)*r, dtype=object)#3 = nmat T
    original = np.zeros((a-nmat-3)*r, dtype=object)#3 = nmat T
    for i in range(a-nmat-3):#(o 3 corresponde ao nmat de T)
        ctefinal[i*r : i*r+r] = cTe[(i+2)*n : (i+2)*n+r]#2 vem de nmatT-1

    nzeros = np.zeros((nmat-1)*r, dtype=object)
    ctefinal = np.append(nzeros, ctefinal)
    ctefinal = np.append(ctefinal, nzeros)

    SS = np.array(SS).reshape(r*(nmat-1), r)


    for i in range(a-nmat-3):
        ctefinal[(nmat-1+i)*r:(nmat+i)*r] = np.dot((ctefinal[(nmat-1+i)*r:(nmat+i)*r] - np.dot(ctefinal[i*r: (nmat-1+i)*r], SS)), IS0)
        original[i*r: (i+1)*r] = ctefinal[(nmat-1+i)*r:(nmat+i)*r]

    return original

#######################################################################
#calc
def invG(G):
    k = G.shape[0]
    invG = np.zeros(k*k, dtype=object).reshape(k,k)

    for i in range(k):
        invG[i] = G[i][:k]
    return invG

#######################################################################
############################ McEliece Poly ############################
#######################################################################
a, z = symbols('a, z')


#########      Edit only this parameters    ################

#primitiv pol g
g = a**3+a+1 #GF(8)
#g = a**6+a+1 #GF(64)

grauS = 3# max order of S

t = 1
k = 4
n = 7

size = k*3 #keysize

#############################################################

G = matgen(k,n,g)

S = SMatrixZ(k, grauS)

T = TMatrixZ(n)
P = T.inv()

G_ = S*G*P

G__ = splitMatrix(G_, 0, grauS + 1)


u = keygen(size, g)

ci = multp4MAT(u, G__, g) #sem erros

cie = sumerr(ci, t, n, g) # com erros
