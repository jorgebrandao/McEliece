from math import factorial
from math import floor
from math import log2
from math import log
import matplotlib.pyplot as plt
from numpy import zeros


def xCy(x, y): #para x > y
    a = factorial(x)
    b = factorial(y)
    c = factorial(x-y)
    div = a//(b*c)
    return div

#exemple from: https://christianepeters.wordpress.com/publications/tools/#isdfq

def fcost(k,n,w,q):
    log2q = log2(q)
    potimo = floor(w/2)
    mincost = 10000000
    x = floor(k/2)

    for p in range(1, potimo):
        anum = xCy(x, p)
        bnum = xCy(k-x, p)

        cte = int(floor(log(anum)/log(q) + p*log(q-1)/log(q)) + 10)

        for l in range(1, cte - 1):
            ops = ((n-k)**2)*(n+k) + ((0.5*k-p+1) + (anum+bnum)*(q-1)**p)*l + q/(q-1)*(w-2*p+1)*2*p*(1+(q-2)/(q-1))*anum*bnum*(q-1)**(2*p)/(q**l)
        
            prob = anum*bnum*xCy(n-k-l, w-2*p)/xCy(n,w)
            cost = log2(ops*log2q/prob)

            if cost < mincost:
                mincost = cost
                bestp = p
                bestl = l
    return print("log2ops=",log2(ops)),print("bestp=",bestp), print("bestl=",bestl),print("mincost=",mincost)


### Meu from https://christianepeters.files.wordpress.com/2012/10/20100526-pqcrypto-isdfq.pdf

def wfact(k, n, w, q):
    x = floor(k/2)
    c = xCy(n, w)

    log2q = log2(q)

    wfi = 100000000

    cost1 = ((n-k)**2)*(n+k)#Cost: Updating the matrix G

    for p in range(1, floor(w/2)):#floor(w/2)floor(k/2)
        a = xCy(x, p)

        cte = floor(log(a)/log(q) + p*log(q-1)/log(q) + 10)

        for l in range(1, cte):
            b = xCy(n-k-l, w-2*p)

            prob = (a**2 * b)/c #The chance that Stern’s algorithm finds e after the first

            cost2 = ((x-p+1)+2*a*(q-1)**p)*l # Cost: Hashing step

            cost3 = q/(q-1)*(w-2*p+1)*2*p*(1+(q-2)/(q-1))*(a**2)*((q-1)**(2*p))/(q**l)#Cost:  Collision handling

            totalcost = cost1 + cost2 + cost3

            wf = log2(totalcost*log2q/prob)

            #wf1 = log2(cost1/prob)
            #wf2 = log2(cost2/prob)
            #wf3 = log2(cost3/prob)

            if wf < wfi:
                wfi = wf
                bestp = p
                bestl = l
    return print("wf=", wfi),print("bestp=",bestp), print("bestl=",bestl), print("totalcost=",log2(totalcost))


def wfi(k, n, w, q):
    x = floor(k/2)
    c = xCy(n, w)

    log2q = log2(q)

    wfi = 100000000

    cost1 = ((n-k)**2)*(n+k)#Cost: Updating the matrix G

    for p in range(1, floor(w/2)):
        a = xCy(x, p)

        cte = floor(log(a)/log(q) + p*log(q-1)/log(q) + 10)

        for l in range(1, cte):
            b = xCy(n-k-l, w-2*p)

            prob = (a**2 * b)/c #The chance that Stern’s algorithm finds e after the first

            cost2 = ((x-p+1)+2*a*(q-1)**p)*l # Cost: Hashing step

            cost3 = q/(q-1)*(w-2*p+1)*2*p*(1+(q-2)/(q-1))*(a**2)*((q-1)**(2*p))/(q**l)#Cost:  Collision handling

            totalcost = cost1 + cost2 + cost3

            wf = log2(totalcost*log2q/prob)

            #wf1 = log2(cost1/prob)
            #wf2 = log2(cost2/prob)
            #wf3 = log2(cost3/prob)

            if wf < wfi:
                wfi = wf
                bestp = p
                bestl = l
    return  wfi



def wfplot(k, n, w, q):
    wfs = zeros(k-64)
    a = list(range(k-64))

    for i in a:
        wfs[i] = wfi(k-i, n, w, q)

    plt.scatter(a, wfs)

    return plt.show()

def wfplot(k, n, w, q):
    wfs = zeros(k-64)
    a = list(range(k-64))

    for i in a:
        wfs[i] = log2(q**(k-i)) + wfi(k-i, n, w, q)

    plt.scatter(a, wfs)

    return plt.show(), print(wfs)