Version: Python 3.5.4



Example 1: (McEliece convolutional codes: with Hamming codes)

$ python
>>> from mceliece import *
>>> s = '65537 is a prime, but 65536 is a power of 2!'
>>> # to encrypt s
>>> encrypted = encrypt(s)
>>> 
>>> # to decrypt encrypted
>>> decrypted = decrypt(encrypted)
>>> 
>>> exit()
$


Example 2: (Mceliece convolutional codes - RS codes: generates the encoder G, the private and the public key and encrypts)
Before open edit the parameters, g, grauS, t, k,n

$ python
>>> from rsencoding import *
>>> #to see all the matrix that are created
>>> G
>>> S
>>> T
>>> P
>>> G_
>>>
>>> # to encrypt u that is random by default, but you can change it
>>> u
>>> ci
>>> cie
>>>
>>> exit()
$


Example3: calculate the Work Factor for ISD attacks

$ python
>>> from workfactor import *
>>> wfact(k, n, w, q)
>>>
>>> exit()
$
