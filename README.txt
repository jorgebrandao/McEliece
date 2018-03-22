Version: Python 3.5.4



Example:

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
