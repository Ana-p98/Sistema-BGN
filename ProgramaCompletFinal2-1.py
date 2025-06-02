import random
from sympy import isprime
import gmpy2

# Punt neutre
O = [None, None]

#Funció que retorna el m.c.d. d'un nombre a i d'un nombre b que és el mòdul i l'invers de a mòdul b amb l'algoritme d'Euclides.
def euclides(a, b):
    if a == 0:
        return [b, 0]
    u0, u1 = 1, 0
    v0, v1 = 0, 1
    a = a % b
    while a != 0:
        q = b // a
        r = b - a * q
        u, v = u0 - q * u1, v0 - q * v1
        b, a = a, r
        u0, u1 = u1, u
        v0, v1 = v1, v
    return [b, v0]

#Funció per calcular l'invers d'un nombre a mòdul p utilitzant l'algoritme d'Euclides.
def inv(a, p):
    if a < 0:
        a %= p
    mcd, v = euclides(a, p)
    if mcd == 1:
        return v%p
    else:
        return None



#Funció que retorna la suma de dos punts P1 i P2: P1+P2
def suma(Punt1, Punt2, primer):
    if Punt1 == O:
        return Punt2
    if Punt2 == O:
        return Punt1
    x1, y1 = Punt1
    x2, y2 = Punt2
    if x1 == x2:
        if (y1 + y2) % primer == 0:
            return O
        inv_denom = inv(2 * y1, primer)
        if inv_denom is None:
            return O
        l = (3 * x1**2) * inv_denom % primer
    else:
        inv_denom = inv((x2 - x1) % primer, primer)
        if inv_denom is None:
            return O
        l = ((y2 - y1) * inv_denom) % primer
    x3 = (l**2 - x1 - x2) % primer
    y3 = (l * (x1 - x3) - y1) % primer
    return [x3, y3]
 
    
#Funció que retorna [m]P, és a dir la suma de P m vegades.
def mSuma(P, p, m):
    R = O
    Q = P
    while m > 0:
        if m % 2 == 1:
            R = suma(R, Q, p)
        Q = suma(Q, Q, p)
        m //= 2
    return R


#Funció que retorna la funció h_{P1,P2}, avaluada a un punt, en funció de si lambda és infinit (False) o no ho és (True).
def h(Punt1, Punt2, Punt_evaluar, primer):
    if Punt1 == O or Punt2 == O:
        return 1
    x, y = Punt_evaluar
    xP, yP = Punt1
    xQ, yQ = Punt2
    if xP == xQ:
        if (yP + yQ) % primer == 0:
            return (x - xP) % primer
        inv_denom = inv(2 * yP, primer)
        if inv_denom is None:
            return 0
        l = ((3 * xP**2) * inv_denom) % primer
    else:
        inv_denom = inv((xQ - xP) % primer, primer)
        if inv_denom is None:
            return 0
        l = ((yQ - yP) * inv_denom) % primer
    denom = inv((x + xP + xQ - l**2) % primer, primer)
    if denom is None:
        return 0
    return ((y - yP - l * (x - xP)) * denom) % primer


#Funció per calcular un nombre en binari
def binari(num):
    bits = []
    while num:
        bits.insert(0, num % 2)
        num //= 2
    return bits


#Algoritme de Miller. 

#Funció f de Miller evaluada a un punt
def miller(Punt, Punt_evaluar, primer, n):
    T = Punt
    f = 1
    bits = binari(n)
    for bit in bits[1:]:
        f = (f * f * h(T, T, Punt_evaluar, primer)) % primer
        T = suma(T, T, primer)
        if bit == 1:
            f = (f * h(T, Punt, Punt_evaluar, primer)) % primer
            T = suma(T, Punt, primer)
    return f


#Funció que retorna el negatiu d'un punt P, -P. 
def punt_negatiu(Punt, primer):
    if Punt!=O:
        return [Punt[0], (-Punt[1]) % primer]
    else:
        return O
    

#Funció que retorna l'aparellament de Weil
def weil(Punt1, Punt2, Punt_evaluar, primer, n):
    Punt_evaluar_neg = punt_negatiu(Punt_evaluar, primer)
    fpQS = miller(Punt1, suma(Punt2, Punt_evaluar, primer), primer, n)
    fpS = miller(Punt1, Punt_evaluar, primer, n)
    fqPS = miller(Punt2, suma(Punt1, Punt_evaluar_neg, primer), primer, n)
    fqS = miller(Punt2, Punt_evaluar_neg, primer, n)
    if 0 in [fpQS, fpS, fqPS, fqS]:
        return 0
    inv_fpS = inv(fpS, primer)
    inv_fqPS = inv(fqPS, primer)
    if inv_fpS is None or inv_fqPS is None:
        return 0
    return (fpQS * inv_fpS % primer) * (fqS * inv_fqPS % primer) % primer


#Funció que retorna un nombre primer de manera aleatòria 
def primer_aleatori(bits):
    while True:
        n = random.getrandbits(bits)
        if isprime(n):
            return n

#Funció que genera el nombre n a partir de la multiplicació de dos nombres primers aleatoris, q i r, i el nombre primer p amb la propietat 4ln-1.
def generar_p(bits):
    q=primer_aleatori(bits)
    r=primer_aleatori(bits)
    n=q*r
    l=1
    while True:
        p = 4*l*n - 1
        if gmpy2.is_prime(p):
            return q,r, n, p
        l += 1       


bits=5
resultat = generar_p(bits)
# q = resultat[0]
# r = resultat[1]
# n = resultat[2]
# p = resultat[3]
q = 2
r = 7
n = q * r  # 35
p = 223

print('q=',q)
print('r=',r)
print('n=',n)
print('p=',p)


#Funció que retorna els punts que pertanyen a la corba y^2=x^3+1
def grup_punts(primer):
    puntX=[]
    puntY=[]
    punt=[]
    for i in range(primer):
        a=(i**3+1)%primer
        for j in range(primer):
            b=(j**2)%primer
            if a == b:
                puntX.append(i)
                puntY.append(j)  
    
    for k in range(len(puntX)):
        punt.append([puntX[k],puntY[k]])
    punt.append(O)
    return punt


#Funció que retorna tots els punts d'un ordre en específic
def punts_amb_ordre(primer, ordre_total):
    Punts = grup_punts(primer)
    punts_valids = []

    for P in Punts:
        if P == O:
            continue #Triam un altre punt
        P_ordre = mSuma(P, primer, ordre_total)
        if P_ordre == O:
            # Comprovam que sigui l'ordre més petit
            ordre_minim = True
            for d in range(1, ordre_total):
                if ordre_total % d == 0 and mSuma(P, primer, d) == O:
                    ordre_minim = False
                    break
            if ordre_minim:
                punts_valids.append(P)

    return punts_valids

#Funció per trobar tots els punts que no tenen un ordre en específic
def punts_sense_ordre(primer, ordre_excloure):
    """
    Retorna els punts de la corba E(F_p) que NO tenen ordre exactament ordre_excloure
    """
    Punts = grup_punts(primer)
    punts_no_valids = []

    for P in Punts:
        if P == O:
            continue
        if mSuma(P, primer, ordre_excloure) != O:
            # Clarament no tenen ordre n
            punts_no_valids.append(P)
        else:
            # Si [n]P = O, potser tenen ordre més petit
            for d in range(1, ordre_excloure):
                if ordre_excloure % d == 0 and mSuma(P, primer, d) == O:
                    punts_no_valids.append(P)
                    break

    return punts_no_valids

punts_n=punts_sense_ordre(p, n)
S=punts_n[0]

#Funció que retorna dos punts aleatoris amb ordre n
def puntPiQ_triat(p, n):
    punts = punts_amb_ordre(p, n)

    if punts:
        P_triat = random.choice(punts) 
        Q_triat = random.choice(punts)
    else:
        raise ValueError("No s'han trobat punts d'ordre n")
        
    return P_triat, Q_triat

#Missatge per encriptar
m = random.randint(1, r-1)
print('m=',m)


punts=puntPiQ_triat(p, n)
P=punts[0]
Q_tilde=punts[1]
Q=mSuma(Q_tilde, p, r)


#Funció que encripta un missate m
def encriptar(m, P, Q, p, n):
    t = random.randint(1, n)
    C=suma(mSuma(P, p, m),mSuma(Q, p, t),p)
    return C

C=encriptar(m, P, Q, p, n)
print('m encriptat=',C)

#Funció que retorna el logaritme
def logaritme(P, C, p, r):
    for m in range(r):
        if mSuma(P, p, m) == C:
            return m
    raise ValueError("No s'ha trobat el logaritme")


#Funció que desencripta el missatge encriptat
def desencriptar(C, P, q, p, r):
    P_tilde = mSuma(P, p, q)         # [q]P
    C_tilde = mSuma(C, p, q)         # [q]C
    m = logaritme(P_tilde, C_tilde, p, r)
    return m

print('m desencriptat=',desencriptar(C, P, q, p, r))

#Trobam tots els betas que siguin arrels cúbiques de la unitat.

def beta(p):
    betas = []
    for i in range(2, p):  # empezar en 2 para excluir el 1 directamente
        if pow(i, 3, p) == 1:
            betas.append(i)
    return betas

print(beta(p))
#Fixam un beta aleatori.
c=random.choice(beta(p))

#Aplicam el mapa distorsionat.
def mapa(P,p,c):
    if P==[None,None]:
        return [None,None]
    else:
        [x,y]=P
        return [(c*x)%p,y]

#Finalment, realitzam l'aparellament de Weil modificat.
def WeilMod(Punt1,Punt2,Punt_eval,primer,n):
    e=weil(Punt1,mapa(Punt2,p,c),Punt_eval,primer,n)
    return e


#C1 és el primer missatge encriptat
#C2 és el segon missatge encriptat
#Funció que retorna la suma de dos missatges encriptats
def add(C1, C2, Q, p, n):
    t = random.randint(1, n)
    C_tilde = suma(C1, C2, p)
    Q_tilde = mSuma(Q, p, t)
    C = suma(C_tilde, Q_tilde, p)
    return C


#Funció que retorna la multiplicació de l'esquema BGN
def mult(C1, C2, Q, S, p, n):
    u = random.randint(1, n)
    e1 = WeilMod(C1, C2, S, p, n)
    e2 = WeilMod(Q, Q, S, p, n)
    if e1 == 0 or e2 == 0:
        return 0
    e2u = pow(e2, u, p)
    return (e1 * e2u) % p







                                               

