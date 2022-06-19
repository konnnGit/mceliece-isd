#!/usr/bin/python3
import random
import sympy
import numpy
from McElieceUtil import *
#from InfoSetUtilities import *


def myReadFromFile(filename):
    with gzip.open(filename, 'rb') as f:
        matrix= sympy.Matrix(pickle.loads(f.read()))
    return matrix
x = sympy.Symbol('x')
alpha = sympy.Symbol('a')
# Generate key pairs
def keygen(m_val, t_val, files):
    global x

    # Random irreducible polynomial of degree m
    k_vec = get_irreducible(m_val)
    
    if args.v: print("k(x) = ", sympy.Poly(k_vec, x, domain='FF(2)'))
    # Random irreducible polynomial of degree t
    g_vec = get_irreducible(t_val)
    if args.v: print("g(x) = ", sympy.Poly(g_vec, x, domain='FF(2)'))
    # Produce k*n generator matrix G for the code
    H_matrix = get_H(m_val, t_val, k_vec, g_vec)
    if args.v: print('\nH ='); sympy.pprint(H_matrix); print(H_matrix.shape)
    G_matrix = get_G(H_matrix[:,:])
    if args.v: print('\nG ='); sympy.pprint(G_matrix); print(G_matrix.shape)
    k_val = G_matrix.shape[0]    
    # Select a random k*k binary non-singular matrix S
    S_matrix = get_nonsingular(k_val)
    if args.v: print('\nS ='); sympy.pprint(S_matrix)
    # Select a random n*n permutation matrix P
    n_val = G_matrix.shape[1]
    permutation = random.sample(range(n_val), n_val)
    P_matrix = sympy.Matrix(n_val, n_val,
    						lambda i, j: int((permutation[i]-j)==0))
    if args.v: print('\nP ='); sympy.pprint(P_matrix)
    # Compute k*n matrix G_pub = SGP
    G_pub = (S_matrix * G_matrix * P_matrix).applyfunc(lambda x: mod(x,2))
    H_pub=(H_matrix*P_matrix).applyfunc(lambda x: mod(x,2))
    '''
    #my code to turn Gpub to systematic
    isSystem, G_pub_Sys=Gauss_Elim(G_pub, G_pub.shape[1]-G_pub.shape[0],G_pub.shape[1])
    if isSystem:
        print('G_pub turned systematic.')
        G_pub=G_pub_Sys 
    else:
        print('Unable to systematic G_pub, exit..')  
        exit() 
    '''
    if args.v: print('\nG_pub ='); sympy.pprint(G_pub); print(G_pub.shape)
    # Public key is (G_pub, t)
    # Private key is (S, G, P)
    # length of the Goppa Code
    if args.v: print("\nn = ", n_val)
    if args.v: print("k =", k_val)
    writeKeys(G_matrix,G_pub, t_val, S_matrix, H_matrix,H_pub, P_matrix, files)
    
# From the parity check matrix, create the Generator matrix
def get_G(bin_H):
    for i in range(bin_H.shape[0]):
        bin_H.col_del(0)
    bin_G = bin_H.T
    ident = sympy.eye(bin_G.shape[0])
    for i in range(bin_G.shape[0]):
        bin_G = bin_G.col_insert(bin_G.shape[1], ident.col(i))
    return(bin_G)
def fixup_H(bin_H, pivot):
    num_removed = 0
    for j in range(bin_H.shape[0]):
        if bin_H.row(j-num_removed) == sympy.zeros(1, bin_H.shape[1]):
            bin_H.row_del(j-num_removed)
            num_removed += 1
    for i in range(bin_H.shape[0]):
        if not i == pivot[i]:
            col_a = bin_H.col(i)
            col_b = bin_H.col(pivot[i])
            bin_H.col_del(i)
            bin_H = bin_H.col_insert(i, col_b)
            bin_H.col_del(pivot[i])
            bin_H = bin_H.col_insert(pivot[i], col_a)
    return(bin_H)
# Create a parity check matrix give then Goppa code information
def get_H(m_val, t_val, k_vec, g_vec): 
    support = []
    g_a = []
    k_poly = sympy.Poly(k_vec, alpha, domain='FF(2)')
    # Store the support of the code into a list
    for i in range(pow(2, m_val)):
        if i == 0:
            # Calculate g(0)
            a_poly = sympy.reduced(sympy.Poly(sympy.Poly(g_vec, x, 
            domain='FF(2)').eval(0), alpha, domain='FF(2)').args[0], 
            [k_poly])[1].set_modulus(2)
        else:
            # Calculate g(a^(i-1))
            a_poly = sympy.reduced(sympy.Poly(g_vec, alpha**(i-1), 
            domain='FF(2)').args[0], [k_poly])[1].set_modulus(2)
        if not a_poly.is_zero:
            # Only store if it is not zero
            if i == 0:
                support.append(sympy.reduced(0, [k_poly])[1].set_modulus(2))    
            else:
                support.append(sympy.reduced(alpha**(i-1),
                			   [k_poly])[1].set_modulus(2))
            g_a.append(a_poly)
    if args.vv: print('\nSupport:\n', support, '\n\ng(a_n):\n', g_a)
    col = []
    # Store the inverses of g(a)
    for element in g_a:
        inverse = sympy.invert(element, k_poly)
        col.append(inverse)
    if args.vv: print('\nInverses:\n', col)
    # Form the Parity check matrix
    poly_H = []
    for i in range(t_val):
        poly_H_row = []
        for j in range(len(support)):
            top = sympy.Poly.pow(support[j], i)
            product =  sympy.reduced(sympy.Poly.mul(top, col[j]),
            						 [k_poly])[1].set_modulus(2)
            poly_H_row.append(product)
        poly_H.append(poly_H_row)
    bin_H = sympy.zeros(t_val * m_val, len(support))
    # Turn the parity check matrix into a binary matrix
    for i in range(t_val):
        for j in range(len(support)):
            current_poly = poly_H[i][j].all_coeffs()
            current_len = len(current_poly)
            for k in range(current_len):
                try:
                    bin_H[(i*(m_val))+k,j] = current_poly[current_len-k-1]
                except:
                    sympy.pprint(bin_H)
                    print('i =', i, ', j =', j, ', k =', k)
                    exit()
    bin_H, pivot = bin_H.rref(iszerofunc=lambda x: x % 2==0)
    bin_H = bin_H.applyfunc(lambda x: mod(x,2))
    bin_H = fixup_H(bin_H, pivot)
    return(bin_H)
def get_nonsingular(k_val):
    # Randomly generate a k*k Matrix
    # If it is singular generate another
    while True:
        matrix = sympy.Matrix(numpy.random.choice([sympy.Integer(0),
        sympy.Integer(1)], size=(k_val, k_val), p=[1./3, 2./3]))
        if matrix.det() % 2:
            return matrix
def get_irreducible(t_val):
    x = sympy.Symbol('x')
    # Randomly generate a polynomial of degree t
    # If it is reducible, generate another
    while True:
        # Coefficient of x^t must be 1
        polylist = [1]
        # Randomly select coefficients
        for i in range(t_val-1):
            polylist.append(random.randint(0, 1))
        # Coeffiecient of x^0 must be 1
        polylist.append(1)
        if (sum(polylist) % 2) == 1:
            # Produce a polynomial from the list
            p = sympy.Poly(polylist, x, domain='FF(2)')
            # Check that this is irreducible
            if p.is_irreducible:
                # Return the coefficients in a list
                return(p.all_coeffs())
# From https://stackoverflow.com/questions/31190182/sympy-solving-matrices
# -in-a-finite-field
def mod(x, modulus):
    numer, denom = x.as_numer_denom()
    try:
        return numer*sympy.mod_inverse(denom,modulus) % modulus
    except:
        print('Error: Unable to apply modulus to matrix')
        exit()
