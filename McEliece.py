#!/usr/bin/python3
import itertools
import math
import datetime
import time
from McElieceKeygen import *

    
def encrypt(mfile, pubkey,cipherFile):    
    G_pub, t_val = readFromFile(pubkey)
    k=G_pub.shape[0]
    n=G_pub.shape[1]
    G_pub=sympy.Matrix(G_pub)
    message = myReadFromFile(mfile)
    #Create a null set of codewords
    cwords=sympy.zeros(message.shape[0],G_pub.shape[1])
    cwordsWE=sympy.zeros(message.shape[0],G_pub.shape[1])
    
    #Create error vector    
    errorv=sympy.zeros(1,G_pub.shape[1])
    positions=[]
    i=0
    while i<t_val:
        pos = random.randrange(G_pub.shape[1])
        if pos not in positions :
                positions.append(pos)
                errorv[0,pos]=1
                i=i+1
    
    #Encrypt the total message with public key
    depth=message.shape[0]
    level=0
    while level<depth:
        cwords[level,:]=(message[level,:] * G_pub).applyfunc(lambda x: mod(x,2))
        level=level+1
    #sympy.pprint(cwords) 
    #print()
    
    #Add the errors   
    level=0    
    while level<depth:
        cwordsWE[level,:]=(cwords[level,:] + errorv).applyfunc(lambda x: mod(x,2))
        level=level+1
        
    print("The error vector is: ")
    sympy.pprint(errorv) 
    print("The encrypted message is: ")   
    sympy.pprint( cwordsWE) 
    
    #Write codeword c to file
    writeCipher(cwordsWE,cipherFile)  
    
    dist=get_dist(cwordsWE[0,:],cwordsWE[1,:])
    
    with open('encrypt_logs.txt','a') as f:
        f.write('\n')
        f.write('\n')
        f.write(str(datetime.datetime.now()))
        f.write('\n')
        f.write('n,k,t='+str(n) +','+ str(k) +','+ str(t_val) +'\n'+ 'distance between cwordsWE[0] and cwordsWE[1] is '+str(dist)+'\nerror is '+sympy.pretty(errorv))
        f.close 
    
   

    
def decrypt(cfile, privkey, outfile):
    S_matrix, H_matrix, P_matrix, t_val = readFromFile(privkey)
    k_val = H_matrix.shape[0]
    n_val = H_matrix.shape[1]
    S_inverse = (S_matrix**-1).applyfunc(lambda x: mod(x,2))
    P_inverse = (P_matrix**-1).applyfunc(lambda x: mod(x,2))
    cipher = myReadFromFile(cfile)
    depth=cipher.shape[0]
    #sympy.pprint(cipher)
    errors_tbit = []
    numbers = []
    for i in range(H_matrix.shape[1]):
        numbers.append(i)
    # Generate t-bit errors and syndromes
    for bits in itertools.combinations(numbers, t_val):
        et = sympy.zeros(1, H_matrix.shape[1])
        for bit in bits:
            et[bit] = 1
        st = (et * H_matrix.T).applyfunc(lambda x: mod(x,2))
        errors_tbit.append([et, st])
    message = []
    s_zeros = sympy.zeros(1,H_matrix[0])
    level=0
    while level<depth:
        print(sympy.pretty(cipher[level,:]), end=' -> ')
        mSG = (cipher[level,:] * P_inverse).applyfunc(lambda x: mod(x,2))
        s_mSG = (mSG * H_matrix.T).applyfunc(lambda x: mod(x,2))
        if not args.f:
            for errors in errors_tbit:
                if errors[1] == s_mSG:
                    recover = (mSG + errors[0]).applyfunc(lambda x: mod(x,2))
                    s_recover = (recover * H_matrix.T).applyfunc(lambda x: mod(x,2))
                    mSG = recover
        mS = mSG.extract([0], list(range(k_val, n_val)))
        m = (mS *S_inverse).applyfunc(lambda x: mod(x,2))
        if(args.v): sympy.pprint(m)
        message.append(m)
        level=level+1
    writePlain(message, outfile)
    
    '''
if args.vv:
    args.v = True
    '''
if args.g and args.m and args.t and args.o:
    logs='logs_mceliece.csv'
    #timestmp=datetime.datetime.now().strftime("%Y%m%d-%H%M%S")    
    pre='m'+str(args.m)+'t'+str(args.t)    
    logDone=pre+' code generated,'
    logError=pre+' NO code generated,'
    try:
        ts1=datetime.datetime.now()
        keygen(args.m, args.t, args.o)
        ts2=datetime.datetime.now()

        with open(logs , 'a') as f:            
            f.write('\n')
            f.write(logDone)
            f.write(str(ts2-ts1))
            f.close() 
    except Exception as ex:
        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)  
        print(message)       
        with open(logs , 'a') as f:
            f.write('\n')            
            f.write(message)
            f.close()
                 
                  
elif args.e  and args.pub and args.o:
    encrypt(args.e, args.pub, args.o)
elif args.d and args.priv and args.o:
    decrypt(args.d, args.priv, args.o)
    
elif args.x and args.par:
    print('x*H.T= ')
    checkCode(args.x,args.par)
else:
    print(parser.format_help())
