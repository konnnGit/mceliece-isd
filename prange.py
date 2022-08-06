    #!/usr/bin/python3
import argparse
import sympy
import gzip
import pickle
import itertools
import math
import datetime
import time
import os
from ISD_Utilities import *
#from McElieceUtil import *

def prange(c, H, t):
    rawH=myReadFromFile(H)
    n=rawH.shape[1]
    k= rawH.shape[1]  - rawH.shape[0]  
    cword_all=myReadFromFile(c)
    #Attempt to work with the first part of the message, which is the first codeword.
    cword=cword_all[0,:]
    syndr=(cword*rawH.T).applyfunc(lambda x: mod(x,2))
    alg=1
    attempts=0
    time.sleep(1)
    #Algorithm inits
    attemptsQ=0
    while alg:
        attempts+=1
        print("Prange attempt number", attempts )
        print("Creating P...")
        #P=permutationMatrix(rawH.shape[1])
        permutation = random.sample(range(n), n)
        P = sympy.Matrix(n, n,
                        lambda i, j: int((permutation[i]-j)==0))
        HP=(rawH*P).applyfunc(lambda x: mod(x,2))
        print("Attempting G.E...")
        #check,primeH=Gauss_Elim(HP,k,n)
        if HP[:,k:n].det()!=0:
            attemptsQ+=1 
            #print('Finding Q(num.', attemptsQ,')...' )
            try:
                Q=HP[:,k:n].inv_mod(2)  
                             
            except ValueError:
                print('Unable to apply G.E, restarting...')
                #time.sleep(1)
                continue
            #time.sleep(1)
            primeSyndr=(syndr*Q.T).applyfunc(lambda x: mod(x,2) )
            zeroVec=sympy.zeros(1,k)
            primeErrorV=zeroVec.row_join(primeSyndr)
            
            #isSyndr=(errorV*rawH.T).applyfunc(lambda x: mod(x,2))
            if int(t)==np.count_nonzero(primeErrorV)  :
                errorV=(primeErrorV*P.T).applyfunc(lambda x: mod(x,2))
                print("Success, wt(e)=w=",t, ", error vector found ",sympy.pretty(errorV))
                
                alg=0
                break
            else: 
                print('Wrong error vector found.')
        else:
            print('Unable to apply G.E, the random H(n-k) submatrix not invertible, restarting...')
    return n,k,attempts         

                 
parser = argparse.ArgumentParser()
parser.add_argument("-c", type=int,	help="Times to repeat the alg")
#parser.add_argument("-l", type=int,	help="size of dimension l(stern)")
parser.add_argument("-t", type=str, help="num of errors")
parser.add_argument("-m", type=str, help="m in GF")
args = parser.parse_args()

if args.m and args.t  and args.c:
    path='stats/'
    logs='prange.csv'   
    logPath=os.path.join(path,logs) 
    mt=args.m+args.t
    with open(logPath , 'a') as f:            
        f.write('\n')
        f.write(str(datetime.datetime.now()))
        f.close()
    cd=mt+'.codeword'
    h=mt+'.Hpub'     
    try:            
        for instant in range(args.c):
            stamp1=time.time()
            n,k,attempt=prange(cd,h,int(args.t))  
            #nomPr=round(math.comb(n-int(args.t),k)/math.comb(n,k),4)  
            #print(nomPr)  
            stamp2=time.time()        
            with open(logPath , 'a') as f:            
                f.write('\n')
                f.write('Prange OK,'+args.m+','+args.t+','+str(n)+','+str(k)+','+str(int(stamp2-stamp1))+','+str(round(1/attempt,4)))           
                f.close()                
    except Exception as ex:  
        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)
        print(message)         
        with open(logPath , 'a') as f:  
            f.write('\n') 
            f.write(message) 
            f.close()  
else :
    print('missing argument')        
