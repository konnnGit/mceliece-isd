#!/usr/bin/python3
import argparse
import sympy
import gzip
import pickle
import itertools
import math
import datetime
import time
import sys
from ISD_Utilities import *
#from McElieceUtil import *




def stern(c,H,t,p,l):
    rawH=myReadFromFile(H) 
    n=rawH.shape[1] 
    k= rawH.shape[1]  - rawH.shape[0] 
    
    if 2*p >t or 2*p>k or p==0:
        print('wrong p distribution or number, exiting..')
        exit()

    cword_all=myReadFromFile(c) 
    cword=cword_all[0,:]
    syndr=(cword*rawH.T).applyfunc(lambda x: mod(x,2))
    attempts=0
    attemptsQ=0
    time.sleep(1)
    #Algorithm inits
    alg=1       
    while alg:
        break1=1
        break2=1
        attempts+=1
        print("Stern attempt number", attempts )
        print("Creating P...")
        #P=permutationMatrix(n)
        permutation = random.sample(range(n), n)
        P = sympy.Matrix(n, n, lambda i, j: int((permutation[i]-j)==0))
        HP=(rawH*P).applyfunc(lambda x: mod(x,2))
        print("Attempting G.E...")
        #check,primeH=Gauss_Elim(HP,k,n)                      
        if HP[:,k:n].det()!=0:
            attemptsQ+=1 
            #print('Finding Q(#', attemptsQ,')...' )
            try:
                Q=HP[:,k:n].inv_mod(2)                
            except ValueError:
                print('Unable to apply G.E, restarting...')
                #time.sleep(1)
                continue
            #time.sleep(1)
            leftPrH= (Q*HP[:,0:k]).applyfunc(lambda x: mod(x,2))            
            kLeft,kRight=getB(k,p)
            lenKLeft=len(kLeft[0])                

            A=leftPrH[0:l,0:lenKLeft]                
            B=leftPrH[0:l,lenKLeft:k]
            C=leftPrH[l:n-k,0:lenKLeft]
            D=leftPrH[l:n-k,lenKLeft:k]   
                        
            primeSyndr=(syndr*Q.T).applyfunc(lambda x: mod(x,2))
            pSyndrl=primeSyndr[:,0:l]
            pSyndrR=primeSyndr[:,l:n-k]        
            #errorV=sympy.zeros(1,n)
            for eL in kLeft:
                eLAT=(eL*A.T).applyfunc(lambda x: mod(x,2))
                for eR in kRight:
                    eRBT=(eR*B.T).applyfunc(lambda x: mod(x,2))
                    if pSyndrl==(eLAT+eRBT).applyfunc(lambda x: mod(x,2))  :
                        eLCTeRDT=((eL*C.T).applyfunc(lambda x: mod(x,2))+(eR*D.T).applyfunc(lambda x: mod(x,2))).applyfunc(lambda x: mod(x,2))
                        prErrR=(pSyndrR+eLCTeRDT).applyfunc(lambda x: mod(x,2))
                        
                        lzero=sympy.zeros(1,l)
                        primeErrorV=eL.row_join(eR).row_join(lzero).row_join(prErrR)  
                          
                        #isSyndr=(errorV*rawH.T).applyfunc(lambda x: mod(x,2))
                        if  int(t)==np.count_nonzero(primeErrorV):
                            errorV=(primeErrorV*P.T).applyfunc(lambda x: mod(x,2))
                            print("Success, wt(e)=w=",t, ", error vector found ",sympy.pretty(errorV))
                            break1=0 
                        else: 
                                print('Wrong error vector found.')                                         
                    if break1==0:
                        break2=0
                        break
                if break2==0:
                    alg=0
                    break
        else:
            print('Unable to apply G.E, the random H(n-k) submatrix not invertible, restarting...')  
    return n,k,attempts         
                 
             



parser = argparse.ArgumentParser()
parser.add_argument("-c", type=int,	help="times to repeat the alg")
parser.add_argument("-p", type=int,	help="num of errors in k/2")
parser.add_argument("-l", type=int,	help="size of dimension l(stern)")
parser.add_argument("-t", type=str, help="num of errors")
parser.add_argument("-m", type=str, help="m in GF")
args = parser.parse_args()

if args.m and args.t  and args.c and args.l and args.p:
    path='stats/'
    logs='stern.csv' 
    logPath=os.path.join(path,logs)   
    mt=args.m+args.t
    cd=mt+'.codeword'
    h=mt+'.Hpub'       
    with open(logPath , 'a') as f:            
        f.write('\n')
        f.write(str(datetime.datetime.now()))
        f.close()
            
    try:    
        for instant in range(args.c):
            stamp1=time.time()
            #def stern(c,H,t,p,l): 
            n,k,attempt=stern(cd,h,int(args.t),args.p,args.l)  
            #nomPr=round(math.comb(n-int(args.t),k)/math.comb(n,k),4)  
            #print(nomPr)  
            stamp2=time.time()        
            with open(logPath , 'a') as f:            
                f.write('\n')
                f.write('Stern OK,'+args.m+','+args.t+','+str(n)+','+str(k)+','+str(args.l)+','+str(args.p)+','+str(int(stamp2-stamp1))+','+str(round(1/attempt,4))      )     
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
     
