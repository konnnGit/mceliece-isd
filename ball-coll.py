#!/usr/bin/python3
import argparse
import sympy
import gzip
import pickle
import itertools
import math
import datetime
from ISD_Utilities import *
#from McElieceUtil import *

def ballcoll(c,H,t,p1,p2,q1,q2,l):
    rawH=myReadFromFile(H) 
    n=rawH.shape[1]
    k= rawH.shape[1]  - rawH.shape[0]
    print('n,k=',n,k)
    #p=p1+p2,q=q1+q2
    if p1+p2+q1+q2 >t or q1+q2>l or p1+p2>k  :
        print('wrong p,q distribution, exiting...')
        exit()
    cword_all=myReadFromFile(c)  
    cword=cword_all[0,:]
    syndr=(cword*rawH.T).applyfunc(lambda x: mod(x,2))
    attempts=0
    attemptsQ=0
    #time.sleep(1)
    #Algorithm inits, 1st loop
    alg=1       
    while alg:
        break1=1
        break2=1
        break3=1
        break4=1
        attempts+=1
        print("Ball Collision attempt= number", attempts )
        print("Creating P...")
        #P=permutationMatrix(n)
        permutation = random.sample(range(n), n)
        P = sympy.Matrix(n, n,lambda i, j: int((permutation[i]-j)==0))
        HP=(rawH*P).applyfunc(lambda x: mod(x,2))
        print("Attempting G.E...")

        if HP[:,k:n].det()!=0:
            attemptsQ+=1 
            #print('Finding Q(#', attemptsQ,')...' )
            try:
                Q=HP[:,k:n].inv_mod(2)
            except ValueError:
                print('Unable to apply G.E, restarting...')
                continue
            
            leftPrH= (Q*HP[:,0:k]).applyfunc(lambda x: mod(x,2)) 
            L1,L2=getBCD3(k,p1,p2) 
            lenk1=len(L1[0])
                        
            A=leftPrH[0:l,0:lenk1]               
            #εδώ ειναι το προβλημα             
            B=leftPrH[0:l,lenk1:k]
            C=leftPrH[l:n-k,0:lenk1]
            D=leftPrH[l:n-k,lenk1:k] 

            primeSyndr=(syndr*Q.T).applyfunc(lambda x: mod(x,2))
            pSyndrl=primeSyndr[:,0:l]
            pSyndrR=primeSyndr[:,l:n-k]
            
            eMlist=[]
            l1=[]
            l2=[]
            '''
            if q1==0 and q2==0:
                #l1,l2=getB(l,0)
                eM=sympy.zeros(1,l)
                eMlist.append(eM)
            else:
                l1,l2=getBCD3(l,q1,q2)             
                for el1 in  l1:
                    for el2 in l2: 
                        eM=el1.row_join(el2)
                        eMlist.append(eM)
            
            '''
            l1,l2=getBCD3(l,q1,q2)
            for el1 in  l1:
                for el2 in l2: 
                    eM=el1.row_join(el2)
                    eMlist.append(eM)
            
            #2nd loop
            for ek1 in L1:
                ek1AT=(ek1*A.T).applyfunc(lambda x: mod(x,2))
                for ek2 in L2:
                    #print('Into the main loop..')
                    #Insert the early abort
                    #if np.count_nonzero(ek1)+np.count_nonzero(ek2)>t-q:
                        #continue
                    ek2BT=(ek2*B.T).applyfunc(lambda x: mod(x,2))
                    ek1ATek2BT=(ek1AT+ek2BT).applyfunc(lambda x: mod(x,2))
                   # print('first mult,ok')
                    for eM in  eMlist:                  
                        sL=(ek1ATek2BT+eM).applyfunc(lambda x: mod(x,2))
                        #print('2nd add,ok')
                        if pSyndrl==sL:
                            #print('equality,ok')                            
                            ek1CT=(ek1*C.T).applyfunc(lambda x: mod(x,2)) 
                            #print('ek1CT')                                                
                            ek2DT=(ek2*D.T).applyfunc(lambda x: mod(x,2))
                            #print('ek2DT') 
                            ek1CTek2DT=(ek1CT+ek2DT).applyfunc(lambda x: mod(x,2)) 
                            #sympy.pprint(ek1CTek2DT) 
                            #print(n-k-l) 
                            #sympy.pprint(pSyndrR)                         
                            prErrR=(ek1CTek2DT+pSyndrR).applyfunc(lambda x: mod(x,2)) 
                            #print('prErrR')                          
                            primeErrorV=ek1.row_join(ek2).row_join(eM).row_join(prErrR)  
                             
                            #isSyndr=(errorV*rawH.T).applyfunc(lambda x: mod(x,2))  and isSyndr==syndr
                            if int(t)==np.count_nonzero(primeErrorV):
                                errorV=(primeErrorV*P.T).applyfunc(lambda x: mod(x,2))
                                print("Success, wt(e)=w=",t, ", error vector found ",sympy.pretty(errorV))
                                break4=0  
                            else: 
                                print('Wrong error vector found.')                                                                   
                        if break4==0:
                            break3=0
                            break
                    if break3==0:
                        break2=0
                        break
                if break2==0:
                    alg=0
                    break
               
            
        else:
            print('Unable to apply G.E, the random H(n-k) submatrix not invertible, restarting...')   
              
    return n,k,lenk1,len(l1[0]),attempts         

                 
parser = argparse.ArgumentParser()
parser.add_argument("-c", type=int,	help="Times to repeat the alg")
parser.add_argument("-l", type=int,	help="size of dimension l(stern)")
parser.add_argument("-p1", type=int,	help="num of errors in k1. At least equals to one.")
parser.add_argument("-p2", type=int,	help="num of errors in k2. At least equals to one.")
parser.add_argument("-q1", type=int,	help="num of errors in l1. <zero if need to declare as zero")
parser.add_argument("-q2", type=int,	help="num of errors in l2. <zero if need to declare as zero")
parser.add_argument("-t", type=str, help="num of errors")
parser.add_argument("-m", type=str, help="m in GF")
args = parser.parse_args()

if args.m and args.t and args.l and args.c and args.p1 and args.p2 and args.q1 and args.q2:
    if args.q1<0:
        q1=0
    else:
        q1=args.q1
    if args.q2<0:
        q2=0
    else:
        q2=args.q2
        
    
    
    path='stats/'
    logs='ball-coll.csv' 
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
            n,k,k1,l1,attempt=ballcoll(cd,h,int(args.t),args.p1,args.p2,q1,q2,args.l) 
            stamp2=time.time()        
            with open(logPath , 'a') as f:            
                f.write('\n')
                f.write('Ball coll OK,'+args.m+','+args.t+','+str(n)+','+str(k)+','+str(k1)+','+str(args.p1)+','+str(args.p2)+','+str(q1)+','+str(q2)+','+str(l1)+','+str(int(stamp2-stamp1))+','+str(round(1/attempt,4)))          
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
