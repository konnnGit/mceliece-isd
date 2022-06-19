#!/python
import sys
import time
import numpy as np
import os
import sympy
import gzip
import pickle
import random
import itertools
import math

'''General utility functions'''

#def συνάρτηση να γραφει απλο κείμενο
def myWriteFile(output, filename):
    with gzip.open(filename, 'wb') as f:
        f.write(pickle.dumps(output))

def mod(x, modulus):
    numer, denom = x.as_numer_denom()
    try:
        return numer*sympy.mod_inverse(denom,modulus) % modulus

    except:
        print('Error: Unable to apply modulus to matrix')

        exit()
def clear():
    os.system('clear')

def myReadFromFile(filename):
    with gzip.open(filename, 'rb') as f:
        matrix= sympy.Matrix(pickle.loads(f.read()))
    return matrix

def permutationMatrix(size):
    P=sympy.zeros(size,size)
    position=[]
    flag=False
    i=0
    while i< size:
        pos=np.random.randint(size)
        if pos not in position:
            position.append(pos)
            i+=1
    #print (position)
    for i in range(size):
        P[i,position[i]]=1
    return P

def doLowZeros(matrix,currentRow, currentCol):
    for row in range(currentRow+1,matrix.shape[0]):
        if  matrix[row,currentCol] :
            matrix=matrix.elementary_row_op(op='n->n+km',k=1,row1=row,row2=currentRow)
            matrix=matrix%2
    return matrix
def doUpperZeros(matrix,currentRow, currentCol):
    for row in range(currentRow-1,-1,-1):
        if  matrix[row,currentCol] :
            matrix=matrix.elementary_row_op(op='n->n+km',k=1,row1=row,row2=currentRow)
            matrix=matrix%2
    return matrix
def swapRows(matrix,currentRow, currentCol):
    for row in range(currentRow,matrix.shape[0]):
        if matrix[row,currentCol]:
            matrix=matrix.elementary_row_op(op='n<->m',row1=currentRow,row2=row)
    return matrix
def isLower0(m,colStart,colStop):
    falseList=[]
    #sympy.pprint(m)
    for r in range(m.shape[0]):
        for c in range(m.shape[0]):
            if c<r and m[r,c]!=0:
                falseList.append(999)
    if 999 in falseList:
        #matrix has 1 in lower triangl
        return False
    else:
        #matrix has 0 in lower triangl,Bingo....
        return True
def getZeroXY(matrix):
    #matrix=M[:,colStart:colStop]
    xList=[]
    yList=[]
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[0]):
            if i==j and matrix[i,j]==0:
                #print('i,j',i,j)
                xList.append(i)
                yList.append(j)
    return xList,yList
#fix the matrix size into the functions
def fixDiagonal(M,colStart,colStop):
    matrix=M[:,colStart:colStop]
    sampleDiag=sympy.ones(1,colStop-colStart)
    #while
    xL,yL=getZeroXY(matrix)
    #print(xL,yL)
    while xL:
        x=xL.pop()
        y=yL.pop()
        #for r in range(matrix.shape[0]):
        r=np.random.randint(matrix.shape[0])
        #print(matrix[r,y])
        if matrix[r,y]:
            M=M.elementary_row_op(op='n->n+km',k=1,row1=x,row2=r)
            M=M%2

    if sampleDiag!=matrix.diagonal():
        return False, M
    else:
        return True,M
'''Gauss elimination with row operations'''
def Gauss_Elim(matrix,colStart, colStop):
    flag1=0
    flag2=0
    while not flag1 and not flag2:
        #sympy.pprint(matrix)
        row=-1
        for col in range(colStart,colStop):
            row+=1
            if matrix[row,col]:
                matrix=doLowZeros(matrix,row,col)
            else:
                matrix=swapRows(matrix,row,col)
        flag2,matrix=fixDiagonal(matrix,colStart,colStop)
        flag1=isLower0(matrix[:,colStart:colStop],colStart,colStop)
    row=matrix.shape[0]
    for col in range(colStop-1,colStart, -1):
        row-=1
        #print('row,col', row,col)
        matrix=doUpperZeros(matrix,row,col)

    if matrix[:,colStart:colStop]==sympy.eye(colStop-colStart):
        return 1,matrix
    else:
        return 0,matrix



#---------- Creating Sets----------------------------------------------------------
def randVec(size,ones):
    vec=sympy.zeros(1,size)
    posList=[]
    i=0
    while i < ones:
        pos=random.randint(0,size-1)
        if pos not in posList:
            posList.append(pos)
            i+=1
    for j in posList:
        vec[0,j]=1
    return vec


def getPossibleVectors(size,ones):     
    veclist=[]
    over=(math.factorial(ones)*math.factorial(size-ones))
    listSize=int(math.factorial(size)/over)
    isSize=0
    #print(listSize)
    while isSize<listSize:
        shuffled=randVec(size,ones)
        if shuffled  in veclist:
            continue
        else:
            veclist.append(shuffled)
            isSize+=1
    return veclist


#--------------BCD version 1---------
def getB_bc(k,p):
    #Seperate randomly the k to k1,k2
    if k%2==0:
        k1=random.randint(int(k/4),int(k/2))#The proble is when k1!=k2 the the possible vectors are asymetrical different in the sets
    else:
        k1=random.randint(int(k-1/4),int(k-1/2))
        
    
    #k1=int(2*k/3)
    #k1=3
    k2=k-k1
    p1=1
    p2=p-1
    return getPossibleVectors(k1,p1), getPossibleVectors(k2,p2)
    
    
#--------------BCD version 2---------   
def getBCD(k,p):
	p1=int(2*p/3)
	p2=p-p1
	#if k is odd seperates the k
    if k%2!=0:
        floor=int((k-1)/2)
        ceil=floor+1
        #floorList=list(itertools.product([int(0), int(1)], repeat=floor))
        #ceilList=list(itertools.product([int(0), int(1)], repeat=ceil))
        return getPossibleVectors(floor,p1), getPossibleVectors(ceil,p2)
    #k is not odd
    else:
    	return getPossibleVectors(int(k/2),p1), getPossibleVectors(int(k/2),p2)
    
	
#------------- Stern ---------------
def getB(k,p):
    #if k is odd seperates the k
    if k%2!=0:
        floor=int((k-1)/2)
        ceil=floor+1
        #floorList=list(itertools.product([int(0), int(1)], repeat=floor))
        #ceilList=list(itertools.product([int(0), int(1)], repeat=ceil))
        return getPossibleVectors(floor,p), getPossibleVectors(ceil,p)
    #k is not odd
    else:

        list=getPossibleVectors(int(k/2),p)
        #left_right=list(itertools.product([int(0), int(1)], repeat=int(k/2)))
        return list ,list
        
        
        
#-------- Creating Sets End-----------------------------------------------------
