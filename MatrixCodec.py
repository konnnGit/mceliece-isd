#!/usr/bin/python3
import argparse
import sympy
import gzip
import pickle
from sympy.combinatorics import Permutation
from McElieceUtil import *

def encode(text, parity):
    h=readFromFile(parity)
    H=sympy.Matrix(h)
    k_val=H.shape[1]-H.shape[0]
    output = []
    current = sympy.Matrix()
    for c in text:
        ascii_val = ord(c)
        for i in range(6, -1, -1):
            if current.shape[1] == k_val:
                output.append(current)
                current = sympy.Matrix()
            bit = pow(2, i)
            if ascii_val >= bit:
                ascii_val -= bit
                if current.shape[1] == 0:
                    current = sympy.Matrix([1])
                else:
                    current = current.col_insert(current.shape[1], 
                    							 sympy.Matrix([1]))
            else:
                if current.shape[1] == 0:
                    current = sympy.Matrix([0])
                else:
                    current = current.col_insert(current.shape[1], 
                    							 sympy.Matrix([0]))
    nonempty = current.shape[1]
    final = sympy.Matrix([1] * nonempty).T
    if final.shape[1] == 0:
        final = sympy.Matrix([0] * nonempty).T
    elif not nonempty == k_val:
        current = current.row_join(sympy.Matrix([0] * (k_val - nonempty)).T)
        final = final.row_join(sympy.Matrix([0] * (k_val - nonempty)).T)
        mix = Permutation.random(final.shape[1]).array_form
        randomize = []
        for i in mix:
            randomize.append(final[i])
        final = sympy.Matrix(randomize).T
    output.append(current)
    output.append(final)
    if args.v:
        print('The encoded text is:')
        for o in output:            
            sympy.pprint(o)
    if args.o:
        with gzip.open(args.o, 'wb') as f:
            f.write(pickle.dumps(output))
def decode(filename):
    with gzip.open(filename, 'rb') as f:
        input = pickle.loads(f.read())
    pad = input[-1]
    used_bits = (pad * sympy.Matrix([1] * pad.shape[1]))[0]
    input = input[0:-1]
    for i in range(pad.shape[1] - used_bits):
        input[-1].col_del(input[-1].shape[1] - 1)
    bit_count = 6
    ascii_val = 0
    output = ''
    for m in input:
        for i in range(m.shape[1]):
            ascii_val += m[i] * pow(2, bit_count)
            bit_count -= 1
            if bit_count < 0:
                bit_count = 6
                output = output + chr(ascii_val)
                ascii_val = 0
    if args.v: print(output)
parser = argparse.ArgumentParser()
parser.add_argument("-o", type=str, help="File to store output")
parser.add_argument("-e", type=str,
					help="String to be encoded into matrices")
parser.add_argument("-d", type=str,
					help="File with matrices to be decoded to ASCII")
parser.add_argument("-k", type=int, help="Length of each matrix")
parser.add_argument("-v", help="Enable verbose mode", action="store_true")

parser.add_argument("-par", type=str, help="Parity matrix")

args = parser.parse_args()
if args.e and args.par:
    encode(args.e, args.par)
elif args.d:
    decode(args.d)
else:
    print(parser.format_help())