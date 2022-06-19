#!/bin/bash
m=$1
t=$2  
   
./McEliece.py -g -m $m -t $t -o $m$t -v
txt='sasa'
./MatrixCodec.py -e $txt -par $m$t.Hpub -v -o $m$t.binMsg    
./McEliece.py -e $m$t.binMsg   -pub $m$t.pub -o $m$t.codeword 
