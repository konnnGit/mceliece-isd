The source code of this attempt is in https://github.com/konnnGit/mceliece-isd/tree/final.

The steps one may follow to deploy the overall application are:

1. Construct a McEliece cryptosystem based on the Galois Field parameters m
and t. For this paradigm m = 7 and t = 9.

		./McEliece.py -g -m 7 -t 9 -o 79 -v

	Note: The flag -o is used to define the name of the necessary binary files
	which contain the keys, etc. The m, t parameters are part of the file names.
	See manual pages A.6 for more information.
	Note: There is a restriction while creating keys. It is recommended not to use
	big values for m and t because it stalls the create process.

2. Encode a plain text message into a binary message using the public version of
the Parity Check Matrix H (Hpub = HP).

		./MatrixCodec.py -e ’Hello World’ -par 79.Hpub -v -o 79.binMsg


3. Encrypt the binary message into a codeword using the McEliece public key.

	McEliece.py -e 79.binMsg -pub 79.pub -o 79.codeword

	Note: One may use the included bash script

		./create_keys_cword.sh 7 9

	in order to execute the steps 1,2 and 3 simultaneously. In that case, the
	message to be sent is included in the script. Editing it one can be change it.

4. Attack to the cipher, having three (3) options:

	a) Prange’s algorithm
		./prange.py -c 5 -m 7 -t 9
		
	Note: In order to run the algorithms in background one may use the nohup command, like

		./nohup prange.py -c 5 -m 7 -t 9 &

b) Stern’s algorithm

		./stern.py -c 35 -m 7 -t 9 -p 1 -l 2

	Note: The flag -p concerns the error bits in each k/2 coordinates.
	Note: The flag -c indicates how many times will be repeated the algorithm.

c) Ball-collision-Decoding algorithm

		./ball-coll.py -c 35 -m 7 -t 9 -p1 2 -p2 1 -q1 -1 -q2 -1 -l 2

	Note: In order to pass zero errors at any of -q1, -q2 coordinates, this
	must declared with a negative value.


5. Decrypt the cipher.

	Decrypt the cipher. In case of decrypting, the flag -d must be selected. We
	coose the correspoding codeword, the private key and define the output palin
	text name.

		./McEliece.py -d 79.codeword -pub 79.priv -o 79.plain

6. Decode the decrypted binary cipher text into the initial plain text.

		./MatrixCodec.py -d 79.plain



Restrictions:
