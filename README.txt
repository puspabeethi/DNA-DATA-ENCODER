Download the repository in your system as a zipped file. 
This file DNA-DATA-ENCODER-main.zip contains this README file, two Python files : DNA-DATA-ENCODER.py and main.py, and two folders : puspabeethis_files and Output_files.

The file DNA-DATA-ENCODER.py contains the program for encoding user data into synthetic DNA molecules.
This file is imported in the main.py file as a python module and the DNA-encoded data is obtained by running this main.py file.

This program facilitates encoding of binary data of maximum size 28,501,824 bits into synthetic DNA molecules called oligos.
The designed oligos have length equal or close to the desired oligo length specified by user, and make use of P primer pairs, where P is also specified by the user.
An oligo is of the form FORWARD PRIMER + ADDRESS BLOCK + PAYLOAD + REVERSE PRIMER (reverse complement).
Note : User data can be encoded into a maximum of 93756 oligos of maximum length 239 nt with primer pair count not exceeding 200.

The main.py file takes the following three arguments, that are fed into the encoder function of DNA-DATA-ENCODER module.
1) argument 1: the .txt file containing the user data in the form of a single binary string of maximum length 28,501,824 bits.
User needs to create a binary text file in the parent folder DNA-DATA-ENCODER-main and store the desired data in the aforementioned format in the text file.
2) argument 2: the desired length of each oligo.
3) argument 3: the desired number of primer pairs P.
User needs to enter these three arguments of DNA-DATA-ENCODER.encoder function in line 6 of main.py file.

The folder puspabeethis_files contains the necessary files required for the operation of the Python module.
The folder Output_files contains a blank fasta file named DNA_encoded_data. After execution of the program, the output files are stored in this folder that comprise
four excel documents: AddressBook, CodeBook, ForwardPrimers and ReversePrimers, and the final output DNA-encoded data stored in the fasta file DNA_encoded_data.fasta