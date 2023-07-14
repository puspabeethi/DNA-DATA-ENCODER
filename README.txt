Download the repository in your system as a zipped file. 
This file DNA-DATA-ENCODER-main.zip contains this README file, two Python files : DNA-DATA-ENCODER.py and main.py, and two folders : puspabeethis_files and Output_files.

The file DNA-DATA-ENCODER.py contains the program for encoding user data into synthetic DNA molecules of length atmost 239 nt.
User data can be encoded into a maximum of 93756 oligos of the form FORWARD PRIMER + ADDRESS BLOCK + PAYLOAD + REVERSE PRIMER (reverse complement)
This file is imported in the main.py file as a python module and the DNA-encoded data is obtained by running this main.py file.

The main.py file takes the following three inputs which are fed into the encoder function of DNA-DATA-ENCODER module:
1) argument 1: the .txt file containing the user data in the form of a single binary string of maximum length 28,501,824 bits.
2) argument 2: the desired length of each oligo.
3) argument 3: the desired number of primer pairs.
