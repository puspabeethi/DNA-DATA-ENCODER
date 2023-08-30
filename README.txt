Download the repository in your system as the zipped folder DNA-DATA-ENCODER-main.zip and unzip it. 
The folder DNA-DATA-ENCODER-main contains this README file, two Python files : DNA-DATA-ENCODER.py and main.py, and two folders : puspabeethis_files and Output_files.

The file DNA-DATA-ENCODER.py contains the program for encoding user data into synthetic DNA molecules.
This file is imported in the main.py file as a python module and the desired DNA-encoded data is obtained by running this main.py file.

The program facilitates encoding binary data of maximum size 28,501,824 bits into short-length synthetic DNA molecules of length upto 239 nt.
The designed oligos have length equal or close to the desired oligo length specified by user, and make use of P primer pairs, where P is also specified by the user.
An oligo is of the form FORWARD PRIMER + ADDRESS BLOCK + PAYLOAD + (reverse complement of) REVERSE PRIMER.
User data can be encoded into a maximum of 93756 oligos of maximum length 239 nt with primer pair count not exceeding 200.

The main.py file takes the following three arguments, that are fed into the encoder function of DNA-DATA-ENCODER module.
1) argument 1: the .txt file containing the user data in the form of a single binary string of maximum length 28,501,824 bits.
(User needs to create a binary text file in the parent folder DNA-DATA-ENCODER-main and store the desired data in the aforementioned format in the text file.
The data size is b bits, where b is the Column no. of the last bit in Line 1.)
2) argument 2: the desired length of each oligo.
(This length should be atleast 60).
3) argument 3: the desired number of primer pairs P.
(The value of P must not exceed 200).
User needs to enter these three arguments of DNA-DATA-ENCODER.encoder function in line 6 of main.py file.

The folder puspabeethis_files contains the necessary files required for the operation of the Python module.

The folder Output_files contains a blank fasta file named DNA_encoded_data. 
After execution of the program, the output files are stored in this folder, that comprise four excel documents: AddressBook, CodeBook, ForwardPrimers, ReversePrimers, and the finally obtained DNA-encoded data which is populated in the aforementioned fasta file DNA_encoded_data.fasta
This fasta file contains all the oligo sequences that collectively store the user data.

The four excel documents in the Output_files folder are necessary for decoding the .txt file, that was encoded using DNA-DATA-ENCODER module.
1) AddressBook contains the ordered address sequences, where the sequence in the 'n'th row of the excel sheet is used as address of 'oligo_(n-1)' in the fasta file.
2) A message 'm' is encoded into the codeword in the '(m+2)'th row of CodeBook.xlsx, where m ranges from 0 to 65535.
(The sequences from row number 65538 of CodeBook.xlsx are used as fillers to achieve the fixed payload length in the last oligo, when the remaining messages for the last oligo fall short of achieving the determined payload length.)
3) ForwardPrimers.xlsx comprises P forward primer sequences, where the sequence in the 'p'th row of the excel sheet corresponds to the forward primer of the '(p-1)'th primer pair in our designed Primer_Pair_List of 200 primer pairs. Similarly, ReversePrimers.xlsx comprises P reverse primer sequences, and the sequence in the 'p'th row of this excel sheet corresponds to the reverse primer of the '(p-1)'th primer pair of Primer_Pair_List, out of which P are required for encoding given user data.

Decoding Procedure:
1) Extract the address block of an oligo and order it using the AddressBook in Output_files folder. An address sequence in the 'n'th row of AddressBook.xlsx refers to the '(n-1)'th oligo. Repeat this to obtain the complete ordered list of oligo sequences. 
Let the total number of oligos be denoted by M. Then the obtained list contains the sequences of 1st oligo through 'M'th oligo.
2) Extract the payload from the 1st oligo, and segregate it into ordered blocks of length 10 each. Let the number of blocks be denoted by x. Then the payload is the concatenation of x 10-length blocks and is of the form Block1_1 + Block1_2 + ...... + Block1_x. 
Prepare the ordered list [Block1_1, Block1_2, ...... , Block1_x].
3) Repeat the above step for the remaining M-1 oligos. For every successive oligo, the prepared ordered list is appended to the end of previously obtained ordered list. 
For example, ordered list for 2nd oligo is appended to that for 1st oligo as [Block1_1, Block1_2, ...... , Block1_x, Block2_1, Block2_2, ...... , Block2_x]. 
The resulting final list is [Block1_1, Block1_2, ...... , Block1_x, Block2_1, Block2_2, ...... , Block2_x, ...... ,BlockM_1, BlockM_2, ...... , BlockM_x], and comprises (x*M) number of 10-length sequences, say y_1 through y_(x*M) (as per the order).
4) Refer CodeBook.xlsx to obtain the message corresponding to y_1. To do this, locate the row number 'r' of y_1 in CodeBook.xlsx; then the message corresponding to y_1 is (r-2). Express (r-2) in 16-bit binary format, as described next.
Obtain the binary representation of (r-2), where the leftmost bit is the MSB and rightmost bit is the LSB. If (r-2) < 2^(15), add suitable number of zeros to the left of MSB to obtain 16-bit binary representation of (r-2). 
For example, if y_1 is located in the 17th row, then it corresponds to message 15, which is 0000000000001111 in 16-bit binary format.
Repeat the same for y_2 through y_(x*M) and obtain the ordered list of messages m_1 through m_(x*M).
5) Concatenate the obtained messages in proper order as m_1 + m_2 + ................ + m_(x*M).
This is the decoded user data.



