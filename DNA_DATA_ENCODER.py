import datetime
import math
import pandas as pd
import random

begin_time = datetime.datetime.now()

def encoder(file_name, desired_oligo_length):
    df = pd.read_excel('puspabeethis_files\PrimerList_400_with_RC.xlsx', sheet_name='200_primer_pairs', header=0)
    primers = df['Primers'].tolist()
    file = open(file_name, "r")
    data = str(file.read())
    data_size = len(data)/8
    m_p = math.ceil((desired_oligo_length - 40)/10) - 1
    M = math.ceil(data_size/(2 * m_p))
    print('data size in Bytes = '+str(data_size))
    print('number of oligos = '+str(M))
    address = addressBook(M)
    codebook = subCode(M)
    oligos = []
    for i in range(M):
        payload = ''
        for j in range(m_p):
            message = 256*(128 * int(data[0]) + 64 * int(data[1]) + 32 * int(data[2]) + 16 * int(data[3]) + 8 * int(
                data[4]) + 4 * int(data[5]) + 2 * int(data[6]) + int(data[7])) + (128 * int(data[8]) + 64 * int(data[9]) + 32 * int(data[10]) + 16 * int(data[11]) + 8 * int(
                data[12]) + 4 * int(data[13]) + 2 * int(data[14]) + int(data[15]))
            data = data[16:]
            print(message, len(data))
            payload = payload + codebook[message]
        print(len(payload))
        n = random.randint(0, 199)
        primerF = primers[2*n + 0]
        primerRC = primers[2*n + 1]
        sequence = primerF + address[i] + payload + primerRC
        oligos.append(sequence)
    index = []
    for i in range(M):
        index.append('>oligo_' + str(i + 1))
    print(index)
    file_text = []
    for j in range(M):
        file_text.append(index[j])
        file_text.append(oligos[j])
    print(file_text)
    print(len(file_text))
    fasta_file = open("DNA_ENCODED_DATA.fasta", "w")
    for k in range(len(file_text)):
        fasta_file.write(file_text[k] + "\n")
    fasta_file.close()

def numberToBase(n,b):
    if n == 0:
        return '0'
    digits = []
    while n:
        digits.append(int(n % b))
        n //= b
    intTemp = digits[::-1]
    s = ""
    s = [s.join(str(j)) for j in intTemp]
    y = ""
    y = y.join(s)
    return y

def CountStringStreakBeginning(charVar):
    if len(charVar) == 1:
        return 1
    for i in range(1,len(charVar)):
        if charVar[i] != charVar[0]:
            return i
    return i + 1

def CountStringStreakEnding(charVar):
    if len(charVar) == 1:
        return 1
    for i in range(len(charVar) - 2,-1,-1):
        if charVar[i] != charVar[len(charVar) - 1]:
            return len(charVar) - (i + 1)
    return len(charVar) - i

def addressBook(M):
    if 1 <= M <= 8:
        n_a = 2
    elif 9 <= M <= 80:
        n_a = 4
    elif 81 <= M <= 504:
        n_a = 5
    elif 505 <= M <= 1008:
        n_a = 6
    elif 1009 <= M <= 6768:
        n_a = 7
    elif 6769 <= M <= 13278:
        n_a = 8
    elif 13279 <= M <= 93756:
        n_a = 9
    GC_LowerBound = math.ceil(.4 * n_a)
    GC_UpperBound = math.floor(.6 * n_a)
    l = math.floor(3 / 2)
    r = 3 - l
    ForbiddenString0 = '0'.rjust(3 + 1, '0')
    ForbiddenString1 = '1'.rjust(3 + 1, '1')
    ForbiddenString2 = '2'.rjust(3 + 1, '2')
    ForbiddenString3 = '3'.rjust(3 + 1, '3')
    TableData = []
    i = 0
    count = 0
    while i < 4 ** n_a :
       charVar = str(numberToBase(i, 4)).rjust(n_a, '0')
       if ((l != 0 and CountStringStreakBeginning(charVar) > l) or (r != 0 and CountStringStreakEnding(charVar) > r)) == False:
          if ForbiddenString0 not in charVar and ForbiddenString1 not in charVar and ForbiddenString2 not in charVar and ForbiddenString3 not in charVar:
             if (charVar.count('2') + charVar.count('3') >= GC_LowerBound) and (charVar.count('2') + charVar.count('3') <= GC_UpperBound):
                charVar = charVar.replace("0", "A").replace("1", "T").replace("2", "G").replace("3", "C")
                TableData.insert(count, charVar)
                count += 1
       i += 1
    if n_a == 8:
        df = pd.read_excel('puspabeethis_files\common_address_len8.xlsx', sheet_name='common_address', header=0)
        common_address = df['common_address'].tolist()
    elif n_a == 9:
        df = pd.read_excel('puspabeethis_files\common_address_len9.xlsx', sheet_name='common_address', header=0)
        common_address = df['common_address'].tolist()
    else:
        common_address = []
    for j in range(len(common_address)):
        TableData.remove(common_address[j])
    writer = pd.ExcelWriter('AddressBook.xlsx', engine='xlsxwriter')
    writer.save()
    df = pd.DataFrame({'Address': TableData})
    writer = pd.ExcelWriter('AddressBook.xlsx', engine='xlsxwriter')
    df.to_excel(writer, sheet_name='BlockLength'+str(n_a), index=False)
    writer.save()
    return TableData

def subCode(M):
    df = pd.read_excel('puspabeethis_files\CodeBook_large.xlsx', sheet_name='Min_HammingDistance_2', header=0)
    codebook_large = df['CodeBook_large'].tolist()
    if 9 <= M <= 80:
        df = pd.read_excel('puspabeethis_files\common_codeword_len4.xlsx', sheet_name='common_codeword', header=0)
        common_codeword = df['common_codeword'].tolist()
    elif 81 <= M <= 504:
        df = pd.read_excel('puspabeethis_files\common_codeword_len5.xlsx', sheet_name='common_codeword', header=0)
        common_codeword = df['common_codeword'].tolist()
    elif 505 <= M <= 1008:
        df = pd.read_excel('puspabeethis_files\common_codeword_len6.xlsx', sheet_name='common_codeword', header=0)
        common_codeword = df['common_codeword'].tolist()
    elif 1009 <= M <= 6768:
        df = pd.read_excel('puspabeethis_files\common_codeword_len7.xlsx', sheet_name='common_codeword', header=0)
        common_codeword = df['common_codeword'].tolist()
    else:
        common_codeword = []
    for j in range(len(common_codeword)):
        codebook_large.remove(common_codeword[j])
    return codebook_large

print(datetime.datetime.now() - begin_time)