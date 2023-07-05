import datetime
import pandas as pd

begin_time = datetime.datetime.now()

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
    return y.rjust(8, '0')

df = pd.read_excel('test_data\FinalPixelList_30000.xlsx', sheet_name='bitDepth_8_bits', header=0)
data = df['PixelList'].tolist()
file_text = ''
for i in range(len(data)):
    file_text = file_text + str(numberToBase(data[i], 2))

file = open("test.txt", "w")
file.write(file_text)
file.close()

print(datetime.datetime.now() - begin_time)