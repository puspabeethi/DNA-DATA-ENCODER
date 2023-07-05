import datetime
import random
import pandas as pd
import math
import numpy as np
from PIL import Image

begin_time = datetime.datetime.now()

img = Image.open(r"C:\Users\User\Desktop\September 2021\Primer and Codebook design\7.Codewords\Final_image_grayscale.PNG")
#img = img.convert("L")
#img.save("Final_image_grayscale.PNG")

pic = img.load()
width, height = img.size
print(str(width) + 'x' + str(height))
pixelList = []
for i in range(width):
    for j in range(height):
        x = pic[i, j]
        print(x)
        pixelList.append(x)

writer = pd.ExcelWriter('FinalPixelList_30000.xlsx', engine='xlsxwriter')
writer.save()
df = pd.DataFrame({'PixelList': pixelList})
writer = pd.ExcelWriter('FinalPixelList_30000.xlsx', engine='xlsxwriter')
df.to_excel(writer, sheet_name='bitDepth_8_bits', index=False)
writer.save()

#coordinate = x, y = 150, 5
#print(img.getpixel(coordinate));

print(datetime.datetime.now() - begin_time)
