import datetime
import DNA_DATA_ENCODER as dna
import DNA_Random_Access_Encoder as ra

begin_time = datetime.datetime.now()

dna.encoder("test.txt", 200)
ra.encoder("C://Users//User//Desktop//Primer and Codebook design//9.Python_module//test_folder", 200)

print(datetime.datetime.now() - begin_time)
