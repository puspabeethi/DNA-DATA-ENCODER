import datetime
import DNA_DATA_ENCODER as dna
import DNA_Random_Access_Encoder as ra

begin_time = datetime.datetime.now()

dna.encoder("test.txt", 200)

print(datetime.datetime.now() - begin_time)
