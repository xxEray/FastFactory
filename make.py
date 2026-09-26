import random
import sys

sys.set_int_max_str_digits(100000)

def rand_bigint(digits):
	return int(str(random.randint(1, 9)) + ''.join(str(random.randint(0, 9)) for _ in range(digits - 1)))

a = rand_bigint(100)
b = rand_bigint(100)
c = a * b

with open('in.txt', 'w') as fin:
	print(a, b, file=fin, sep='\n')
with open('ans.txt', 'w') as fout:
	print(c, file=fout)
