import sys
import os
import numpy
import format_output as fo
def Hamming74Encode(input_data):
	data_count = 4
	output = 0
	parity_count = 3
	data_bits = numpy.zeros(data_count)
	parity_bits = numpy.zeros(parity_count)
	for index in range(data_count):
		data_bits[index] = input_data & 1
		input_data = input_data >> 1
	
	even = 0
	
	parity_bits[2] = data_bits[0] + data_bits[1] + data_bits[2]                + even
	parity_bits[1] = data_bits[0] + data_bits[1]                + data_bits[3] + even
	parity_bits[0] = data_bits[0]                + data_bits[2] + data_bits[3] + even
	
	interleaved_data = numpy.zeros(data_count + parity_count)
	sequential_data = numpy.zeros(data_count + parity_count)
	
	interleaved_data[0] = parity_bits[0] % 2
	interleaved_data[1] = parity_bits[1] % 2
	interleaved_data[2] = data_bits[0]
	interleaved_data[3] = parity_bits[2] % 2
	interleaved_data[4] = data_bits[1]
	interleaved_data[5] = data_bits[2]
	interleaved_data[6] = data_bits[3]
	
	interleaved_codeword = 0
	for index in range(data_count + parity_count):
		interleaved_codeword = interleaved_codeword + int(interleaved_data[index] * (2**index))
	
	sequential_data[0] = data_bits[0]
	sequential_data[1] = data_bits[1]
	sequential_data[2] = data_bits[2]
	sequential_data[3] = data_bits[3]
	sequential_data[4] = parity_bits[0] % 2
	sequential_data[5] = parity_bits[1] % 2
	sequential_data[6] = parity_bits[2] % 2
	
	sequential_codeword = 0
	for index in range(data_count + parity_count):
		sequential_codeword = sequential_codeword + int(sequential_data[index] * (2**index))
	
	return sequential_codeword
	
def Hamming74Decode(input_data):
	valid_codewords = numpy.zeros(16)
	for index in range(16):
		valid_codewords[index] = Hamming74Encode(index)
	
	result = numpy.zeros(16)
	for index in range(16):
		product = int(input_data) ^ int(valid_codewords[index])
		distance = 0
		for bit_index in range(8):
			if product & 1:
				distance = distance + 1
			product = product >> 1
		result[index] = distance
	
	bit_distance = numpy.amin(result)	
	result = numpy.argmin(result) & 0xF
	
	return [result, bit_distance]

encode_table = numpy.zeros(16)

for data in range(16):
	encode_table[data] = Hamming74Encode(data)
	
decode_table = numpy.zeros(128)
distance_table = numpy.zeros(128)
for data in range(128):
	result = Hamming74Decode(data)
	decode_table[data] = result[0]
	distance_table[data] = result[1]
	
#generate a new directory for the report
run_number = 0
print('trying to make a new directory')
while True:
	run_number = run_number + 1
	dirname = f'./run{run_number}/'
	try:
		os.mkdir(dirname)
	except:
		print(dirname + ' exists')
		continue
	break
print(f'made directory {dirname}')
# Generate and save report file
report_file_name = f'run{run_number}_report.txt'
try:
	report_file = open(dirname + report_file_name, 'w+')
except:
	print('Unable to create report file.')
with report_file:
	report_file.write('# Command line: ')
	for argument in sys.argv:
		report_file.write(f'{argument} ')


	report_file.write('\n\n# Hamming(7,4) Encoding Table')
	report_file.write(fo.GenInt16ArrayHexC(f'Hamming74Encode', encode_table, 1))

	report_file.write('\n\n# Hamming(7,4) Decoding Table')
	report_file.write(fo.GenInt16ArrayHexC(f'Hamming74Decode', decode_table, 8))

	report_file.write('\n\n# Hamming(7,4) Distance Table')
	report_file.write(fo.GenInt16ArrayC(f'Hamming74Distance', distance_table, 8))


	report_file.close()