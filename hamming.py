import numpy
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
	
	parity_bits[0] = data_bits[0] + data_bits[1]                + data_bits[3] + even
	parity_bits[1] = data_bits[0]                + data_bits[2] + data_bits[3] + even
	parity_bits[2] =                data_bits[1] + data_bits[2] + data_bits[3] + even
	
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
	
	return sequential_codeword, interleaved_codeword


for data in range(16):
	codewords = Hamming74Encode(data)
	print(f's: {hex(codewords[0])}    i:{hex(codewords[1])}')