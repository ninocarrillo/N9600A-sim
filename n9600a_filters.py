import struct
import numpy as np
import sys
import os
import n9600a_strings as strings
import matplotlib as mpl
import matplotlib.pyplot as plt

def GetIIRConfig(config, num, id_string):
	this = {}

	key_string = "iir order"
	try:
		this[f'{key_string}'] = int(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "iir scale bits"
	try:
		this[f'{key_string}'] = int(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "iir saturation bits"
	try:
		this[f'{key_string}'] = int(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "iir b coefs"
	try:
	 	this[f'{key_string}'] = strings.StringToIntArray(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "iir a coefs"
	try:
	 	this[f'{key_string}'] = strings.StringToIntArray(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	return this

def InitFIR(this):
	this['Buffer'] = np.zeros(len(this['Taps']))
	return this

def UpdateFIR(this, sample):
	this['Buffer'] = this['Buffer'][1:]
	this['Buffer'] = np.append(this['Buffer'], np.array([sample]))
	this['Output'] = np.rint(np.convolve(this['Buffer'], this['Taps'], 'valid') / pow(2, (16 + this['OutputShift'])))
	return this

def InitCanonicIIR(this):
	this['NegativeSaturation'] = -pow(2,this['iir saturation bits'] - 1)
	this['PositiveSaturation'] = pow(2,this['iir saturation bits'] - 1) - 1
	this['Output'] = 0
	this['W'] = np.zeros(this['iir order'] + 1)
	return this

def InitIIR(this):
	this['NegativeSaturation'] = -pow(2,this['iir saturation bits'] - 1)
	this['PositiveSaturation'] = pow(2,this['iir saturation bits'] - 1) - 1
	this['Output'] = 0
	this['X'] = np.zeros(this['iir order'] + 1)
	this['Y'] = np.zeros(this['iir order'] + 1)
	return this

def MultiplySaturateScale(val1, val2, pos_sat, neg_sat, scale_bits):
	result = val1 * val2
	if result > pos_sat:
		result = pos_sat
	elif result < neg_sat:
		result = neg_sat
	try:
		result = int(result) >> int(scale_bits)
	except:
		print(f'scale_bits {scale_bits}')
		print(f'result {result}')
		sys.exit(3)
	return result

def MultiplySaturate(val1, val2, pos_sat, neg_sat):
	result = val1 * val2
	if result > pos_sat:
		result = pos_sat
	elif result < neg_sat:
		result = neg_sat
	return result

def UpdateCanonicIIR(this, sample):
	# Multiply, scale, accumulate the 'a' factors with the delayed intermediate values
	accumulator = sample
	for index in range(1,this['iir order'] + 1):
		accumulator += MultiplySaturateScale(this['iir a coefs'][index], this['W'][index], this['PositiveSaturation'], this['NegativeSaturation'], this['iir scale bits'])
	# Save this intermediate sum
	this['W'][0] = accumulator
	accumulator = 0
	for index in range(this['iir order'] + 1):
		accumulator += MultiplySaturateScale(this['iir b coefs'][index], this['W'][index], this['PositiveSaturation'], this['NegativeSaturation'], this['iir scale bits'])
	# update the delay registers
	for index in range(this['iir order'],0,-1):
		this['W'][index] = this['W'][index - 1]
	this['Output'] = accumulator
	return this

def UpdateIIR(this, sample):
	# Update the input delay registers
	for index in range(this['iir order'], 0, -1):
		this['X'][index] = this['X'][index - 1]
	this['X'][0] = sample
	# Calculate the intermediate sum
	v = 0
	for index in range(this['iir order'] + 1):
		v += MultiplySaturateScale(this['X'][index], this['iir b coefs'][index], this['PositiveSaturation'], this['NegativeSaturation'], this['iir scale bits'])
	# Update the output delay registers
	for index in range(this['iir order'], 0, -1):
		this['Y'][index] = this['Y'][index - 1]
	# Calculate the final sum
	for index in range(1, this['iir order'] + 1):
		v += MultiplySaturateScale(this['Y'][index], this['iir a coefs'][index], this['PositiveSaturation'], this['NegativeSaturation'], this['iir scale bits'])
	this['Y'][0] = v
	this['Output'] = v
	return this

def IIRTest(state):
	argv = state['argv']
	config = state['config']

	#generate a new directory for the reports
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


	IIR = InitIIR(GetIIRConfig(config, 1, "IIR "))

	print(IIR)

	# generate a square wave
	count = 300
	mag = 20000
	y = np.zeros(count)
	z = np.zeros(count)

	j = 0
	data = True
	for i in range(count):
		if data == True:
			y[i] = mag
		else:
			y[i] = 0
		j += 1
		if j == 24:
			j = 0
			if data == True:
				data = False
			else:
				data = True
		IIR = UpdateIIR(IIR, y[i])
		z[i] = IIR['Output']

	print(IIR)

	fig1,ax1 = plt.subplots()
	ax1.plot(y)
	ax1.plot(z)
	plt.show()

	return
