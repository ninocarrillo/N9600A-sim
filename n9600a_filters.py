import struct
import numpy as np

def InitFIR(this):
	this['Buffer'] = np.zeros(len(this['Taps']))
	return this

def UpdateFIR(this, sample):
	this['Buffer'] = this['Buffer'][1:]
	this['Buffer'] = np.append(this['Buffer'], np.array([sample]))
	this['Output'] = np.rint(np.convolve(this['Buffer'], this['Taps'], 'valid') / pow(2, (16 + this['OutputShift'])))
	return this

def InitCanonicIIR(this):
	return this

def MultiplySaturateScale(val1, val2, saturation_bits, scale_bits):
	pos_sat = pow(2,saturation_bits - 1) - 1
	neg_sat = -pow(2,saturation_bits - 1)
	result = val1 * val2
	if result > pos_sat:
		result = pos_sat
	elif result < neg_sat:
		result = neg_sat
	result = result >> scale_bits
	return result

def UpdateCanonicIIR(this, sample):
	# Multiply, scale, accumulate the 'a' factors with the delayed intermediate values
	accumulator = sample
	for index in range(1,this['Order']):
		accumulator += MultiplySaturateScale(this['a'][index], this['w'][index], this['saturation bits'], this['scale bits'])
	# Save this intermediate sum
	this['w'][0] = accumulator
	for index in range(this['Order']):
		accumulator += MultiplySaturateScale(this['b'][index], this['w'][index], this['saturation bits'], this['scale bits'])
	# update the delay registers
	for index in range(this['Order'],1,-1):
		this['w'][index] = this['w'][index - 1]
	this['Output'] = accumulator
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

		
	return
