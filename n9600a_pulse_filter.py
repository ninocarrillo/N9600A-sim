import numpy as np
import sys
import n9600a_strings as strings
import math

def GetSymbolMapConfig(state):
	argv = state['argv']
	config = state['config']
	this = {}
	id_string = "Symbol Map"

	key_string = "expander"
	try:
		this[key_string] = config[f'{id_string}'][f'{key_string}']
	except:
		print(f'{sys.argv[1]} [{id_string}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "symbol bits"
	try:
		this[key_string] = int(config[f'{id_string}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "symbol map"
	try:
		this[key_string] = strings.StringToIntArray(config[f'{id_string}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	return this

def GetRRCFilterConfig(state):
	argv = state['argv']
	config = state['config']
	this = {}
	id_string = "Pulse Filter"

	key_string = "sample rate"
	try:
		this[f'{key_string}'] = int(config[f'{id_string}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "symbol rate"
	try:
		this[f'{key_string}'] = int(config[f'{id_string}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "rolloff rate"
	try:
		this[f'{key_string}'] = float(config[f'{id_string}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "symbol span"
	try:
		this[f'{key_string}'] = int(config[f'{id_string}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	return this

def InitRRCFilter(this):
	this['Oversample'] = this['sample rate'] // this['symbol rate']
	this['TapCount'] = this['symbol span'] * this['Oversample']
	this['TimeStep'] = 1 / this['sample rate']
	this['SymbolTime'] = 1 / this['symbol rate']
	this['Time'] = np.arange(0, this['TapCount'] * this['TimeStep'], this['TimeStep']) - (this['TapCount'] * this['TimeStep'] / 2) + (this['TimeStep'] / 2)
	this['SymbolTicks'] = np.arange(this['Time'][0] - (this['TimeStep'] / 2), this['Time'][this['TapCount'] - 1], this['SymbolTime'])
	this['Taps'] = np.zeros(this['TapCount'])

	# discontinuity:
	# print(this['TimeStep'] / (4 * this['rolloff rate']))
	index = 0
	try:
		asymptote = this['SymbolTime'] / (4 * this['rolloff rate'])
	except:
		asymptote = False
	for time in this['Time']:
		if math.isclose(time,-asymptote) or math.isclose(time, asymptote):
			numerator = this['rolloff rate'] * ((1 + 2 / np.pi) * np.sin(np.pi/(4 * this['rolloff rate'])) + (1 - (2 / np.pi)) * np.cos(np.pi / (4 * this['rolloff rate'])))
			denominator = this['SymbolTime'] * pow(2, 0.5)
			this['Taps'][index] = numerator / denominator
		else:
			numerator = np.sin(np.pi * time * (1 - this['rolloff rate']) / this['SymbolTime']) + 4 * this['rolloff rate'] * time * np.cos(np.pi * time * (1 + this['rolloff rate']) / this['SymbolTime']) / this['SymbolTime']
			denominator = np.pi * time * (1 - pow(4 * this['rolloff rate'] * time / this['SymbolTime'], 2)) / this['SymbolTime']
			try:
				this['Taps'][index] = numerator / (denominator * this['SymbolTime'])
			except:
				pass
		index += 1
	this['Taps'] = this['Taps'] / np.linalg.norm(this['Taps'])
	return this

def ExpandSampleStream(data, filter):
	bit_count = len(data) * 8
	#print('BitCount', bit_count)
	symbol_count = bit_count // filter['SymbolMap']['symbol bits']
	#print('SymbolCount', symbol_count)
	sample_count = symbol_count * filter['Oversample']
	#print('SampleCount', sample_count)
	flush_count = filter['symbol span'] * filter['Oversample']
	samples = np.zeros(sample_count + flush_count)
	sample_index = 0
	symbols_per_byte = 8 // filter['SymbolMap']['symbol bits']
	if filter['SymbolMap']['expander'] == 'impulse':
		for byte in data:
			byte = int(byte)
			for byte_index in range(symbols_per_byte):
				symbol = np.bitwise_and(byte, 3)
				byte = np.right_shift(byte, filter['SymbolMap']['symbol bits'])
				samples[sample_index] = filter['SymbolMap']['symbol map'][symbol]
				sample_index += 1
				for extend_index in range(filter['Oversample'] - 1):
					samples[sample_index] = 0
					sample_index += 1
	elif filter['SymbolMap']['expander'] == 'step':
		for byte in data:
			byte = int(byte)
			for byte_index in range(symbols_per_byte):
				symbol = np.bitwise_and(byte, 3)
				byte = np.right_shift(byte, filter['SymbolMap']['symbol bits'])
				samples[sample_index] = filter['SymbolMap']['symbol map'][symbol]
				sample_index += 1
				for extend_index in range(filter['Oversample'] - 1):
					samples[sample_index] = filter['SymbolMap']['symbol map'][symbol]
					sample_index += 1
	for flush_index in range(flush_count):
		samples[sample_index] = 0
		sample_index += 1
	return samples

def GenEyeData(samples, oversample, delay):
	samples = samples[delay:]
	samples = samples[:delay]
	y_data = np.zeros(len(samples))
	x_data = np.zeros(len(samples))
	index = 0
	offset = oversample // 2
	for sample in samples:
		for sub_index in range(oversample):
			try:
				x_data[index] = (sub_index + offset) % oversample
				y_data[index] = samples[index]
				index += 1
			except:
				pass

	data = {}
	data['x'] = x_data
	data['y'] = y_data
	return(data)
