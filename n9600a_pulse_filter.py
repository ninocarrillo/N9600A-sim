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

def GetAGCConfig(state):
	this = {}
	argv = state['argv']
	config = state['config']
	id_string = "AGC"

	key_string = "agc attack rate"
	try:
		this[key_string] = int(config[f'{id_string}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "agc sustain period"
	try:
		this[key_string] = int(config[f'{id_string}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "agc decay rate"
	try:
		this[key_string] = int(config[f'{id_string}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "agc high thresh"
	try:
		this[key_string] = int(config[f'{id_string}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	id_string = "Decimator"
	key_string = "decimation"
	try:
		this[key_string] = int(config[f'{id_string}'][f'{key_string}'])
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

	key_string = "bit count"
	try:
		this[f'{key_string}'] = int(config[f'{id_string}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "undersample"
	try:
		this[f'{key_string}'] = int(config[f'{id_string}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "amplitude"
	try:
		this[f'{key_string}'] = int(config[f'{id_string}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)
		
	key_string = "window"
	try:
		this[f'{key_string}'] = config[f'{id_string}'][f'{key_string}']
	except:
		print(f'{sys.argv[1]} [{id_string}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	return this

def GetGaussFilterConfig(state):
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

	# key_string = "pulse high"
	# try:
	# 	this[f'{key_string}'] = float(config[f'{id_string}'][f'{key_string}'])
	# except:
	# 	print(f'{sys.argv[1]} [{id_string}] \'{key_string}\' is missing or invalid')
	# 	sys.exit(-2)
	#
	# key_string = "pulse low"
	# try:
	# 	this[f'{key_string}'] = float(config[f'{id_string}'][f'{key_string}'])
	# except:
	# 	print(f'{sys.argv[1]} [{id_string}] \'{key_string}\' is missing or invalid')
	# 	sys.exit(-2)

	key_string = "BT"
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

	key_string = "amplitude"
	try:
		this[f'{key_string}'] = int(config[f'{id_string}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "undersample"
	try:
		this[f'{key_string}'] = int(config[f'{id_string}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	return this

def InitFilterDecimator(filter_decimator):
	filter_decimator['FilterBuffer'] = np.zeros(len(filter_decimator['Filter']))
	filter_decimator['DataBuffer'] = np.array([])
	filter_decimator['FilterShift'] = -5
	filter_decimator['DecimationCounter'] = 0
	filter_decimator['NewSample'] = 0
	filter_decimator['PeakDetector'] = {'AttackRate':filter_decimator['agc attack rate'], 'SustainPeriod': filter_decimator['agc sustain period'], 'DecayRate':filter_decimator['agc decay rate'], 'SustainCount':0, 'Envelope':0}
	filter_decimator['OutputSampleRate'] = filter_decimator['InputSampleRate'] // filter_decimator['decimation']
	return filter_decimator

def InitGaussFilter(this):
	this['Oversample'] = this['sample rate'] // (this['symbol rate'] * this['undersample'])
	this['TapCount'] = this['symbol span'] * this['Oversample']
	this['TimeStep'] = 1 / this['sample rate']
	this['SymbolTime'] = 1 / this['symbol rate']
	# generate normalized time
	this['Time'] = np.arange(0, this['symbol span'], 1 / this['Oversample']) - (this['symbol span'] / 2)
	this['SymbolTicks'] = np.arange(-this['symbol span'] / 2, this['symbol span'] / 2, 1)
	this['Taps'] = np.zeros(this['TapCount'])

	alpha = pow(np.log(2) / 2, 0.5) / this['BT']
	print('alpha ', alpha)
	index = 0
	for time in this['Time']:
		#numerator = np.exp(-this['alpha'] * pow(time, 2))
		numerator = np.exp(-pow(np.pi,2) * pow(time, 2) / pow(alpha, 2))
		try:
			this['Taps'][index] = numerator
		except:
			pass
		index += 1
	this['Taps'] = this['Taps'] / np.linalg.norm(this['Taps'])
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
	this['RC'] = np.convolve(this['Taps'], this['Taps'], 'same')
	
	this['FilterWindow'] = np.zeros(this['TapCount'])
	
	N = this['TapCount'] - 1
	if this['window'] == 'hann':
		for index in range(this['TapCount']):
			this['FilterWindow'][index] = np.power(np.sin(np.pi * index / N),2)
	elif this['window'] == 'rect':
		for index in range(this['TapCount']):
			this['FilterWindow'][index] = 1
	elif this['window'] == 'blackmann':
		a0 = 0.355768
		a1 = 0.487396
		a2 = 0.144232
		a3 = 0.012604
		for index in range(this['TapCount']):
			this['FilterWindow'][index] = a0 - (a1 * np.cos(2 * np.pi * index / N)) + (a2 * np.cos(4 * np.pi * index / N)) - (a3 * np.cos(6 * np.pi * index / N))
	elif this['window'] == 'blackmann-harris':
		a0 = 0.35875
		a1 = 0.48829
		a2 = 0.14128
		a3 = 0.01168
		for index in range(this['TapCount']):
			this['FilterWindow'][index] = a0 - (a1 * np.cos(2 * np.pi * index / N)) + (a2 * np.cos(4 * np.pi * index / N)) - (a3 * np.cos(6 * np.pi * index / N))
	elif this['window'] == 'flattop':
		a0 = 0.21557895
		a1 = 0.41663158
		a2 = 0.277263158
		a3 = 0.083578947
		a4 = 0.006947368
		for index in range(this['TapCount']):
			this['FilterWindow'][index] = a0 - (a1 * np.cos(2 * np.pi * index / N)) + (a2 * np.cos(4 * np.pi * index / N)) - (a3 * np.cos(6 * np.pi * index / N)) + (a4  * np.cos(8 * np.pi * index / N))
	

	this['Taps'] = np.multiply(this['Taps'], this['FilterWindow'])
	this['WindowedRC'] = np.convolve(this['Taps'], this['Taps'], 'same')
	this['WindowedRC'] = this['WindowedRC'] * max(this['Taps']) / max(this['WindowedRC'])
	return this

def BytesToSymbols(data, filter):
	bit_count = len(data) * 8
	symbol_count = bit_count // filter['SymbolMap']['symbol bits']
	sample_count = symbol_count
	samples = np.zeros(sample_count)
	sample_index = 0
	symbols_per_byte = 8 // filter['SymbolMap']['symbol bits']
	symbol_mask = int(pow(2,filter['SymbolMap']['symbol bits']) - 1)
	for byte in data:
		byte = int(byte)
		for byte_index in range(symbols_per_byte):
			symbol = np.bitwise_and(byte, symbol_mask)
			byte = np.right_shift(byte, filter['SymbolMap']['symbol bits'])
			samples[sample_index] = filter['SymbolMap']['symbol map'][symbol]
			sample_index += 1
	return samples

def GenFilterPhases(filter):
	filter['PhaseTaps'] = np.zeros(len(filter['Taps']))
	for i0 in range(filter['Oversample']):
		i0i = (filter['Oversample'] - i0) - 1
		for i1 in range (filter['symbol span']):
			try:
				filter['PhaseTaps'][i0 * filter['symbol span'] + i1] = filter['Taps'][i0 + (i1 * filter['Oversample'])]
			except:
				print(i0, i1)
	return filter

def ImpulseOversample(data, filter):
	undersample_fir = np.zeros([filter['Oversample'], filter['symbol span']])
	for i0 in range(filter['Oversample']):
		undersample_fir[i0] = filter['Taps'][i0::filter['Oversample']]
	# data is a list of 8-bit integer values
	symbol_stream = np.zeros(len(data) * 8 // filter['SymbolMap']['symbol bits'])
	x3 = 0
	for x1 in data:
		x1 = int(x1)
		for x2 in range(8 // filter['SymbolMap']['symbol bits']):
			symbol_stream[x3] = filter['SymbolMap']['symbol map'][np.right_shift(x1, (8 - filter['SymbolMap']['symbol bits']))]
			x1 = np.left_shift(x1, filter['SymbolMap']['symbol bits'])
			x1 = np.bitwise_and(x1, 255)
			x3 += 1
	sample_buffer = np.zeros(filter['symbol span'])
	oversample_index = 0
	sample_stream = np.zeros(len(symbol_stream) * filter['Oversample'])
	sample_index = 0
	for x1 in symbol_stream:
		for x2 in range(filter['symbol span'] - 1):
			sample_buffer[x2] = sample_buffer[x2 + 1]
		sample_buffer[filter['symbol span']-1] = x1
		for x3 in range(filter['Oversample']):
			sample_stream[sample_index] = np.convolve(sample_buffer, undersample_fir[x3], mode='valid')
			sample_index += 1

	return sample_stream

def ImpulseOversample2(data, filter):
	YMem = np.zeros([filter['Oversample'], filter['symbol span']])
	for i0 in range(filter['Oversample']):
		YMem[i0] = filter['Taps'][i0::filter['Oversample']]
	# data is a list of 8-bit integer values
	symbol_stream = np.zeros(len(data) * 8 // filter['SymbolMap']['symbol bits'])
	x3 = 0
	for x1 in data:
		x1 = int(x1)
		for x2 in range(8 // filter['SymbolMap']['symbol bits']):
			symbol_stream[x3] = filter['SymbolMap']['symbol map'][np.right_shift(x1, (8 - filter['SymbolMap']['symbol bits']))]
			x1 = np.left_shift(x1, filter['SymbolMap']['symbol bits'])
			x1 = np.bitwise_and(x1, 255)
			x3 += 1
	sample_buffer = np.zeros(filter['symbol span'])
	sample_stream = np.zeros(len(symbol_stream) * filter['Oversample'])
	sample_index = 0
	YMemIndex = 0
	x0 = 0
	while x0 < len(symbol_stream):
		if YMemIndex == 0:
			for x2 in range(filter['symbol span']-1):
				sample_buffer[x2] = sample_buffer[x2 + 1]
			sample_buffer[filter['symbol span']-1] = symbol_stream[x0]
			x0 += 1
		sample_stream[sample_index] = np.convolve(sample_buffer, YMem[YMemIndex], mode='valid')
		sample_index += 1
		YMemIndex += 1
		if YMemIndex >= filter['Oversample']:
			YMemIndex = 0

	return sample_stream

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
				symbol = np.right_shift(byte, (8 - filter['SymbolMap']['symbol bits']))
				byte = np.left_shift(byte, filter['SymbolMap']['symbol bits'])
				byte = np.bitwise_and(byte, 255)
				samples[sample_index] = filter['SymbolMap']['symbol map'][symbol]
				sample_index += 1
				for extend_index in range(filter['Oversample'] - 1):
					samples[sample_index] = 0
					sample_index += 1
	elif filter['SymbolMap']['expander'] == 'step':
		for byte in data:
			byte = int(byte)
			for byte_index in range(symbols_per_byte):
				symbol = np.bitwise_and(byte, symbol_mask)
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

def GenEyeData2(samples, oversample, delay):
	trace_count = int(np.floor((len(samples) - delay) / oversample))
	#print(trace_count)
	samples = samples[delay:]
	traces = np.zeros([oversample, trace_count])
	index = 0
	for trace_number in range(trace_count):
		for sub_index in range(oversample):
			try:
				traces[sub_index,trace_number] = samples[index]
			except:
				pass
			index += 1
	return traces

def GenEyeData(samples, oversample, delay):
	#samples = samples[delay:]
	#samples = samples[:delay]
	y_data = np.zeros(len(samples))
	x_data = np.zeros(len(samples))
	index = 0
	offset = oversample / 2
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

def PeakDetect2(signal_value, detector):
	signal_value = abs(signal_value)
	delta = signal_value - detector['Envelope']
	if delta > 0:
		detector['Envelope'] = detector['Envelope'] + np.rint((delta * detector['AttackRate'] / 128))
		if detector['Envelope'] > signal_value:
			detector['Envelope'] = signal_value
		detector['SustainCount'] = 0
	if detector['SustainCount'] >= detector['SustainPeriod']:
		detector['Envelope'] = detector['Envelope'] - detector['DecayRate']
		if detector['Envelope'] < 0:
			detector['Envelope'] = 0
			detector['SustainCount'] = 0
	detector['SustainCount'] = detector['SustainCount'] + 1
	return detector

def FilterDecimate2(filter):
	filter['FilterBuffer'] = np.rint(np.convolve(filter['FilterBuffer'], filter['Filter'], 'valid'))
	filter['FilterBuffer'] = filter['FilterBuffer'][::filter['DecimationRate']]
	filter['EnvelopeBuffer'] = np.zeros(len(filter['FilterBuffer']))
	filter['ThreshBuffer'] = np.zeros(len(filter['FilterBuffer']))
	index = 0
	for data in filter['FilterBuffer']:
		filter['EnvelopeBuffer'][index] = filter['PeakDetector']['Envelope']
		filter['ThreshBuffer'][index] = np.rint(filter['PeakDetector']['Envelope'] * filter['agc high thresh'] / 128)
		data = data // pow(2, (16 + filter['FilterShift']))
		filter['FilterBuffer'][index] = data
		index += 1
		if filter['InputAGCEnabled'] == True:
			filter['PeakDetector'] = PeakDetect2(data, filter['PeakDetector'])
			filter['GainChange'] = 0
			if filter['PeakDetector']['Envelope'] > 24576:
				filter['FilterShift'] = filter['FilterShift'] + 1
				filter['PeakDetector']['Envelope'] = filter['PeakDetector']['Envelope'] / 2
				if filter['FilterShift'] > 16:
					filter['FilterShift'] = 16
			if filter['PeakDetector']['Envelope'] < 8192:
				filter['FilterShift'] = filter['FilterShift'] - 1
				filter['PeakDetector']['Envelope'] = filter['PeakDetector']['Envelope'] * 2
				if filter['FilterShift'] < -16:
					filter['FilterShift'] = -16
	filter['FilterBuffer'] = np.clip(filter['FilterBuffer'], -32768, 32767)
	return filter

def GenPulseFilterPatterns(this):
	bit_width = this['symbol span'] * this['SymbolMap']['symbol bits']
	symbol_count = (2**bit_width) // this['SymbolMap']['symbol bits']
	samples_per_symbol = this['Oversample']
	symbol_mask = np.power(2, this['SymbolMap']['symbol bits']) - 1
	symbol_tap = symbol_mask << ((this['symbol span'] - 1) * this['SymbolMap']['symbol bits'])
	shift_mask = np.power(2,this['symbol span'] * this['SymbolMap']['symbol bits'])
	this['FilterPatterns'] = np.zeros(symbol_count * samples_per_symbol)
	for i0 in range(symbol_count):
		#i0 is the symbol to be factored
		y = np.zeros(samples_per_symbol * this['symbol span']*this['SymbolMap']['symbol bits'])
		factor_me = i0
		# bit-reverse i0
		shift_register = 0
		for ix in range(bit_width):
			shift_register <<= 1
			if factor_me & 1:
				shift_register |= 1
			factor_me >>= 1
		factor_me = shift_register
		bit_index = 0
		for i2 in range(this['symbol span'] * this['SymbolMap']['symbol bits']):
			shift_register <<= this['SymbolMap']['symbol bits']
			factor = factor_me & symbol_mask
			factor_me = factor_me >> this['SymbolMap']['symbol bits']
			level = this['SymbolMap']['symbol map'][factor]

			for i3 in range(samples_per_symbol):
				if i3 == 0:
					y[(i2 * samples_per_symbol) + i3] = level
				else:
					if this['SymbolMap']['expander'] == 'step':
						y[(i2 * samples_per_symbol) + i3] = level
					else:
						y[(i2 * samples_per_symbol) + i3] = 0

			z = np.convolve(y, this['Taps'], 'full')
			# trim the invalid results:
			z = z[len(this['Taps'])//2:]
			z = z[:samples_per_symbol*this['symbol span']]

			# select the center symbol length
			x_offset = ((this['symbol span'] * samples_per_symbol) // 2) - (samples_per_symbol // 2)
			xz = np.arange(x_offset, x_offset + samples_per_symbol)
			z = z[x_offset:x_offset+samples_per_symbol]
			i4 = 0
			for x in z:
				this['FilterPatterns'][((i0) * samples_per_symbol) + i4] = x
				i4 += 1
	try:
		if (this['amplitude']) != 0:
			this['FilterPatterns'] = np.rint(this['amplitude'] * this['FilterPatterns'] / max(abs(this['FilterPatterns'])))
	except:
		pass
	return this
