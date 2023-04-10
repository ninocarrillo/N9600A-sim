import numpy as np
import sys
import n9600a_strings as strings

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
	this['TapCount'] = this['symbol span'] * this['sample rate'] // this['symbol rate']
	this['TimeStep'] = 1 / this['sample rate']
	this['SymbolTime'] = 1 / this['symbol rate']
	this['Time'] = np.arange(0, this['TapCount'] * this['TimeStep'], this['TimeStep']) - (this['TapCount'] * this['TimeStep'] / 2) + (this['TimeStep'] / 2)
	this['SymbolTicks'] = np.arange(this['Time'][0] - (this['TimeStep'] / 2), this['Time'][this['TapCount'] - 1], this['SymbolTime'])
	this['Taps'] = np.zeros(this['TapCount'])

	# discontinuity:
	# print(this['TimeStep'] / (4 * this['rolloff rate']))
	index = 0
	for time in this['Time']:
		if time == -this['SymbolTime'] / (4 * this['rolloff rate']) or time == this['SymbolTime'] / (4 * this['rolloff rate']):
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

	return this