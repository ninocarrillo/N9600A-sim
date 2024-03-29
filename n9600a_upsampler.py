import numpy as np
import sys
import n9600a_strings as strings

def GetUpsamplerConfig(config):
	this = {}

	section = 'Upsampler'

	key = 'enabled'
	try:
		this[f'{key}'] = config[f'{section}'].getboolean(f'{key}')
	except:
		print(f'{sys.argv[1]} [{section}] \'{key}\' is missing or invalid')
		sys.exit(-2)

	key = 'rate'
	try:
		this[f'{key}'] = int(config[f'{section}'][f'{key}'])
	except:
		print(f'{sys.argv[1]} [{section}] \'{key}\' is missing or invalid')
		sys.exit(-2)

	key = 'taps'
	try:
	 	this[f'{key}'] = strings.StringToIntArray(config[f'{section}'][f'{key}'])
	except:
		print(f'{sys.argv[1]} [{section}] \'{key}\' is missing or invalid')
		sys.exit(-2)

	return this

def InitUpsampler(this):
	this['OutputSampleRate'] = this['InputSampleRate'] * this['rate']
	return this

def Upsample(this):
	this['OutputBuffer'] = np.zeros(len(this['InputBuffer']) * this['rate'])
	this['OutputBuffer'][::this['rate']] = this['InputBuffer']
	this['OutputBuffer'] = np.rint(np.convolve(this['OutputBuffer'], this['taps'], 'valid')) // 32768
	return this

def UpsampleSlow(this):
	this['OutputBuffer'] = np.zeros(len(this['InputBuffer']) * this['rate'])
	input_sample_count = len(this['InputBuffer'])
	input_sample_index = 0
	filter_tap_index = 0
	#while input_sample_index < input_sample_count:

	return this
