import numpy as np
import sys
import n9600a_strings as strings

def GetInputFilterConfig(state):
	argv = state['argv']
	config = state['config']
	this = {}
	try:
	 	this['Filter'] = strings.StringToIntArray(config['Input Filter']['taps'])
	except:
		print(f'{argv[1]} [Input Filter] \'taps\' is missing or invalid')
		sys.exit(-2)
	# Get input filter AGC option:
	try:
		this['InputAGCEnabled'] = config['Input Filter'].getboolean('agc enabled')
	except:
		print(f'{argv[1]} [Input Filter] \'agc enabled\' is missing or invalid')
		sys.exit(-2)

	# Get Input Filter AGC Attack Rate:
	try:
		this['InputAGCAttackRate'] = int(float(config['Input Filter']['agc attack rate']))
	except:
		print(f'{argv[1]} [Input Filter] \'agc attack rate\' is missing or invalid')
		sys.exit(-2)

	# Get Input Filter AGC Sustain Period:
	try:
		this['InputAGCSustainPeriod'] = int(float(config['Input Filter']['agc sustain period']))
	except:
		print(f'{argv[1]} [Input Filter] \'agc sustain period\' is missing or invalid')
		sys.exit(-2)

	# Get Input Filter AGC Decay Rate:
	try:
		this['InputAGCDecayRate'] = int(float(config['Input Filter']['agc decay rate']))
	except:
		print(f'{argv[1]} [Input Filter] \'agc decay rate\' is missing or invalid')
		sys.exit(-2)

	# Get Input Filter ADC bit count:
	try:
		this['InputBitCount'] = int(float(config['Input Filter']['input bit count']))
	except:
		print(f'{argv[1]} [Input Filter] \'input bit count\' is missing or invalid')
		sys.exit(-2)
	# Get input filter Sample Rate:
	try:
		this['InputSampleRate'] = int(float(config['Input Filter']['sample rate']))
	except:
		print(f'{argv[1]} [Input Filter] \'sample rate\' is missing or invalid')
		sys.exit(-2)

	# Get input filter Decimation Rate:
	try:
		this['DecimationRate'] = int(config['Input Filter']['decimation'])
	except:
		print(f'{argv[1]} [Input Filter] \'decimation\' is missing or invalid')
		sys.exit(-2)

	try:
		this['MaxAGCShift'] = int(float(config['Input Filter']['agc max shift']))
	except:
		print(f'{argv[1]} [Input Filter] \'agc max shift\' is missing or invalid')
		sys.exit(-2)

	try:
		this['MinAGCShift'] = int(float(config['Input Filter']['agc min shift']))
	except:
		print(f'{argv[1]} [Input Filter] \'agc min shift\' is missing or invalid')
		sys.exit(-2)

	return this

def GetInputFilterConfig2(state):
	argv = state['argv']
	config = state['config']
	this = {}
	try:
	 	this['Filter'] = strings.StringToIntArray(config['Input Filter']['taps'])
	except:
		print(f'{argv[1]} [Input Filter] \'taps\' is missing or invalid')
		sys.exit(-2)
	# Get input filter AGC option:
	try:
		this['InputAGCEnabled'] = config['Input Filter'].getboolean('agc enabled')
	except:
		print(f'{argv[1]} [Input Filter] \'agc enabled\' is missing or invalid')
		sys.exit(-2)

	# Get Input Filter AGC Attack Rate:
	try:
		this['InputAGCAttackRate'] = int(config['Input Filter']['agc attack rate'])
	except:
		print(f'{argv[1]} [Input Filter] \'agc attack rate\' is missing or invalid')
		sys.exit(-2)

	# Get Input Filter AGC Sustain Period:
	try:
		this['InputAGCSustainPeriod'] = int(config['Input Filter']['agc sustain period'])
	except:
		print(f'{argv[1]} [Input Filter] \'agc sustain period\' is missing or invalid')
		sys.exit(-2)

	# Get Input Filter AGC Decay Rate:
	try:
		this['InputAGCDecayRate'] = int(config['Input Filter']['agc decay rate'])
	except:
		print(f'{argv[1]} [Input Filter] \'agc decay rate\' is missing or invalid')
		sys.exit(-2)

	# Get Input Filter Threshold:
	try:
		this['agc high thresh'] = int(config['Input Filter']['agc high thresh'])
	except:
		print(f'{argv[1]} [Input Filter] \'agc high thresh\' is missing or invalid')
		sys.exit(-2)

	# Get Input Filter ADC bit count:
	try:
		this['InputBitCount'] = int(config['Input Filter']['input bit count'])
	except:
		print(f'{argv[1]} [Input Filter] \'input bit count\' is missing or invalid')
		sys.exit(-2)
	# Get input filter Sample Rate:
	try:
		this['InputSampleRate'] = int(config['Input Filter']['sample rate'])
	except:
		print(f'{argv[1]} [Input Filter] \'sample rate\' is missing or invalid')
		sys.exit(-2)

	# Get input filter Decimation Rate:
	try:
		this['DecimationRate'] = int(config['Input Filter']['decimation'])
	except:
		print(f'{argv[1]} [Input Filter] \'decimation\' is missing or invalid')
		sys.exit(-2)
	return this

def GetInputFilterConfig3(state):
	argv = state['argv']
	config = state['config']
	this = {}
	try:
	 	this['Filter'] = strings.StringToIntArray(config['Input Filter']['taps'])
	except:
		print(f'{argv[1]} [Input Filter] \'taps\' is missing or invalid')
		sys.exit(-2)
	# Get input filter AGC option:
	try:
		this['InputAGCEnabled'] = config['Input Filter'].getboolean('agc enabled')
	except:
		print(f'{argv[1]} [Input Filter] \'agc enabled\' is missing or invalid')
		sys.exit(-2)
	key_string = "agc iir order"
	try:
		this[f'{key_string}'] = int(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "agc iir scale bits"
	try:
		this[f'{key_string}'] = int(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "agc iir saturation bits"
	try:
		this[f'{key_string}'] = int(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "agc iir b coefs"
	try:
	 	this[f'{key_string}'] = strings.StringToIntArray(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "agc iir a coefs"
	try:
	 	this[f'{key_string}'] = strings.StringToIntArray(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	# Get Input Filter ADC bit count:
	try:
		this['InputBitCount'] = int(config['Input Filter']['input bit count'])
	except:
		print(f'{argv[1]} [Input Filter] \'input bit count\' is missing or invalid')
		sys.exit(-2)
	# Get input filter Sample Rate:
	try:
		this['InputSampleRate'] = int(config['Input Filter']['sample rate'])
	except:
		print(f'{argv[1]} [Input Filter] \'sample rate\' is missing or invalid')
		sys.exit(-2)

	# Get input filter Decimation Rate:
	try:
		this['DecimationRate'] = int(config['Input Filter']['decimation'])
	except:
		print(f'{argv[1]} [Input Filter] \'decimation\' is missing or invalid')
		sys.exit(-2)
	return this

def GetInputFilterConfig5(state):
	argv = state['argv']
	config = state['config']
	this = {}
	# Get input filter AGC option:
	try:
		this['InputAGCEnabled'] = config['Input Filter'].getboolean('agc enabled')
	except:
		print(f'{argv[1]} [Input Filter] \'agc enabled\' is missing or invalid')
		sys.exit(-2)

	# Get Input Filter AGC Attack Rate:
	try:
		this['InputAGCAttackRate'] = int(float(config['Input Filter']['agc attack rate']))
	except:
		print(f'{argv[1]} [Input Filter] \'agc attack rate\' is missing or invalid')
		sys.exit(-2)

	# Get Input Filter AGC Sustain Period:
	try:
		this['InputAGCSustainPeriod'] = int(float(config['Input Filter']['agc sustain period']))
	except:
		print(f'{argv[1]} [Input Filter] \'agc sustain period\' is missing or invalid')
		sys.exit(-2)

	# Get Input Filter AGC Decay Rate:
	try:
		this['InputAGCDecayRate'] = int(float(config['Input Filter']['agc decay rate']))
	except:
		print(f'{argv[1]} [Input Filter] \'agc decay rate\' is missing or invalid')
		sys.exit(-2)

	# Get Input Filter ADC bit count:
	try:
		this['InputBitCount'] = int(float(config['Input Filter']['input bit count']))
	except:
		print(f'{argv[1]} [Input Filter] \'input bit count\' is missing or invalid')
		sys.exit(-2)
	# Get input filter Sample Rate:
	try:
		this['InputSampleRate'] = int(float(config['Input Filter']['sample rate']))
	except:
		print(f'{argv[1]} [Input Filter] \'sample rate\' is missing or invalid')
		sys.exit(-2)

	# Get input filter Decimation Rate:
	try:
		this['DecimationRate'] = int(float(config['Input Filter']['decimation']))
	except:
		print(f'{argv[1]} [Input Filter] \'decimation\' is missing or invalid')
		sys.exit(-2)


	try:
		this['MaxAGCShift'] = int(float(config['Input Filter']['agc max shift']))
	except:
		print(f'{argv[1]} [Input Filter] \'agc max shift\' is missing or invalid')
		sys.exit(-2)


	try:
		this['MinAGCShift'] = int(float(config['Input Filter']['agc min shift']))
	except:
		print(f'{argv[1]} [Input Filter] \'agc min shift\' is missing or invalid')
		sys.exit(-2)

	return this
