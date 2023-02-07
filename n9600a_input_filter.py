import numpy as np
import sys
import n9600a_strings as strings

def GetInputFilterConfig(state):
	argv = state['argv']
	config = state['config']
	FilterDecimator = {}
	try:
	 	FilterDecimator['Filter'] = strings.StringToIntArray(config['Input Filter']['taps'])
	except:
		print(f'{argv[1]} [Input Filter] \'taps\' is missing or invalid')
		sys.exit(-2)
	# Get input filter AGC option:
	try:
		FilterDecimator['InputAGCEnabled'] = config['Input Filter'].getboolean('agc enabled')
	except:
		print(f'{argv[1]} [Input Filter] \'agc enabled\' is missing or invalid')
		sys.exit(-2)

	# Get Input Filter AGC Attack Rate:
	try:
		FilterDecimator['InputAGCAttackRate'] = int(config['Input Filter']['agc attack rate'])
	except:
		print(f'{argv[1]} [Input Filter] \'agc attack rate\' is missing or invalid')
		sys.exit(-2)

	# Get Input Filter AGC Sustain Period:
	try:
		FilterDecimator['InputAGCSustainPeriod'] = int(config['Input Filter']['agc sustain period'])
	except:
		print(f'{argv[1]} [Input Filter] \'agc sustain period\' is missing or invalid')
		sys.exit(-2)

	# Get Input Filter AGC Decay Rate:
	try:
		FilterDecimator['InputAGCDecayRate'] = int(config['Input Filter']['agc decay rate'])
	except:
		print(f'{argv[1]} [Input Filter] \'agc decay rate\' is missing or invalid')
		sys.exit(-2)

	# Get Input Filter ADC bit count:
	try:
		FilterDecimator['InputBitCount'] = int(config['Input Filter']['input bit count'])
	except:
		print(f'{argv[1]} [Input Filter] \'input bit count\' is missing or invalid')
		sys.exit(-2)
	# Get input filter Sample Rate:
	try:
		FilterDecimator['InputSampleRate'] = int(config['Input Filter']['sample rate'])
	except:
		print(f'{argv[1]} [Input Filter] \'sample rate\' is missing or invalid')
		sys.exit(-2)

	# Get input filter Decimation Rate:
	try:
		FilterDecimator['DecimationRate'] = int(config['Input Filter']['decimation'])
	except:
		print(f'{argv[1]} [Input Filter] \'decimation\' is missing or invalid')
		sys.exit(-2)
	return FilterDecimator
