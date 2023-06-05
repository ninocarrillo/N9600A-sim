import numpy
import math

def doit(config, audio):
	if int(config['sample_rate']) != config['SampleRate']:
		print("Warning! Sample rate mismatch.")
	this = {}
	config['Oversample'] = int(config['sample_rate']) // int(config['symbol_rate'])
	config['TapCount'] = int(config['symbol_span']) * config['Oversample']
	config['TimeStep'] = 1 / int(config['sample_rate'])
	config['SymbolTime'] = 1 / int(config['symbol_rate'])
	config['Time'] = numpy.arange(0, config['TapCount'] * config['TimeStep'], config['TimeStep']) - (config['TapCount'] * config['TimeStep'] / 2) + (config['TimeStep'] / 2)
	config['SymbolTicks'] = numpy.arange(config['Time'][0] - (config['TimeStep'] / 2), config['Time'][config['TapCount'] - 1], config['SymbolTime'])
	config['Taps'] = numpy.zeros(config['TapCount'])
	index = 0
	try:
		asymptote = config['SymbolTime'] / (4 * float(float(config['rolloff_rate'])))
	except:
		asymptote = False
	for time in config['Time']:
		if math.isclose(time,-asymptote) or math.isclose(time, asymptote):
			numerator = float(float(config['rolloff_rate'])) * ((1 + 2 / numpy.pi) * numpy.sin(numpy.pi/(4 * float(float(config['rolloff_rate'])))) + (1 - (2 / numpy.pi)) * numpy.cos(numpy.pi / (4 * float(float(config['rolloff_rate'])))))
			denominator = config['SymbolTime'] * pow(2, 0.5)
			config['Taps'][index] = numerator / denominator
		else:
			numerator = numpy.sin(numpy.pi * time * (1 - float(float(config['rolloff_rate']))) / config['SymbolTime']) + 4 * float(config['rolloff_rate']) * time * numpy.cos(numpy.pi * time * (1 + float(float(config['rolloff_rate']))) / config['SymbolTime']) / config['SymbolTime']
			denominator = numpy.pi * time * (1 - pow(4 * float(float(config['rolloff_rate'])) * time / config['SymbolTime'], 2)) / config['SymbolTime']
			try:
				config['Taps'][index] = numerator / (denominator * config['SymbolTime'])
			except:
				pass
		index += 1
	config['Taps'] = config['Taps'] / numpy.linalg.norm(config['Taps'])

	audio = numpy.convolve(audio, numpy.rint(config['Taps'] * int(config['amplitude'])), 'valid') // 65536


	return audio
