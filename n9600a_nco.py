import sys
import configparser
import struct
import scipy.io.wavfile
import numpy as np
import os
import n9600a_progdemod as demod
import format_output as fo
import n9600a_strings as strings
import n9600a_input_filter as input_filter
import n9600a_nco as nco

def GetNCOConfig(config, num, id_string):
	this = {}
	key_string = "enabled"
	try:
		this[f'{key_string}'] = config[f'{id_string}{num}'].getboolean(f'{key_string}')
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "design sample rate"
	try:
		this[f'{key_string}'] = int(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "wavetable size"
	try:
		this[f'{key_string}'] = int(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "set frequency"
	try:
		this[f'{key_string}'] = int(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "amplitude bits"
	try:
		this[f'{key_string}'] = int(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "phase dither bits"
	try:
		this[f'{key_string}'] = int(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)
	return this

def InitNCO(this):
	this['NormalizationFactor'] = int(np.ceil(this['design sample rate'] / this['wavetable size']))
	this['InPhase'] = 0
	this['QuadraturePhaseOffset'] = int(np.rint(this['design sample rate'] / 4))
	this['Control'] = 0
	this['Amplitude'] = pow(2,this['amplitude bits'] - 1) - 1
	this['WaveTable'] = np.zeros(this['wavetable size'])
	for i in range(this['wavetable size']):
		this['WaveTable'][i] = np.rint(this['Amplitude'] * np.sin(i * 2 * np.pi / this['wavetable size']))
	return this

def UpdateNCO(this):
	if this['phase dither bits'] > 0:
		this['Dither'] = np.floor(np.random.rand() * pow(2,this['phase dither bits'])) - pow(2,this['phase dither bits'] - 1)
	else:
		this['Dither'] = 0

	this['InPhase'] += this['set frequency'] + this['Control']
	while this['InPhase'] < 0:
		this['InPhase'] += this['design sample rate']
	while this['InPhase'] >= this['design sample rate']:
		this['InPhase'] -= this['design sample rate']

	inphase = this['InPhase'] + this['Dither']
	quadraturephase = inphase + this['QuadraturePhaseOffset']

	while inphase < 0:
		inphase += this['design sample rate']
	while inphase >= this['design sample rate']:
		inphase -= this['design sample rate']

	while quadraturephase < 0:
		quadraturephase += this['design	sample rate']
	while quadraturephase >= this['design sample rate']:
		quadraturephase -= this['design sample rate']

	this['Sine'] = this['WaveTable'][int(inphase // this['NormalizationFactor'])]
	this['Cosine'] = this['WaveTable'][int(quadraturephase // this['NormalizationFactor'])]
	return this

def Test(state):
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


	NCO = []

	NCO.append({})
	NCO.append({})
	NCO[1] = GetNCOConfig(config, 1, "NCO ")
	NCO[1] = InitNCO(NCO[1])
	print(NCO[1])

	scipy.io.wavfile.write(dirname+"Wavetable.wav", NCO[1]['design sample rate'], NCO[1]['WaveTable'].astype(np.int16))


	# generate 3 seconds of set frequency
	samples = NCO[1]['design sample rate'] * 3
	AudioSine = np.zeros(samples)
	AudioCosine = np.zeros(samples)
	AudioDither = np.zeros(samples)
	for i in range(samples):
		NCO[1] = UpdateNCO(NCO[1])
		AudioSine[i] = NCO[1]['Sine']
		AudioCosine[i] = NCO[1]['Cosine']
		AudioDither[i] = NCO[1]['Dither']

	scipy.io.wavfile.write(dirname+"AudioSine.wav", NCO[1]['design sample rate'], AudioSine.astype(np.int16))
	scipy.io.wavfile.write(dirname+"AudioCosine.wav", NCO[1]['design sample rate'], AudioCosine.astype(np.int16))
	scipy.io.wavfile.write(dirname+"AudioDither.wav", NCO[1]['design sample rate'], AudioDither.astype(np.int16))

	return
