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

def GetDPSK3Config(config, num, id_string):
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

def InitDPSK3(this):
	return this

def DemodulateDPSK3(this):
	return this


def FullProcess(state):
	argv = state['argv']
	config = state['config']

	print(f'Started DPSK3 process')

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

	print(f'Reading settings for Filter Decimator')
	FilterDecimator = input_filter.GetInputFilterConfig(state)
	FilterDecimator = demod.InitFilterDecimator(FilterDecimator)

	print(f'Reading settings for DPSK Demodulators')
	DPSKDemodulator = []
	DemodulatorCount = 0
	id_string = "DPSK3 Demodulator "
	for DemodulatorNumber in range(4):
		DPSKDemodulator.append({})
		if config.has_section(f'{id_string}{DemodulatorNumber}'):
			print(f'Reading settings for {id_string}{DemodulatorNumber}')
			DemodulatorCount += 1
			DPSKDemodulator[DemodulatorNumber] = GetDPSK3Config(config, DemodulatorNumber, id_string)
			DPSKDemodulator[DemodulatorNumber]['InputSampleRate'] = FilterDecimator['OutputSampleRate']
			DPSKDemodulator[DemodulatorNumber] = InitDPSK3(DPSKDemodulator[DemodulatorNumber])


	AX25Decoder = [{}]
	Descrambler = [{}]
	for index in range(1,DemodulatorCount+1):
		print(f'Initializing AX25Decoder {index} and Descrambler {index}')
		AX25Decoder.append({})
		Descrambler.append({})
		AX25Decoder[index] = demod.InitAX25Decoder()
		Descrambler[index]['Polynomial'] = int('b101',16) # double differential decoding
		Descrambler[index] = demod.InitDescrambler(Descrambler[index])

	try:
		samplerate, audio = scipy.io.wavfile.read(argv[2])
		# Take two bits of resolution away
		audio = audio >> (16 - FilterDecimator['InputBitCount'])
	except:
		print('Unable to open wave file.')
		sys.exit(-2)

	print("Opened file. \r\nSample rate:", samplerate, "\r\nLength:", len(audio))

	print(f'\nFiltering and decimating audio. ')
	FilterDecimator['FilterBuffer'] = audio
	FilterDecimator = demod.FilterDecimate(FilterDecimator)
	print(f'Done.')

	scipy.io.wavfile.write(dirname+"FilteredSignal.wav", FilterDecimator['OutputSampleRate'], FilterDecimator['FilterBuffer'].astype(np.int16))

	print(f'\nDemodulating audio. ')
	DPSKDemodulator[1]['InputBuffer'] = FilterDecimator['FilterBuffer']
	DPSKDemodulator[1] = DemodulateDPSK3(DPSKDemodulator[1])
	print(f'Done.')

	return
