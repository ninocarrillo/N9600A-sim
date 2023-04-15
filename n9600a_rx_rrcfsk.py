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
import n9600a_pulse_filter as pulse_filter
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq

def GetRRCDemodulatorConfig(config, num):
	this = {}
	id_string = "Demodulator "
	key_string = "enabled"
	try:
		this[f'{key_string}'] = config[f'{id_string}{num}'].getboolean(f'{key_string}')
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)
	return this

def GetSlicerConfig(config, num):
	this = {}
	id_string = "Data Slicer"

	key_string = "symbol rate"
	try:
		this[key_string] = int(config[f'{id_string} {num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string} {num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "lock rate"
	try:
		this[key_string] = float(config[f'{id_string} {num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string} {num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "sample rate"
	try:
		this[key_string] = int(config[f'{id_string} {num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string} {num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	id_string = "Symbol Map"
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

def FullProcess(state):
	argv = state['argv']
	config = state['config']

	print(f'Started RRC FSK process')
	print(f'Reading settings for Pulse Filter')
	PulseFilter = pulse_filter.GetRRCFilterConfig(state)
	PulseFilter = pulse_filter.InitRRCFilter(PulseFilter)
	
	# Read only AGC settings into FilterDecimator
	FilterDecimator = pulse_filter.GetAGCConfig(state)
	# Transfer the calculated RRC taps to filter
	FilterDecimator['Filter'] = np.rint(PulseFilter['Taps'] * 32768)
	FilterDecimator['InputSampleRate'] = PulseFilter['sample rate']
	FilterDecimator['InputAGCEnabled'] = True
	FilterDecimator['DecimationRate'] = FilterDecimator['decimation']
	FilterDecimator = pulse_filter.InitFilterDecimator(FilterDecimator)

	Demodulator = []
	DataSlicer = []
	DemodulatorCount = 0
	for DemodulatorNumber in range(4):
		Demodulator.append({})
		DataSlicer.append({})
		if config.has_section(f'Demodulator {DemodulatorNumber}'):
			print(f'Reading settings for Demodulator {DemodulatorNumber}')
			DemodulatorCount += 1
			Demodulator[DemodulatorNumber] = GetRRCDemodulatorConfig(config, DemodulatorNumber)
			DataSlicer[DemodulatorNumber] = GetSlicerConfig(config, DemodulatorNumber)
			DataSlicer[DemodulatorNumber] = demod.InitDataSlicerN(DataSlicer[DemodulatorNumber])
			print(DataSlicer[DemodulatorNumber])

			#DataSlicer[DemodulatorNumber] = demod.InitDataSlicer(DataSlicer[DemodulatorNumber])

	try:
		samplerate, audio = scipy.io.wavfile.read(sys.argv[2])
		print(max(audio))
		# Take two bits of resolution away
		audio = audio >> (16 - PulseFilter['bit count'])
		print(max(audio))
	except:
		print(f'Unable to open wave file {sys.argv[2]}.')
		sys.exit(-2)

	print("Opened file. \r\nSample rate:", samplerate, "\r\nLength:", len(audio))

	# filter input data
	print(f'\nFiltering and decimating audio. ')
	print(FilterDecimator)
	FilterDecimator['FilterBuffer'] = audio
	FilterDecimator = pulse_filter.FilterDecimate2(FilterDecimator)
	print(f'Done.')


	#filtered_audio = np.convolve(audio, PulseFilter['Taps'], 'valid')
	filtered_audio = FilterDecimator['FilterBuffer']
	plt.figure()
	plt.suptitle(f"RRC 4FSK Rolloff Rate:{PulseFilter['rolloff rate']}, Span:{PulseFilter['symbol span']}, Sample Rate:{PulseFilter['sample rate']}")
	plt.subplot(221)
	plt.plot(PulseFilter['Time'], PulseFilter['Taps'], 'b')
	plt.plot(PulseFilter['Time'], PulseFilter['RC'], 'r')
	plt.xticks(PulseFilter['SymbolTicks'])
	plt.xticks(color='w')
	#plt.xlabel("Symbol Intervals")
	plt.title("Impulse Response")
	plt.legend(["RRC", "RC"])
	plt.grid(True)

	plt.subplot(222)
	plt.plot(filtered_audio)
	plt.plot(FilterDecimator['EnvelopeBuffer'])
	plt.plot(FilterDecimator['ThreshBuffer'])
	plt.plot(-FilterDecimator['ThreshBuffer'])
	plt.plot(-FilterDecimator['EnvelopeBuffer'])
	plt.grid(True)
	plt.show()

	# Slice symbols
	for demod_index in range(1, DemodulatorCount + 1):
		for sample in filtered_audio:
			DataSlicer[demod_index]['NewSample'] = sample
			DataSlicer[demod_index] = demod.ProgSliceDataN(DataSlicer[demod_index])
			for data_bit in DataSlicer[demod_index]['Result']:
				print(data_bit)

	return
