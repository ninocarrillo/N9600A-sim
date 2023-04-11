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

			#DataSlicer[DemodulatorNumber] = demod.InitDataSlicer(DataSlicer[DemodulatorNumber])

	try:
		samplerate, audio = scipy.io.wavfile.read(sys.argv[2])
		# Take two bits of resolution away
		audio = audio >> (16 - FilterDecimator['InputBitCount'])
	except:
		print(f'Unable to open wave file {sys.argv[2]}.')
		sys.exit(-2)

	print("Opened file. \r\nSample rate:", samplerate, "\r\nLength:", len(audio))




	return
