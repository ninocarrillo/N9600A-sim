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

def GetFSK4DemodulatorConfig(config, num):
	this = {}
	
	id_string = "FSK4 Demodulator "
	
	key_string = "enabled"
	try:
		this[f'{key_string}'] = config[f'{id_string}{num}'].getboolean(f'{key_string}')
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)
		
	key_string = "invert"
	try:
		this[f'{key_string}'] = config[f'{id_string}{num}'].getboolean(f'{key_string}')
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)
		
	key_string = "samples per symbol"
	try:
		this[f'{key_string}'] = int(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)
		
	key_string = "threshold"
	try:
		this[f'{key_string}'] = int(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)
		
	key_string = "symbol map"
	try:
		this[f'{key_string}'] = strings.StringToIntArray(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)
		
	return this

def InitFSK4Demod(this):
	return this

def FullProcess(state):
	argv = state['argv']
	config = state['config']

	print(f'Started RRC 4FSK Demodulation process')

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
	FilterDecimator = input_filter.GetInputFilterConfig5(state)
	PulseFilter = pulse_filter.GetRRCFilterConfig(state)
	PulseFilter = pulse_filter.InitRRCFilter(PulseFilter)
	FilterDecimator['Filter'] = np.rint(PulseFilter['Taps'] * PulseFilter['amplitude'])
	FilterDecimator = demod.InitFilterDecimator(FilterDecimator)


	
	DemodulatorCount = 0
	FSK4Demodulator = [{}]
	for DemodulatorNumber in range(4):
		
		if config.has_section(f'FSK4 Demodulator {DemodulatorNumber}'):
			print(f'Reading settings for FSK4 Demodulator {DemodulatorNumber}')
			DemodulatorCount += 1
			FSK4Demodulator.append({})
			FSK4Demodulator[DemodulatorNumber] = GetFSK4DemodulatorConfig(config, DemodulatorNumber)
			FSK4Demodulator[DemodulatorNumber]['InputSampleRate'] = FilterDecimator['OutputSampleRate']
			FSK4Demodulator[DemodulatorNumber] = InitFSK4Demod(FSK4Demodulator[DemodulatorNumber])


	AX25Decoder = [{}]
	Descrambler = [{}]
	for index in range(1,DemodulatorCount+1):
		print(f'Initializing AX25Decoder {index}, Descrambler {index}')
		AX25Decoder.append({})
		Descrambler.append({})
		AX25Decoder[index] = demod.InitAX25Decoder()
		Descrambler[index]['Polynomial'] = int('0x21001',16) # G3RUH poly
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


	plt.figure()
	plt.subplot(231)
	plt.plot(FilterDecimator['FilterBuffer'])
	plt.title('Input Signal')
	FilterDecimator = demod.FilterDecimate(FilterDecimator)
	plt.subplot(232)
	plt.plot(FilterDecimator['FilterBuffer'])
	plt.title('Filtered Signal')
	plt.subplot(233)
	plt.plot(FilterDecimator['Filter'])
	plt.title('Filter Kernel')
	#plt.plot(QPSKDemodulator[1]['SamplePulse'])
	plt.subplot(234)
	plt.plot(FilterDecimator['EnvelopeBuffer'])
	plt.subplot(235)
	#correlator_buffer = [-3,0,0,0,0,0,-3,0,0,0,0,0,1,0,0,0,0,0,3,0,0,0,0,0,3,0,0,0,0,0,3,0,0,0,0,0,-3,0,0,0,0,0,-1,0,0,0,0,0,3,0,0,0,0,0,1,0,0,0,0,0,-1,0,0,0,0,0,1]
	correlator_buffer = [+3,0,0,0,0,0,+3,0,0,0,0,0,+3,0,0,0,0,0,+3,0,0,0,0,0,-3,0,0,0,0,0,-3,0,0,0,0,0,-3,0,0,0,0,0,+3,0,0,0,0,0,-3,0,0,0,0,0,+3,0,0,0,0,0,-3,0,0,0,0,0,+3,0,0,0,0,0,+3,0,0,0,0,0,+3,0,0,0,0,0,+3,0,0,0,0,0,-3,0,0,0,0,0,-3,0,0,0,0,0,+3,0,0,0,0,0,-3,0,0,0,0,0,-3,0,0,0,0,0,+3,0,0,0,0,0,-3,0,0,0,0,0,-3,0,0,0,0,0,-3]
	correlator_buffer = np.flip(np.convolve(correlator_buffer, FilterDecimator['Filter'], 'same'))
	#correlator_buffer = np.flip(np.convolve(correlator_buffer, FilterDecimator['Filter'], 'same'))
	correlation = np.rint(np.convolve(correlator_buffer, FilterDecimator['FilterBuffer'], 'valid'))
	#plt.plot(correlation)
	plt.plot(correlator_buffer)
	plt.subplot(236)
	plt.plot(correlation)
	plt.show()


	scipy.io.wavfile.write(dirname+"FilteredSignal.wav", FilterDecimator['OutputSampleRate'], FilterDecimator['FilterBuffer'].astype(np.int16))

	# Generate and save report file
	report_file_name = f'run{run_number}_report.txt'
	try:
		report_file = open(dirname + report_file_name, 'w+')
	except:
		print('Unable to create report file.')
	with report_file:
		report_file.write('# Command line: ')
		for argument in sys.argv:
			report_file.write(f'{argument} ')
		report_file.write('\n#\n########## Begin Transcribed .ini file: ##########\n')
		try:
			ini_file = open(sys.argv[1])
		except:
			report_file.write('Unable to open .ini file.')
		with ini_file:
			for character in ini_file:
				report_file.write(character)

		report_file.write('\n\n########## End Transcribed .ini file: ##########\n')



		report_file.write('\n')
		report_file.write(fo.GenInt16ArrayC(f'AGCScaleTable', FilterDecimator['AGCScaleTable'], 16))
		report_file.write('\n\n')

		report_file.close()

	return
