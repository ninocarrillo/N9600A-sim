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
		this[f'{key_string}'] = int(float(config[f'{id_string}{num}'][f'{key_string}']))
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "threshold"
	try:
		this[f'{key_string}'] = int(float(config[f'{id_string}{num}'][f'{key_string}']))
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "symbol map"
	try:
		this[f'{key_string}'] = strings.StringToIntArray(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "sample phase"
	try:
		this[f'{key_string}'] = int(float(config[f'{id_string}{num}'][f'{key_string}']))
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "sync rate 1"
	try:
		this[f'{key_string}'] = int(float(config[f'{id_string}{num}'][f'{key_string}']))
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "sync rate 2"
	try:
		this[f'{key_string}'] = int(float(config[f'{id_string}{num}'][f'{key_string}']))
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "sync step"
	try:
		this[f'{key_string}'] = int(float(config[f'{id_string}{num}'][f'{key_string}']))
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "sync period"
	try:
		this[f'{key_string}'] = int(float(config[f'{id_string}{num}'][f'{key_string}']))
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "sync offset 1"
	try:
		this[f'{key_string}'] = int(float(config[f'{id_string}{num}'][f'{key_string}']))
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "sync offset 2"
	try:
		this[f'{key_string}'] = int(float(config[f'{id_string}{num}'][f'{key_string}']))
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "sync delay 1"
	try:
		this[f'{key_string}'] = int(float(config[f'{id_string}{num}'][f'{key_string}']))
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "sync delay 2"
	try:
		this[f'{key_string}'] = int(float(config[f'{id_string}{num}'][f'{key_string}']))
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	return this

def InitFSK4Demod(this):
	this['FirstBitDemap'] = np.zeros(4)
	this['SecondBitDemap'] = np.zeros(4)
	this['Perfect1'] = int(this['threshold'] / 2)
	this['Perfect3'] = int(this['Perfect1'] * 3);
	# index position 0 = symbol value -3, index position 3 = symbol value 3
	for index in range(4):
		value_number = 0
		if this['symbol map'][index] == -3:
			value_number = 0
		if this['symbol map'][index] == -1:
			value_number = 1
		if this['symbol map'][index] == 1:
			value_number = 2
		if this['symbol map'][index] == 3:
			value_number = 3
		this['FirstBitDemap'][value_number] = np.right_shift(index & 2, 1)
		this['SecondBitDemap'][value_number] = index & 1
	this['SyncClock'] = 0
	this['LastSyncSample'] = 0
	this['SampleHistory'] = np.zeros(this['samples per symbol'])
	this['SampleIndex'] = 0
	this['DelayIndex1'] = 1
	this['DelayIndex2'] = 1
	#print(this)
	return this

def Demodulate4FSK(this):
	this['Result'] = np.zeros(len(this['InputAudio'] * 2.2 / this['samples per symbol']))
	this['SampleAudio'] = np.zeros(len(this['InputAudio']))
	this['MissDistance'] = np.zeros(len(this['InputAudio']))
	bit_index = 0
	phase = 0
	symbol = 0
	sample_index = 0
	absolute_sample_index = 0
	for sample in this['InputAudio']:

		# Determine the miss distance
		miss1 = sample - this['Perfect3']
		miss2 = sample - this['Perfect1']
		miss3 = sample + this['Perfect1']
		miss4 = sample + this['Perfect3']

		if abs(miss1) > abs(miss2):
			miss1 = miss2
		if abs(miss4) > abs(miss3):
			miss4 = miss3
		#if abs(miss1) > abs(miss4):
		#	miss1 = miss4
		this['MissDistance'][sample_index] = miss1

		this['SyncClock'] += this['sync step']
		if this['SyncClock'] > ((this['sync period'] // 2) - 1):
			this['SyncClock'] -= this['sync period']
		#if phase == this['sample phase']:
			this['SampleAudio'][sample_index] = sample
			if sample > 0:
				if sample >= this['threshold']:
					symbol = 3
				else:
					symbol = 2
			else:
				if sample <= (-this['threshold']):
					symbol = 0
				else:
					symbol = 1
			# demap the symbol
			this['Result'][bit_index] = this['FirstBitDemap'][symbol]
			bit_index += 1
			this['Result'][bit_index] = this['SecondBitDemap'][symbol]
			bit_index += 1

		phase += 1
		if phase >= this['samples per symbol']:
			phase = 0


		if 1 == 0:
			if this['LastSyncSample'] > 0:
				if miss1 <= 0:
					# Zero Crossing
					this['SyncClock'] *= this['sync rate 1']
					this['SyncClock'] //= this['sync step']
			else:
				if miss1 > 0:
					# Zero Crossing
					this['SyncClock'] *= this['sync rate 1']
					this['SyncClock'] //= this['sync step']
			this['LastSyncSample'] = miss1

		if 2 == 0:
			if this['LastSyncSample'] > 0:
				if sample <= 0:
					# Zero Crossing
					this['SyncClock'] -= this['sync offset 1']
					this['SyncClock'] *= this['sync rate 1']
					this['SyncClock'] //= this['sync step']
					this['SyncClock'] += this['sync offset 1']
			else:
				if sample > 0:
					# Zero Crossing
					this['SyncClock'] -= this['sync offset 1']
					this['SyncClock'] *= this['sync rate 1']
					this['SyncClock'] //= this['sync step']
					this['SyncClock'] += this['sync offset 1']
			this['LastSyncSample'] = sample

		if 3 == 3:
			if this['SampleHistory'][this['DelayIndex1']] < -this['threshold']:
				if sample >= this['threshold']:
					# Zero Crossing
					this['SyncClock'] -= this['sync offset 1']
					this['SyncClock'] *= this['sync rate 1']
					this['SyncClock'] //= this['sync step']
					this['SyncClock'] += this['sync offset 1']


					this['MissDistance'][sample_index] = this['threshold']

			elif this['SampleHistory'][this['DelayIndex1']] >= this['threshold']:
				if sample < -this['threshold']:
					# Zero Crossing
					this['SyncClock'] -= this['sync offset 1']
					this['SyncClock'] *= this['sync rate 1']
					this['SyncClock'] //= this['sync step']
					this['SyncClock'] += this['sync offset 1']

					this['MissDistance'][sample_index] = this['threshold']

		if 4 == 4:
			if this['SampleHistory'][this['DelayIndex2']] < 0:
				if sample >= 0 and sample < this['threshold']:
					# Zero Crossing
					this['SyncClock'] -= this['sync offset 2']
					this['SyncClock'] *= this['sync rate 2']
					this['SyncClock'] //= this['sync step']
					this['SyncClock'] += this['sync offset 2']


					this['MissDistance'][sample_index] = this['threshold']

			else:
				if sample < 0 and sample >= -this['threshold']:
					# Zero Crossing
					this['SyncClock'] -= this['sync offset 2']
					this['SyncClock'] *= this['sync rate 2']
					this['SyncClock'] //= this['sync step']
					this['SyncClock'] += this['sync offset 2']

					this['MissDistance'][sample_index] = this['threshold']

		sample_index += 1
		this['SampleHistory'][this['SampleIndex']] = sample
		this['SampleIndex'] += 1
		if this['SampleIndex'] >= this['samples per symbol']:
			this['SampleIndex'] = 0
		this['DelayIndex1'] = this['SampleIndex'] + this['sync delay 1']
		if this['DelayIndex1'] >= this['samples per symbol']:
			this['DelayIndex1'] -= this['samples per symbol']
		this['DelayIndex2'] = this['SampleIndex'] + this['sync delay 2']
		if this['DelayIndex2'] >= this['samples per symbol']:
			this['DelayIndex2'] -= this['samples per symbol']
		absolute_sample_index += 1
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
		if 1 == 0:
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
	FilterDecimator = demod.FilterDecimate(FilterDecimator)

	FSK4Demodulator[1]['InputAudio'] = FilterDecimator['FilterBuffer']
	if FSK4Demodulator[1]['invert'] == True:
		FSK4Demodulator[1]['InputAudio'] = -FSK4Demodulator[1]['InputAudio']
	FSK4Demodulator[1] = Demodulate4FSK(FSK4Demodulator[1])

	total_packets = 0
	absolute_bit_index = 0
	last_absolute_bit_index = 0
	for data_bit in FSK4Demodulator[1]['Result']:
		Descrambler[1]['NewBit'] = data_bit
		Descrambler[1] = demod.ProgUnscramble(Descrambler[1])
		AX25Decoder[1]['NewBit'] = Descrambler[1]['Result']
		AX25Decoder[1] = demod.ProgDecodeAX25(AX25Decoder[1])
		if AX25Decoder[1]['OutputTrigger'] == True:
			AX25Decoder[1]['OutputTrigger'] = False
			total_packets += 1
			CRC = AX25Decoder[1]['CRC'][0]
			decodernum = '1'
			filename = f'Packet-{total_packets}_CRC-{format(CRC,"#06x")}_decoder-{decodernum}_Index-{absolute_bit_index}'
			bit_difference = absolute_bit_index - last_absolute_bit_index
			bit_delta = f'-{bit_difference}'
			print(f'{dirname+filename+bit_delta}')
			last_absolute_bit_index = absolute_bit_index
			#for byte in AX25Decoder[1]['Output'][16:]:
			#	print(str(byte & 0x7F,"ascii"), end='')
			print(str(AX25Decoder[1]['Output'][16:] & 0x7F, "ascii"), end='\n\n')
			#try:
			#	bin_file = open(dirname + filename + '.bin', '+wb')
			#except:
			#	pass
			#with bin_file:
			#	for byte in AX25Decoder[1]['Output']:
			#		bin_file.write(byte.astype('uint8'))
			#	bin_file.close()
		absolute_bit_index += 1

	if 1 == 0:
		plt.figure()
		plt.subplot(221)
		plt.plot(audio)
		plt.title('Input Signal')
		plt.subplot(222)
		plt.plot(FSK4Demodulator[1]['InputAudio'])
		plt.plot(FSK4Demodulator[1]['SampleAudio'])
		#plt.plot(FSK4Demodulator[1]['MissDistance'])
		plt.title('Filtered Signal')
		plt.subplot(223)
		plt.plot(FilterDecimator['Filter'])
		plt.title('Filter Kernel')
		#plt.plot(QPSKDemodulator[1]['SamplePulse'])
		plt.subplot(224)
		plt.plot(FilterDecimator['EnvelopeBuffer'])
		plt.title('Envelope Buffer')
		plt.show()

	try:
		scipy.io.wavfile.write(dirname+"FilteredSignal.wav", FilterDecimator['OutputSampleRate'], FilterDecimator['FilterBuffer'].astype(np.int16))
	except:
		pass

	if 1 == 0:
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

	return total_packets
