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
import n9600a_filters as filters
import matplotlib as mpl
import matplotlib.pyplot as plt

def GetDPSK3Config(config, num, id_string):
	this = {}
	this['NCO'] = {}
	this['LoopFilter'] = {}
	this['OutputFilter'] = {}
	key_string = "enabled"
	try:
		this[f'{key_string}'] = config[f'{id_string}{num}'].getboolean(f'{key_string}')
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "nco design sample rate"
	try:
		this['NCO'][f'{key_string}'] = int(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "nco wavetable size"
	try:
		this['NCO'][f'{key_string}'] = int(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "nco set frequency"
	try:
		this['NCO'][f'{key_string}'] = int(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "nco amplitude bits"
	try:
		this['NCO'][f'{key_string}'] = int(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "nco phase dither bits"
	try:
		this['NCO'][f'{key_string}'] = int(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "loop filter shift"
	try:
		this['LoopFilter']['OutputShift'] = int(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "output filter shift"
	try:
		this['OutputFilter']['OutputShift'] = int(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "loop filter taps"
	try:
		this['LoopFilter']['Taps'] = strings.StringToIntArray(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "output filter taps"
	try:
		this['OutputFilter']['Taps'] = strings.StringToIntArray(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	return this

def InitDPSK3(this):
	this['NCO'] = nco.InitNCO(this['NCO'])
	this['LoopFilter'] = filters.InitFIR(this['LoopFilter'])
	this['OutputFilter'] = filters.InitFIR(this['OutputFilter'])
	return this

def DemodulateDPSK3(this):
	this['FirstMixer'] = np.zeros(len(this['InputBuffer']))
	this['SecondMixer'] = np.zeros(len(this['InputBuffer']))
	this['ThirdMixer'] = np.zeros(len(this['InputBuffer']))
	this['PhaseAccumulator'] = np.zeros(len(this['InputBuffer']))
	this['LoopFilterOutput'] = np.zeros(len(this['InputBuffer']))
	this['DataFilterOutput'] = np.zeros(len(this['InputBuffer']))
	this['SamplePulse'] = np.zeros(len(this['InputBuffer']))
	this['Result'] = np.zeros(int(len(this['InputBuffer']) * 1.1 * this['NCO']['nco set frequency'] / this['NCO']['nco design sample rate']))
	index = 0
	databit_index = 0
	if this['enabled'] == True:
		for sample in this['InputBuffer']:
			this['NCO'] = nco.UpdateNCO(this['NCO'])


			# mix sample stream with NCO negative sin output to create carrier error
			this['FirstMixer'][index] = np.rint(sample * (-this['NCO']['Sine']) / 65536)

			# LPF carrier error
			this['LoopFilter'] = filters.UpdateFIR(this['LoopFilter'], this['FirstMixer'][index])
			this['LoopFilterOutput'][index] = this['LoopFilter']['Output']

			# mix sample stream with NCO cosine to create data output
			this['SecondMixer'][index] = np.rint(sample * this['NCO']['Cosine'] / 65536)

			# LPF data output
			this['OutputFilter'] = filters.UpdateFIR(this['OutputFilter'], this['SecondMixer'][index])
			this['DataFilterOutput'][index] = this['OutputFilter']['Output']

			# mix data output with carrier error to create NCO control signal
			this['ThirdMixer'][index] = np.rint(this['OutputFilter']['Output'] * this['LoopFilter']['Output']) / 4096

			# scale the NCO control signal
			this['NCO']['Control'] = np.rint(this['ThirdMixer'][index])


			# Downsample and threshold data output and save as Result
			if this['NCO']['QuadraturePhaseRollover'] == True:
				this['NCO']['QuadraturePhaseRollover'] = False
				if this['OutputFilter']['Output'] > 0:
					this['Result'][databit_index] =  1
					this['SamplePulse'][index] = 1000
				else:
					this['Result'][databit_index] = 0
					this['SamplePulse'][index] = -1000
				databit_index += 1
			else:
				this['SamplePulse'][index] = 0


			this['PhaseAccumulator'][index] = this['NCO']['InPhase']



			index += 1
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
		Descrambler[index]['Polynomial'] = int('0x5',16) # double differential decoding
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

	print(f'\nDifferential decoding and AX25 decoding data. ')

	total_packets = 0

	for data_bit in DPSKDemodulator[1]['Result']:
		Descrambler[1]['NewBit'] = data_bit
		Descrambler[1] = demod.ProgUnscramble(Descrambler[1])
		AX25Decoder[1]['NewBit'] = Descrambler[1]['Result']
		AX25Decoder[1] = demod.ProgDecodeAX25(AX25Decoder[1])
		if AX25Decoder[1]['OutputTrigger'] == True:
			AX25Decoder[1]['OutputTrigger'] = False
			# Check for uniqueness
			total_packets += 1
			CRC = AX25Decoder[1]['CRC'][0]
			decodernum = '1'
			filename = f'Packet-{total_packets}_CRC-{format(CRC,"#06x")}_decoder-{decodernum}_Index-{index}'
			print(f'{dirname+filename}')
			# try:
			# 	bin_file = open(dirname + filename + '.bin', '+wb')
			# except:
			# 	pass
			# with bin_file:
			# 	for byte in AX25Decoder[2]['Output']:
			# 		bin_file.write(byte.astype('uint8'))
			# 	bin_file.close()

	scipy.io.wavfile.write(dirname+"DemodSignal.wav", FilterDecimator['OutputSampleRate'], DPSKDemodulator[1]['Result'].astype(np.int16))
	scipy.io.wavfile.write(dirname+"LoopFilter.wav", FilterDecimator['OutputSampleRate'],DPSKDemodulator[1]['LoopFilterOutput'].astype(np.int16))
	scipy.io.wavfile.write(dirname+"DataFilter.wav", FilterDecimator['OutputSampleRate'], DPSKDemodulator[1]['DataFilterOutput'].astype(np.int16))
	scipy.io.wavfile.write(dirname+"ThirdMixer.wav", FilterDecimator['OutputSampleRate'], DPSKDemodulator[1]['ThirdMixer'].astype(np.int16))
	scipy.io.wavfile.write(dirname+"PhaseAccumulator.wav", FilterDecimator['OutputSampleRate'], DPSKDemodulator[1]['PhaseAccumulator'].astype(np.int16))

	scipy.io.wavfile.write(dirname+"SamplePulse.wav", FilterDecimator['OutputSampleRate'], DPSKDemodulator[1]['SamplePulse'].astype(np.int16))

	scipy.io.wavfile.write(dirname+"SecondMixer.wav", FilterDecimator['OutputSampleRate'], DPSKDemodulator[1]['SecondMixer'].astype(np.int16))
	scipy.io.wavfile.write(dirname+"FirstMixer.wav", FilterDecimator['OutputSampleRate'], DPSKDemodulator[1]['FirstMixer'].astype(np.int16))

	x = np.arange(len(DPSKDemodulator[1]['DataFilterOutput']))
	fig1,ax1 = plt.subplots()
	ax1.plot(DPSKDemodulator[1]['DataFilterOutput'])
	ax1.plot(DPSKDemodulator[1]['SamplePulse'])
	plt.show()
	fig2,ax2 = plt.subplots()
	ax2.plot(DPSKDemodulator[1]['ThirdMixer'])
	plt.show()
	return
