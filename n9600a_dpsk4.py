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
import time

def GetDPSK4Config(config, num, id_string):
	this = {}
	this['NCO'] = {}
	this['LoopFilter'] = {}
	this['I_LPF'] = {}
	this['Q_LPF'] = {}
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

	key_string = "i lpf iir order"
	try:
		this['I_LPF']['iir order'] = int(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "i lpf iir scale bits"
	try:
		this['I_LPF']['iir scale bits'] = int(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "i lpf iir saturation bits"
	try:
		this['I_LPF']['iir saturation bits'] = int(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "i lpf iir b coefs"
	try:
		this['I_LPF']['iir b coefs'] = strings.StringToIntArray(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "i lpf iir a coefs"
	try:
		this['I_LPF']['iir a coefs'] = strings.StringToIntArray(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "q lpf iir order"
	try:
		this['Q_LPF']['iir order'] = int(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "q lpf iir scale bits"
	try:
		this['Q_LPF']['iir scale bits'] = int(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "q lpf iir saturation bits"
	try:
		this['Q_LPF']['iir saturation bits'] = int(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "q lpf iir b coefs"
	try:
		this['Q_LPF']['iir b coefs'] = strings.StringToIntArray(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "q lpf iir a coefs"
	try:
		this['Q_LPF']['iir a coefs'] = strings.StringToIntArray(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "loop filter iir order"
	try:
		this['LoopFilter']['iir order'] = int(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "loop filter iir scale bits"
	try:
		this['LoopFilter']['iir scale bits'] = int(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "loop filter iir saturation bits"
	try:
		this['LoopFilter']['iir saturation bits'] = int(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "loop filter iir b coefs"
	try:
		this['LoopFilter']['iir b coefs'] = strings.StringToIntArray(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "loop filter iir a coefs"
	try:
		this['LoopFilter']['iir a coefs'] = strings.StringToIntArray(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	return this

def InitDPSK4(this):
	this['NCO'] = nco.InitNCO(this['NCO'])
	this['LoopFilter'] = filters.InitIIR(this['LoopFilter'])
	this['I_LPF'] = filters.InitIIR(this['I_LPF'])
	this['Q_LPF'] = filters.InitIIR(this['Q_LPF'])
	this['BitAccumulator'] = 0
	return this

def DemodulateDPSK4(this):
	this['I_Mixer'] = np.zeros(len(this['InputBuffer']))
	this['Q_Mixer'] = np.zeros(len(this['InputBuffer']))
	this['LoopMixer'] = np.zeros(len(this['InputBuffer']))
	this['I_LPFOutput'] = np.zeros(len(this['InputBuffer']))
	this['Q_LPFOutput'] = np.zeros(len(this['InputBuffer']))
	this['LoopFilterOutput'] = np.zeros(len(this['InputBuffer']))
	this['PhaseAccumulator'] = np.zeros(len(this['InputBuffer']))
	this['SamplePulse'] = np.zeros(len(this['InputBuffer']))
	this['SineOutput'] = np.zeros(len(this['InputBuffer']))
	this['CosineOutput'] = np.zeros(len(this['InputBuffer']))
	this['Result'] = np.zeros(int(len(this['InputBuffer']) * 1.1 * this['NCO']['nco set frequency'] / this['NCO']['nco design sample rate']))
	index = 0
	last_progress_print = -1
	progress_steps = 10000
	sample_count = len(this['InputBuffer']) // progress_steps
	databit_index = 0
	start_time = time.time()
	if this['enabled'] == True:
		for sample in this['InputBuffer']:
			progress = index // sample_count
			if progress != last_progress_print:
				interval_time = time.time() - start_time
				start_time = time.time()
				time_remaining = np.rint(interval_time * (progress_steps - progress))
				last_progress_print = progress
				print(f'Interval: {progress}, Time Remaining: {time_remaining}')
				print(this['NCO']['Control'])
			# print(f'{index}')
			this['NCO'] = nco.UpdateNCO(this['NCO'])
			this['SineOutput'][index] = this['NCO']['Sine']
			this['CosineOutput'][index] = this['NCO']['Cosine']
			# mix sample stream with NCO sinewave to create I branch
			this['I_Mixer'][index] = np.rint(sample * (2*this['NCO']['Sine']) / 65536)
			# low-pass filter the mix product
			this['I_LPF'] = filters.UpdateIIR(this['I_LPF'], this['I_Mixer'][index])
			# print(this)
			this['I_LPFOutput'][index] = this['I_LPF']['Output']

			# mix sample stream with NCO cosine to create Q branch
			this['Q_Mixer'][index] = np.rint(sample * (2*this['NCO']['Cosine']) / 65536)
			# low-pass filter the mix product
			this['Q_LPF'] = filters.UpdateIIR(this['Q_LPF'], this['Q_Mixer'][index])
			this['Q_LPFOutput'][index] = this['Q_LPF']['Output']

			# mix the I and Q branch
			this['LoopMixer'][index] = np.rint(this['Q_LPF']['Output'] * this['I_LPF']['Output'] / 140000) # proportional feedback
			# this['LoopMixer'][index] = np.rint(this['Q_LPF']['Output'] * this['I_LPF']['Output'] / 280000) # proportional feedback
			# this['LoopMixer'][index] = np.rint(this['Q_LPF']['Output'] * this['I_LPF']['Output'] / 34000)
			this['LoopFilter'] = filters.UpdateIIR(this['LoopFilter'], this['LoopMixer'][index])
			this['LoopFilterOutput'][index] = np.rint(this['LoopFilter']['Output'])
			# scale the NCO control signal
			this['NCO']['Control'] = this['LoopFilterOutput'][index]
			# if this['NCO']['Control'] > 250:
				# this['NCO']['Control'] = 250
			# elif this['NCO']['Control'] < -250:
				# this['NCO']['Control'] = -250
			# print(this['NCO']['Control'])


			# save the data signal
			this['I_LPFOutput'][index] = this['I_LPF']['Output']


			# Downsample and threshold data output and save as Result
			if this['NCO']['QuadraturePhaseRollover'] == True:
				this['NCO']['QuadraturePhaseRollover'] = False
				if this['I_LPF']['Output'] > 0:
					try:
						this['Result'][databit_index] =	 1
					except:
						pass
					this['SamplePulse'][index] = 3000
				else:
					try:
						this['Result'][databit_index] = 0
					except:
						pass
					this['SamplePulse'][index] = -3000
				databit_index += 1
			else:
				this['SamplePulse'][index] = 0

			# this['BitAccumulator'] += this['I_LPF']['Output']
			#
			# if this['NCO']['InPhaseRollover'] == True:
			#	this['NCO']['InPhaseRollover'] = False
			#	if this['BitAccumulator'] > 0:
			#		try:
			#			this['Result'][databit_index] =	 1
			#		except:
			#			pass
			#		this['SamplePulse'][index] = this['BitAccumulator']
			#	else:
			#		try:
			#			this['Result'][databit_index] = 0
			#		except:
			#			pass
			#		this['SamplePulse'][index] = this['BitAccumulator']
			#	databit_index += 1
			#	this['BitAccumulator'] = 0
			# else:
			#	this['SamplePulse'][index] = 0


			this['PhaseAccumulator'][index] = this['NCO']['InPhase']



			index += 1
	return this


def FullProcess(state):
	argv = state['argv']
	config = state['config']

	print(f'Started DPSK4 process')

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
	id_string = "DPSK4 Demodulator "
	for DemodulatorNumber in range(4):
		DPSKDemodulator.append({})
		if config.has_section(f'{id_string}{DemodulatorNumber}'):
			print(f'Reading settings for {id_string}{DemodulatorNumber}')
			DemodulatorCount += 1
			DPSKDemodulator[DemodulatorNumber] = GetDPSK4Config(config, DemodulatorNumber, id_string)
			DPSKDemodulator[DemodulatorNumber]['InputSampleRate'] = FilterDecimator['OutputSampleRate']
			DPSKDemodulator[DemodulatorNumber] = InitDPSK4(DPSKDemodulator[DemodulatorNumber])


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
	DPSKDemodulator[1] = DemodulateDPSK4(DPSKDemodulator[1])
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
			try:
				bin_file = open(dirname + filename + '.bin', '+wb')
			except:
				pass
			with bin_file:
				for byte in AX25Decoder[1]['Output']:
					bin_file.write(byte.astype('uint8'))
				bin_file.close()

	scipy.io.wavfile.write(dirname+"DemodSignal.wav", FilterDecimator['OutputSampleRate'], DPSKDemodulator[1]['Result'].astype(np.int16))
	scipy.io.wavfile.write(dirname+"LoopFilter.wav", FilterDecimator['OutputSampleRate'],DPSKDemodulator[1]['LoopFilterOutput'].astype(np.int16))
	scipy.io.wavfile.write(dirname+"I_LPF.wav", FilterDecimator['OutputSampleRate'], DPSKDemodulator[1]['I_LPFOutput'].astype(np.int16))
	scipy.io.wavfile.write(dirname+"Q_LPF.wav", FilterDecimator['OutputSampleRate'], DPSKDemodulator[1]['Q_LPFOutput'].astype(np.int16))

	scipy.io.wavfile.write(dirname+"SamplePulse.wav", FilterDecimator['OutputSampleRate'], DPSKDemodulator[1]['SamplePulse'].astype(np.int16))

	scipy.io.wavfile.write(dirname+"I_Mixer.wav", FilterDecimator['OutputSampleRate'], DPSKDemodulator[1]['I_Mixer'].astype(np.int16))
	scipy.io.wavfile.write(dirname+"Q_Mixer.wav", FilterDecimator['OutputSampleRate'], DPSKDemodulator[1]['Q_Mixer'].astype(np.int16))

	plt.figure()
	plt.subplot(221)
	plt.plot(FilterDecimator['FilterBuffer'])
	plt.plot(DPSKDemodulator[1]['I_LPFOutput'])
	plt.plot(DPSKDemodulator[1]['Q_LPFOutput'])
	plt.title('I and Q LPF Outputs')
	plt.subplot(222)
	plt.plot(DPSKDemodulator[1]['LoopFilterOutput'])
	plt.title('Loop Filter Output')
	plt.subplot(223)
	plt.plot(DPSKDemodulator[1]['I_LPFOutput'])
	plt.plot(DPSKDemodulator[1]['SamplePulse'])
	plt.title('I Output and Sample Pulse')
	plt.subplot(224)
	plt.plot(DPSKDemodulator[1]['SineOutput'])
	plt.plot(DPSKDemodulator[1]['CosineOutput'])
	plt.title('NCO Outputs')
	plt.show()
	return
