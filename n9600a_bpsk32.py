import sys
import configparser
import struct
import scipy.io.wavfile
import numpy as np
import os
import format_output as fo
import n9600a_strings as strings
import n9600a_input_filter as input_filter
import n9600a_pulse_filter as pulse_filter
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq
import random
import n9600a_nco as nco
import n9600a_progdemod as demod
import n9600a_filters as filters
import time
from n9600a_fm import ModulateFM
from n9600a_analysis import AnalyzeSpectrum
from scipy.signal import firwin


def GetBPSKDemodConfig(config, num, id_string):
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
		this['NCO'][f'{key_string}'] = int(float(config[f'{id_string}{num}'][f'{key_string}']))

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

	key_string = "loop filter p"
	try:
		this['LoopFilter'][f'{key_string}'] = float(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "loop filter i"
	try:
		this['LoopFilter'][f'{key_string}'] = float(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "loop filter i max"
	try:
		this['LoopFilter'][f'{key_string}'] = int(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "loop integrator trim"
	try:
		this['LoopFilter'][f'{key_string}'] = int(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "intermediate decimation"
	try:
		this[f'{key_string}'] = int(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "post decimation"
	try:
		this[f'{key_string}'] = int(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	return this

def InitBPSKDemod(this):
	this['NCO'] = nco.InitNCO(this['NCO'])
	this['LoopFilter'] = filters.InitIIR(this['LoopFilter'])
	this['LoopFilter2'] = this['LoopFilter'].copy()
	this['LoopFilter3'] = this['LoopFilter'].copy()
	this['LoopFilter2'] = filters.InitIIR(this['LoopFilter2'])
	this['LoopFilter3'] = filters.InitIIR(this['LoopFilter3'])
	this['I_LPF'] = filters.InitIIR(this['I_LPF'])
	this['Q_LPF'] = filters.InitIIR(this['Q_LPF'])
	this['BitAccumulator'] = 0
	this['OutputSampleRate'] = this['InputSampleRate'] // this['intermediate decimation']
	return this

def DemodulateBPSK(this):
	this['I_Mixer'] = np.zeros(len(this['InputBuffer']))
	this['Q_Mixer'] = np.zeros(len(this['InputBuffer']))
	this['LoopMixer'] = np.zeros(len(this['InputBuffer']))
	this['I_LPFOutput'] = np.zeros(len(this['InputBuffer']))
	this['Q_LPFOutput'] = np.zeros(len(this['InputBuffer']))
	this['LoopFilterOutput'] = np.zeros(len(this['InputBuffer']))
	this['PhaseAccumulator'] = np.zeros(len(this['InputBuffer']))
	this['SamplePulse'] = np.zeros(len(this['InputBuffer']))
	this['SineOutput'] = np.zeros(len(this['InputBuffer']))
	this['LoopIntegral'] = np.zeros(len(this['InputBuffer']))

	this['NCOControlOutput'] = np.zeros(len(this['InputBuffer']))
	this['CosineOutput'] = np.zeros(len(this['InputBuffer']))
	this['Result'] = np.zeros(int(len(this['InputBuffer']) * 1.1 * this['NCO']['nco set frequency'] / this['NCO']['nco design sample rate']))
	index = 0
	last_progress_print = -1
	progress_steps = 100
	sample_count = len(this['InputBuffer']) // progress_steps
	databit_index = 0
	start_time = time.time()
	if this['enabled'] == True:
		integral = 0
		for sample in this['InputBuffer']:
			progress = index // sample_count
			if progress != last_progress_print:
				last_progress_print = progress
				interval_time = time.time() - start_time
				if (progress > 0):
					time_remaining = np.rint(interval_time * (progress_steps - progress) / progress)
					print(f'Interval: {progress}, Time Remaining: {time_remaining}')
					print(this['NCO']['Control'])
			this['NCO'] = nco.UpdateNCO(this['NCO'])
			this['SineOutput'][index] = this['NCO']['Sine']
			this['CosineOutput'][index] = this['NCO']['Cosine']
			# mix sample stream with NCO sinewave to create I branch
			this['I_Mixer'][index] = np.rint(sample * (this['NCO']['Sine']) // 32768) # simulate 15 bit fractional integer multiplication
			# low-pass filter the mix product
			this['I_LPF'] = filters.UpdateIIR(this['I_LPF'], this['I_Mixer'][index])
			this['I_LPFOutput'][index] = this['I_LPF']['Output']

			# mix sample stream with NCO cosine to create Q branch
			this['Q_Mixer'][index] = np.rint(sample * (this['NCO']['Cosine']) // 32768) # simulate 15 bit fractional integer multiplication
			# low-pass filter the mix product
			this['Q_LPF'] = filters.UpdateIIR(this['Q_LPF'], this['Q_Mixer'][index])
			this['Q_LPFOutput'][index] = this['Q_LPF']['Output']

			# mix the I and Q branch
			this['LoopMixer'][index] = np.rint((this['Q_LPF']['Output'] * this['I_LPF']['Output']) // 32768) # simulate 15 bit fractional integer multiplication
			this['LoopFilter'] = filters.UpdateIIR(this['LoopFilter'], this['LoopMixer'][index])
			this['LoopFilterOutput'][index] = np.rint(this['LoopFilter']['Output'])
			# scale the NCO control signal
			p = np.rint(this['LoopFilterOutput'][index] * this['LoopFilter']['loop filter p'])
			integral += this['LoopFilterOutput'][index] + this['LoopFilter']['loop integrator trim']
			if abs(integral) > this['LoopFilter']['loop filter i max']:
				integral = 0

			this['LoopIntegral'][index] = integral

			i = np.rint(integral * this['LoopFilter']['loop filter i'])

			this['NCO']['Control'] = p + i

			this['NCOControlOutput'][index] = this['NCO']['Control']

			# save the data signal
			this['I_LPFOutput'][index] = this['I_LPF']['Output']

			this['PhaseAccumulator'][index] = this['NCO']['InPhase']



			index += 1

		this['I_LPFOutput'] = this['I_LPFOutput'][0::this['intermediate decimation']]
		this['Q_LPFOutput'] = this['Q_LPFOutput'][0::this['intermediate decimation']]
	return this



def FullProcess(state):
	argv = state['argv']
	config = state['config']

	if state['reports'] == True:
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

	print(f'Started BPSK Demodulation process')
	print(f'Reading settings for Filter Decimator')

	FilterDecimator = input_filter.GetInputFilterConfig(state)
	FilterDecimator = demod.InitFilterDecimator(FilterDecimator)

	DemodulatorNumber = 1
	id_string = "BPSK Demodulator "
	print(f'Reading settings for {id_string}{DemodulatorNumber}')
	BPSKDemodulator = GetBPSKDemodConfig(config, DemodulatorNumber, id_string)
	BPSKDemodulator['InputSampleRate'] = FilterDecimator['OutputSampleRate']
	BPSKDemodulator = InitBPSKDemod(BPSKDemodulator)

	ReceivePulseFilter = pulse_filter.GetRRCFilterConfig(state)
	ReceivePulseFilter['sample rate'] = BPSKDemodulator['OutputSampleRate']
	ReceivePulseFilter = pulse_filter.InitRRCFilter(ReceivePulseFilter)
	ReceivePulseFilter['SymbolMap'] = pulse_filter.GetSymbolMapConfig(state)

	DataSlicer = {}
	try:
		DataSlicer['BitRate'] = int(config[f'Data Slicer {DemodulatorNumber}']['slicer bit rate'])
	except:
		print(f'{sys.argv[1]} [Data Slicer {DemodulatorNumber}] \'slicer bit rate\' is missing or invalid')
		sys.exit(-2)
	try:
		DataSlicer['Rate'] = float(config[f'Data Slicer {DemodulatorNumber}']['slicer lock rate'])
	except:
		print(f'{sys.argv[1]} [Data Slicer {DemodulatorNumber}] \'slicer lock rate\' is missing or invalid')
		sys.exit(-2)

	DataSlicer['InputSampleRate'] = BPSKDemodulator['OutputSampleRate'] // BPSKDemodulator['post decimation']
	print('Data Slicer Sample Rate: ', DataSlicer['InputSampleRate'])
	DataSlicer = demod.InitDataSlicer(DataSlicer)

	print(f'Initializing Descrambler and AX25Decoder')
	AX25Decoder = demod.InitAX25Decoder()
	Descrambler = {}
	Descrambler['Polynomial'] = int('0x63003',16) # G3RUH poly * differential decoding
	Descrambler= demod.InitDescrambler(Descrambler)

	try:
		samplerate, audio = scipy.io.wavfile.read(argv[2])
		# Adjust bit resolution as specified
		audio = audio >> (16 - FilterDecimator['InputBitCount'])
	except:
		print('Unable to open wave file.')
		sys.exit(-2)

	print("Opened file. \r\nSample rate:", samplerate, "\r\nLength:", len(audio))

	print(f'\nFiltering and decimating audio. ')
	FilterDecimator['FilterBuffer'] = audio
	FilterDecimator = demod.FilterDecimate(FilterDecimator)
	print(f'Done.')

	if state['reports'] == True:
		scipy.io.wavfile.write(dirname+"FilteredSignal.wav", FilterDecimator['OutputSampleRate'], FilterDecimator['FilterBuffer'].astype(np.int16))

	print(f'\nDemodulating audio. ')
	BPSKDemodulator['InputBuffer'] = FilterDecimator['FilterBuffer']
	BPSKDemodulator = DemodulateBPSK(BPSKDemodulator)
	print(f'Done.')

	ReceivePulseFilter['Taps'] = np.rint(ReceivePulseFilter['Taps'] * 8191)

	FilteredIOutput = np.convolve(BPSKDemodulator['I_LPFOutput'], ReceivePulseFilter['Taps'], 'valid') // 65536
	FilteredQOutput = np.convolve(BPSKDemodulator['Q_LPFOutput'], ReceivePulseFilter['Taps'], 'valid') // 65536

	FilteredIOutput = FilteredIOutput[0::BPSKDemodulator['post decimation']]
	FilteredQOutput = FilteredQOutput[0::BPSKDemodulator['post decimation']]

	print(f'\nSlicing, differential decoding, and AX25 decoding data. ')

	total_packets = 0

	loop_count = len(FilteredIOutput)
	for index in range(loop_count):
		DataSlicer['NewSample'] = FilteredIOutput[index]
		DataSlicer = demod.ProgSliceData(DataSlicer)
		for data_bit in DataSlicer['Result']:
			Descrambler['NewBit'] = data_bit
			Descrambler = demod.ProgUnscramble(Descrambler)
			AX25Decoder['NewBit'] = Descrambler['Result']
			AX25Decoder = demod.ProgDecodeAX25(AX25Decoder)
			if AX25Decoder['OutputTrigger'] == True:
				AX25Decoder['OutputTrigger'] = False
				# Check for uniqueness
				total_packets += 1
				CRC = AX25Decoder['CRC'][0]
				filename = f'Packet-{total_packets}_CRC-{format(CRC,"#06x")}_Index-{index}'
				print(f'{filename}')
				if state['reports'] == True:
					try:
						bin_file = open(dirname + filename + '.txt', '+wb')
					except:
						pass
					with bin_file:
						for byte in AX25Decoder['Output']:
							bin_file.write(byte.astype('uint8'))
						bin_file.close()

	print(f'Total packets decoded: {total_packets}')

	if state['reports'] == True:
		scipy.io.wavfile.write(dirname+"DemodSignal.wav", FilterDecimator['OutputSampleRate'], BPSKDemodulator['Result'].astype(np.int16))
		scipy.io.wavfile.write(dirname+"LoopFilter.wav", FilterDecimator['OutputSampleRate'],BPSKDemodulator['LoopFilterOutput'].astype(np.int16))
		scipy.io.wavfile.write(dirname+"I_LPF.wav", FilterDecimator['OutputSampleRate'], BPSKDemodulator['I_LPFOutput'].astype(np.int16))
		scipy.io.wavfile.write(dirname+"Q_LPF.wav", FilterDecimator['OutputSampleRate'], BPSKDemodulator['Q_LPFOutput'].astype(np.int16))

		scipy.io.wavfile.write(dirname+"SamplePulse.wav", FilterDecimator['OutputSampleRate'], BPSKDemodulator['SamplePulse'].astype(np.int16))

		scipy.io.wavfile.write(dirname+"I_Mixer.wav", FilterDecimator['OutputSampleRate'], BPSKDemodulator['I_Mixer'].astype(np.int16))
		scipy.io.wavfile.write(dirname+"Q_Mixer.wav", FilterDecimator['OutputSampleRate'], BPSKDemodulator['Q_Mixer'].astype(np.int16))


	#print(PulseFilter['Taps'])

	if state['plots'] == True:
		plt.figure()
		plt.subplot(221)
		#plt.plot(FilterDecimator['FilterBuffer'])
		plt.plot(FilterDecimator['FilterBuffer'])
		plt.plot(BPSKDemodulator['I_LPFOutput'])
		plt.plot(FilterDecimator['EnvelopeBuffer'])
		#plt.plot(BPSKDemodulator[1]['Q_LPFOutput'])
		plt.title('I LPF Output')
		plt.legend(['Filtered Input','I_LPF Output','Envelope'])
		plt.subplot(222)
		plt.plot(BPSKDemodulator['LoopMixer'])
		plt.plot(BPSKDemodulator['LoopIntegral'])
		plt.plot(BPSKDemodulator['LoopFilterOutput'])
		plt.legend(['LoopMixer','LoopIntegral', 'LoopFilter'])
		plt.title('Loop Filter Output')
		plt.subplot(223)
		#plt.plot(BPSKDemodulator[1]['I_LPFOutput'])
		#plt.plot(BPSKDemodulator[1]['SamplePulse'])
		plt.plot(FilteredIOutput)
		plt.title('Filtered I Output')
		plt.subplot(224)
		plt.plot(BPSKDemodulator['NCOControlOutput'])
		plt.title('NCO Control')
		plt.show()



	if state['reports'] == True:
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

			report_file.write('\n')
			report_file.write(fo.GenInt16ArrayC(f'ReceiveFilter', ReceivePulseFilter['Taps'], ReceivePulseFilter['Oversample']))
			report_file.write('\n\n')

			report_file.write(f'Total packets decoded: {total_packets}')

			report_file.close()

	return total_packets


def ModulateGauss(state):
	argv = state['argv']
	config = state['config']
	print(f'Started Shaped BPSK Modulator process')
	print(f'Reading settings for Gauss Pulse Shaping Filter')
	PulseFilter = pulse_filter.GetGaussFilterConfig(state)
	PulseFilter = pulse_filter.InitGaussFilter(PulseFilter)
	PulseFilter['SymbolMap'] = pulse_filter.GetSymbolMapConfig(state)
	NCO = nco.GetNCOConfig(config, 1, "TX NCO ")
	NCO['Amplitude'] = PulseFilter['amplitude']
	NCO = nco.InitNCO(NCO)

	PulseFilter = pulse_filter.GenPulseFilterPatterns(PulseFilter)
	print(max(PulseFilter['FilterPatterns']))
	plt.plot(PulseFilter['FilterPatterns'])
	plt.show()

	#BitStream = pulse_filter.ExpandSampleStream(state['InputData'], PulseFilter)
	BitStream = pulse_filter.BytesToSymbols(state['InputData'], PulseFilter)
	plt.figure()
	plt.plot(BitStream)
	plt.show()
	#ModulatingWaveform = np.convolve(PulseFilter['Taps'], BitStream)
	#ModAmplitude = 64
	#ModulatingWaveform = np.rint(ModAmplitude * ModulatingWaveform / max(ModulatingWaveform))
	ModulatingWaveform = np.zeros(len(BitStream) * PulseFilter['Oversample'] * PulseFilter['undersample'])

	i = 0
	shift_register = int(0)
	filter_mask = (2**PulseFilter['symbol span']) - 1
	for bit in BitStream:
		shift_register <<= 1
		if bit == 1:
			shift_register |= 1
		shift_register = shift_register & filter_mask
		for phase in range(PulseFilter['Oversample']):
			for subphase in range(PulseFilter['undersample']):
				ModulatingWaveform[i] = PulseFilter['FilterPatterns'][(shift_register * PulseFilter['Oversample']) + phase]
				i += 1

	plt.figure()
	plt.plot(ModulatingWaveform)
	plt.show()

	Baseband = np.zeros(len(ModulatingWaveform))
	i = 0
	for Amplitude in ModulatingWaveform:
		NCO = nco.UpdateNCO(NCO)
		Baseband[i] = Amplitude * NCO['Sine'] // PulseFilter['amplitude']
		i += 1

	plt.figure()
	plt.plot(Baseband)
	plt.show()

	plt.figure()
	plt.subplot(221)
	plt.plot(PulseFilter['Time'], PulseFilter['Taps'], 'b')
	#plt.plot(PulseFilter['Time'], PulseFilter['RC'], 'r')
	plt.xticks(PulseFilter['SymbolTicks'])
	plt.xticks(color='w')
	plt.grid(True)

	plt.subplot(222)
	plt.plot(ModulatingWaveform, 'b')

	eye_data = pulse_filter.GenEyeData2(ModulatingWaveform, PulseFilter['Oversample'], 0)
	plt.subplot(223)
	plt.plot(eye_data)


	fft_n = len(ModulatingWaveform)
	x = np.linspace(0.0, fft_n * PulseFilter['TimeStep'], fft_n, endpoint = False)
	x_fft = fftfreq(fft_n, PulseFilter['TimeStep'])[:fft_n//2]
	ModulatingWaveform_fft = fft(ModulatingWaveform)
	ModulatingWaveform_fft = fft(ModulatingWaveform)
	fft_max = max(abs(ModulatingWaveform_fft))
	ModulatingWaveform_fft = ModulatingWaveform_fft / fft_max
	plt.subplot(224)
	plt.plot(x_fft, 10*np.log(np.abs(ModulatingWaveform_fft[0:fft_n//2])))
	plt.xlim(0,PulseFilter['symbol rate'] * 4)
	plt.ylim(-100,10)
	plt.show()

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

	ModulatingWaveform = ModulatingWaveform / max(ModulatingWaveform)
	ModulatingWaveform = ModulatingWaveform * 32767
	scipy.io.wavfile.write(dirname+"ModSignal.wav", PulseFilter['sample rate'], ModulatingWaveform.astype(np.int16))
	scipy.io.wavfile.write(dirname+"Baseband.wav", PulseFilter['sample rate'], Baseband.astype(np.int16))


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
		report_file.write(fo.GenInt16ArrayC(f'BPSKFilterPatterns', PulseFilter['FilterPatterns'], PulseFilter['Oversample']))
		report_file.write('\n\n')


		report_file.write('\n')
		report_file.write(fo.GenInt16ArrayC(f'SineTable', NCO['WaveTable'], 8))
		report_file.write(fo.GenInt16ArrayC(f'QuarterSineTable', NCO['WaveTable'][0:(len(NCO['WaveTable']) // 4) + 1], 8))
		report_file.write('\n\n')


		report_file.close()


	return



def ModulateRRC(state):
	argv = state['argv']
	config = state['config']
	print(f'Started Shaped BPSK Modulator process')
	print(f'Reading settings for RRC Pulse Shaping Filter')
	#PulseFilter = pulse_filter.GetGaussFilterConfig(state)
	#PulseFilter = pulse_filter.InitGaussFilter(PulseFilter)
	PulseFilter = pulse_filter.GetRRCFilterConfig(state)
	PulseFilter = pulse_filter.InitRRCFilter(PulseFilter)



	#plt.figure()
	#plt.plot(PulseFilter['FilterWindow'])
	#plt.title('Filter Window')
	#plt.show()

	#plt.figure()
	#plt.plot(PulseFilter['Taps'])
	#plt.plot(PulseFilter['WindowedRC'])
	#plt.title('Windowed RRC')
	#plt.show()

	PulseFilter['SymbolMap'] = pulse_filter.GetSymbolMapConfig(state)
	NCO = nco.GetNCOConfig(config, 1, "TX NCO ")
	NCO['Amplitude'] = PulseFilter['amplitude']
	NCO = nco.InitNCO(NCO)


	channel_filter = firwin(
				(PulseFilter['Oversample'] * 8) + 1,
				400,
				pass_zero='highpass',
				fs=PulseFilter['sample rate']
			)


	#PulseFilter['Taps'] = np.convolve(PulseFilter['Taps'], channel_filter, 'same')

	#PulseFilter['Taps'] = firwin(
	#			(PulseFilter['Oversample'] * PulseFilter['symbol span']),
	#			995,
	#			pass_zero='lowpass',
	#			fs=PulseFilter['sample rate']
	#		)
	PulseFilter['Taps'] /= max(PulseFilter['Taps'])

	PulseFilter = pulse_filter.GenPulseFilterPatterns(PulseFilter)

	#BitStream = pulse_filter.ExpandSampleStream(state['InputData'], PulseFilter)
	BitStream = pulse_filter.BytesToSymbols(state['InputData'], PulseFilter)
	#plt.figure()
	#plt.plot(BitStream)
	#plt.title('BitStream')
	#plt.show()
	#ModulatingWaveform = np.convolve(PulseFilter['Taps'], BitStream)
	#ModAmplitude = 64
	#ModulatingWaveform = np.rint(ModAmplitude * ModulatingWaveform / max(ModulatingWaveform))
	ModulatingWaveform = np.zeros(len(BitStream) * PulseFilter['Oversample'] * PulseFilter['undersample'])

	i = 0
	shift_register = int(0)
	filter_mask = (2**PulseFilter['symbol span']) - 1
	for bit in BitStream:
		shift_register <<= 1
		if bit == 1:
			shift_register |= 1
		shift_register = shift_register & filter_mask
		for phase in range(PulseFilter['Oversample']):
			ModulatingWaveform[i] = PulseFilter['FilterPatterns'][(shift_register * PulseFilter['Oversample']) + phase]
			i += 1


	#plt.figure()
	#plt.plot(ModulatingWaveform)
	#plt.title('Modulating Waveform')
	#plt.show()

	Baseband = np.zeros(len(ModulatingWaveform))
	for i in range(len(ModulatingWaveform)):
		NCO = nco.UpdateNCO(NCO)
		Baseband[i] = ModulatingWaveform[i] * NCO['Sine']


	Baseband = np.convolve(Baseband, channel_filter)

	Baseband = Baseband / max(Baseband)

	plt.figure()
	plt.suptitle(f"RRC Rolloff Rate:{PulseFilter['rolloff rate']}, Span:{PulseFilter['symbol span']}, Sample Rate:{PulseFilter['sample rate']}")
	plt.subplot(221)
	plt.plot(PulseFilter['Taps'], 'b')
	#plt.plot(channel_filter, 'r')
	plt.xticks(PulseFilter['SymbolTicks'])
	plt.xticks(color='w')
	plt.grid(True)
	plt.title("RRC Filter Taps")


	DemodulatedWaveform = np.convolve(PulseFilter['Taps'], ModulatingWaveform, 'valid')
	mod_eye_data = pulse_filter.GenEyeData2(ModulatingWaveform / PulseFilter['amplitude'], PulseFilter['Oversample'], 0)
	demod_eye_data = pulse_filter.GenEyeData2(DemodulatedWaveform / PulseFilter['amplitude'] * 16, PulseFilter['Oversample'], 0)
	plt.subplot(223)
	plt.plot(mod_eye_data)
	plt.title("Transmit Eye")
	plt.xlabel("Sample Index")
	plt.subplot(224)
	plt.plot(demod_eye_data)
	plt.title("Receive Eye")
	plt.xlabel("Sample Index")


	plt.subplot(222)
	baseband_psd = AnalyzeSpectrum(Baseband, PulseFilter['sample rate'], 0.99)
	plt.plot(baseband_psd[0], baseband_psd[1], '.', markersize=1)
	plt.plot(baseband_psd[2], baseband_psd[3])
	plt.grid(True)
	plt.xlim(0,baseband_psd[4])
	plt.ylim(-60,10)
	plt.title(f"TX Audio Spectrum, 99% Power Bandwidth: {round(baseband_psd[4] / 2000, 1)} kHz")
	plt.xlabel("Hz")
	plt.ylabel("dBFS")
	plt.show()

	fm_waveform = ModulateFM(Baseband / max(Baseband), PulseFilter['inner deviation'], PulseFilter['sample rate'])
	plt.figure()
	plt.plot(fm_waveform)
	plt.show()

	plt.figure()
	psd_9999 = AnalyzeSpectrum(fm_waveform, PulseFilter['sample rate'], 0.9999)
	psd_999 = AnalyzeSpectrum(fm_waveform, PulseFilter['sample rate'], 0.999)
	psd_99 = AnalyzeSpectrum(fm_waveform, PulseFilter['sample rate'], 0.99)
	# returns [x_fft, waveform_psd, obw_x, obw_mask, obw]
	plt.plot(psd_99[2], psd_99[3],'green')
	plt.plot(psd_999[2], psd_999[3], 'orange')
	plt.plot(psd_9999[2], psd_9999[3], 'gray')
	plt.legend([f'99%: {round(psd_99[4]/1000,1)} kHz', f'99.9%: {round(psd_999[4]/1000,1)} kHz', f'99.99%: {round(psd_9999[4]/1000,1)} kHz'])
	plt.plot(psd_999[0], psd_999[1], '.', markersize=5)
	#plt.plot(psd_999[0], psd_999[1])
	plt.xlim(-6*PulseFilter['symbol rate'],6*PulseFilter['symbol rate'])
	plt.ylim(-100,10)
	plt.ylabel("dBFS")
	plt.xlabel("Deviation from Carrier Frequency, Hz")
	plt.grid(True)
	plt.title(f"Power Spectrum, {len(state['InputData'])} Random Bytes\nSymbol Rate: {PulseFilter['symbol rate']}, Inner Deviation: {PulseFilter['inner deviation']}")
	plt.show()

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

	ModulatingWaveform = ModulatingWaveform / max(ModulatingWaveform)
	ModulatingWaveform = ModulatingWaveform * 32767
	scipy.io.wavfile.write(dirname+"ModSignal.wav", PulseFilter['sample rate'], ModulatingWaveform.astype(np.int16))
	scipy.io.wavfile.write(dirname+"Baseband.wav", PulseFilter['sample rate'], Baseband.astype(np.float64))


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

		report_file.write(fo.GenInt16ArrayC('RRCFilter', PulseFilter['Taps'] * np.rint(65536 / np.sum(PulseFilter['Taps'])), PulseFilter['Oversample']))
		print(np.sum(PulseFilter['Taps']))


		report_file.write('\n')
		report_file.write(f'\nMax in table: {max(PulseFilter["FilterPatterns"])}')
		report_file.write(f'\nMin in table: {min(PulseFilter["FilterPatterns"])}\n')
		report_file.write(fo.GenInt16ArrayC(f'BPSKFilterPatterns', PulseFilter['FilterPatterns'], PulseFilter['Oversample']))
		report_file.write('\n\n')

		report_file.write(fo.GenInt16ArrayC(f'HalfBPSKFilterPatterns', PulseFilter['FilterPatterns'][0:len(PulseFilter['FilterPatterns'])//2], PulseFilter['Oversample']))


		report_file.write('\n\n')
		report_file.write('\n\n')

		report_file.write(fo.GenInt16ArrayC(f'Filter Window', np.rint(PulseFilter['FilterWindow'] * 32767), PulseFilter['Oversample']))
		report_file.write('\n\n')
		report_file.write('\n\n')

		report_file.write(fo.GenInt16ArrayC(f'Channel Filter', np.rint(channel_filter * 32767), PulseFilter['Oversample']))


		report_file.write('\n\n')
		report_file.write(fo.GenInt16ArrayC(f'SineTable', NCO['WaveTable'], 8))
		report_file.write(fo.GenInt16ArrayC(f'QuarterSineTable', NCO['WaveTable'][0:(len(NCO['WaveTable']) // 4) + 1], 8))

		report_file.write('\n\n')


		report_file.close()


	return
