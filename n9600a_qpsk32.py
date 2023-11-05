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


def GetQPSKDemodConfig(config, num, id_string):
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

	key_string = "loop filter p"
	try:
		this['LoopFilter'][f'{key_string}'] = float(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)

	key_string = "loop filter i1"
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

	key_string = "loop filter gain"
	try:
		this['LoopFilter'][f'{key_string}'] = float(config[f'{id_string}{num}'][f'{key_string}'])
	except:
		print(f'{sys.argv[1]} [{id_string}{num}] \'{key_string}\' is missing or invalid')
		sys.exit(-2)



	return this

def InitQPSKDemod(this):
	this['NCO'] = nco.InitNCO(this['NCO'])
	this['LoopFilter'] = filters.InitIIR(this['LoopFilter'])
	this['LoopFilter2'] = this['LoopFilter'].copy()
	this['LoopFilter3'] = this['LoopFilter'].copy()
	this['LoopFilter2'] = filters.InitIIR(this['LoopFilter2'])
	this['LoopFilter3'] = filters.InitIIR(this['LoopFilter3'])
	this['I_LPF'] = filters.InitIIR(this['I_LPF'])
	this['Q_LPF'] = filters.InitIIR(this['Q_LPF'])
	this['BitAccumulator'] = 0
	return this

def DemodulateQPSK(this):
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
			# print(f'{index}')
			this['NCO'] = nco.UpdateNCO(this['NCO'])
			this['SineOutput'][index] = this['NCO']['Sine']
			this['CosineOutput'][index] = this['NCO']['Cosine']
			#branch_gain = 0.125
			# mix sample stream with NCO sinewave to create I branch
			this['I_Mixer'][index] = np.rint(sample * (this['NCO']['Cosine']) / 32768) # simulate 15 bit fractional integer multiplication
			# low-pass filter the mix product
			this['I_LPF'] = filters.UpdateIIR(this['I_LPF'], this['I_Mixer'][index])
			this['I_LPFOutput'][index] = this['I_LPF']['Output']

			# mix sample stream with NCO cosine to create Q branch
			this['Q_Mixer'][index] = np.rint(sample * (this['NCO']['Sine']) / 32768) # simulate 15 bit fractional integer multiplication
			# low-pass filter the mix product
			this['Q_LPF'] = filters.UpdateIIR(this['Q_LPF'], this['Q_Mixer'][index])
			this['Q_LPFOutput'][index] = this['Q_LPF']['Output']

			# Apply the signum function to I and Q.
			if this['I_LPFOutput'][index] > 0:
				I_Signum = 1
			else:
				I_Signum = -1

			if this['Q_LPFOutput'][index] > 0:
				Q_Signum = 1
			else:
				Q_Signum = -1

			# Cross-mix the branches
			I = np.rint(this['I_LPFOutput'][index] * Q_Signum)
			Q = np.rint(this['Q_LPFOutput'][index] * I_Signum)


			this['LoopMixer'][index] = I - Q


			#this['LoopMixer'][index] = np.rint(this['Q_LPF']['Output'] * this['I_LPF']['Output'] / 32768) # simulate 15 bit fractional integer multiplication
			this['LoopFilter'] = filters.UpdateIIR(this['LoopFilter'], this['LoopMixer'][index])
			this['LoopFilterOutput'][index] = np.rint(this['LoopFilter']['Output'])
			# scale the NCO control signal
			p = this['LoopFilterOutput'][index] * this['LoopFilter']['loop filter p']
			#integral += np.rint(this['LoopFilterOutput'][index] * this['LoopFilter']['loop filter i1'])
			integral += this['LoopFilterOutput'][index]
			#integral += this['LoopFilterOutput'][index]
			if abs(integral) > this['LoopFilter']['loop filter i max']:
				integral = 0
			this['LoopIntegral'][index] = integral
			this['NCO']['Control'] = np.rint((p + (np.rint(integral * this['LoopFilter']['loop filter i1']))) * this['LoopFilter']['loop filter gain'])

			this['NCOControlOutput'][index] = this['NCO']['Control']


			# save the data signal
			this['I_LPFOutput'][index] = this['I_LPF']['Output']





			this['PhaseAccumulator'][index] = this['NCO']['InPhase']



			index += 1
	return this

def SliceIQData(slicer):
	slicer['Midpoint'] = 0
	slicer['LastISample'] = 0
	slicer['LastQSample'] = 0
	slicer['IResult'] = np.zeros(int((len(slicer['IInput']) / slicer['Oversample']) * 1.1))
	slicer['QResult'] = np.zeros(int((len(slicer['QInput']) / slicer['Oversample']) * 1.1))
	output_index = 0
	LastISample = 0
	LastQSample = 0
	for input_index in range(len(slicer['IInput'])):
		ThisISample = slicer['IInput'][input_index]
		ThisQSample = slicer['QInput'][input_index]
		slicer['PLLClock'] += slicer['PLLStep']
		if slicer['PLLClock'] > ((slicer['PLLPeriod'] // 2) - 1):
			slicer['PLLClock'] -= slicer['PLLPeriod']
			slicer['IResult'][output_index] = ThisISample
			slicer['QResult'][output_index] = ThisQSample
			output_index += 1

		if (LastISample < 0 and ThisISample > 0) or (LastISample > 0 and ThisISample < 0) or (LastQSample < 0 and ThisQSample > 0) or (LastQSample > 0 and ThisQSample < 0):
			slicer['PLLClock'] = np.rint(slicer['Rate'] * slicer['PLLClock'])
		LastISample = ThisISample
		LastQSample = ThisQSample
	return slicer

def FullProcess(state):
	argv = state['argv']
	config = state['config']

	print(f'Started QPSK Demodulation process')

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

	PulseFilter = pulse_filter.GetRRCFilterConfig(state)
	#PulseFilter['sample rate'] = FilterDecimator['OutputSampleRate']
	PulseFilter = pulse_filter.InitRRCFilter(PulseFilter)

	print(f'Reading settings for QPSK Demodulators')
	QPSKDemodulator = []
	DemodulatorCount = 0
	id_string = "QPSK Demodulator "
	for DemodulatorNumber in range(4):
		QPSKDemodulator.append({})
		if config.has_section(f'{id_string}{DemodulatorNumber}'):
			print(f'Reading settings for {id_string}{DemodulatorNumber}')
			DemodulatorCount += 1
			QPSKDemodulator[DemodulatorNumber] = GetQPSKDemodConfig(config, DemodulatorNumber, id_string)
			QPSKDemodulator[DemodulatorNumber]['InputSampleRate'] = FilterDecimator['OutputSampleRate']
			QPSKDemodulator[DemodulatorNumber] = InitQPSKDemod(QPSKDemodulator[DemodulatorNumber])



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
	QPSKDemodulator[1]['InputBuffer'] = FilterDecimator['FilterBuffer']
	QPSKDemodulator[1] = DemodulateQPSK(QPSKDemodulator[1])

	FilteredIOutput = np.convolve(QPSKDemodulator[1]['I_LPFOutput'], np.rint(PulseFilter['Taps'] * 8191), 'valid') // 32768
	FilteredQOutput = np.convolve(QPSKDemodulator[1]['Q_LPFOutput'], np.rint(PulseFilter['Taps'] * 8191), 'valid') // 32768

	print(f'Done.')
	DataSlicer = {}
	DataSlicer['Oversample'] = 12
	DataSlicer['PLLStep'] = 128
	DataSlicer['PLLClock'] = 0
	DataSlicer['Rate'] = 0.9
	DataSlicer['PLLPeriod'] = DataSlicer['PLLStep'] * DataSlicer['Oversample']
	DataSlicer['IInput'] = FilteredIOutput[::4]
	DataSlicer['QInput'] = FilteredQOutput[::4]
	#DataSlicer['IInput'] = QPSKDemodulator[1]['I_LPFOutput'][::4]
	#DataSlicer['QInput'] = QPSKDemodulator[1]['Q_LPFOutput'][::4]

	DataSlicer = SliceIQData(DataSlicer)



	scipy.io.wavfile.write(dirname+"DemodSignal.wav", FilterDecimator['OutputSampleRate'], QPSKDemodulator[1]['Result'].astype(np.int16))
	scipy.io.wavfile.write(dirname+"LoopFilter.wav", FilterDecimator['OutputSampleRate'],QPSKDemodulator[1]['LoopFilterOutput'].astype(np.int16))
	scipy.io.wavfile.write(dirname+"I_LPF.wav", FilterDecimator['OutputSampleRate'], QPSKDemodulator[1]['I_LPFOutput'].astype(np.int16))
	scipy.io.wavfile.write(dirname+"Q_LPF.wav", FilterDecimator['OutputSampleRate'], QPSKDemodulator[1]['Q_LPFOutput'].astype(np.int16))

	scipy.io.wavfile.write(dirname+"I_Mixer.wav", FilterDecimator['OutputSampleRate'], QPSKDemodulator[1]['I_Mixer'].astype(np.int16))
	scipy.io.wavfile.write(dirname+"Q_Mixer.wav", FilterDecimator['OutputSampleRate'], QPSKDemodulator[1]['Q_Mixer'].astype(np.int16))

	#print(PulseFilter['Taps'])

	plt.figure()
	plt.subplot(221)
	#plt.plot(FilterDecimator['FilterBuffer'])
	plt.plot(FilterDecimator['FilterBuffer'])
	plt.plot(QPSKDemodulator[1]['I_LPFOutput'])
	plt.plot(QPSKDemodulator[1]['Q_LPFOutput'])
	plt.plot(FilterDecimator['EnvelopeBuffer'])
	#plt.plot(QPSKDemodulator[1]['Q_LPFOutput'])
	plt.title('I LPF Output')
	plt.legend(['Filtered Input','I_LPF Output','Q_LPF Output', 'Envelope'])
	plt.subplot(222)
	plt.plot(QPSKDemodulator[1]['LoopMixer'])
	plt.plot(QPSKDemodulator[1]['LoopIntegral'])
	plt.plot(QPSKDemodulator[1]['LoopFilterOutput'])
	plt.legend(['LoopMixer','LoopIntegral', 'LoopFilter'])
	plt.title('Loop Filter Output')
	plt.subplot(223)
	#plt.plot(QPSKDemodulator[1]['I_LPFOutput'])
	#plt.plot(QPSKDemodulator[1]['SamplePulse'])
	plt.plot(FilteredIOutput)
	plt.plot(FilteredQOutput)
	plt.legend(['Filtered I', 'Filtered Q'])
	plt.title('Filtered I/Q Output')
	plt.subplot(224)
	plt.plot(QPSKDemodulator[1]['NCOControlOutput'])
	plt.title('NCO Control')
	plt.show()

	plt.figure()
	plt.title('I/Q Demodulator Constellation')
	plt.scatter(QPSKDemodulator[1]['I_LPFOutput'], QPSKDemodulator[1]['Q_LPFOutput'], 0.1, c='tab:gray')
	plt.scatter(FilteredIOutput, FilteredQOutput, 0.1, c='tab:blue')
	plt.scatter(DataSlicer['IResult'], DataSlicer['QResult'], 5, c='tab:red')
	plt.legend(['Unfiltered I/Q', 'Filtered I/Q', 'Sample Points'])
	plt.xlim([-16000,16000])
	plt.ylim([-12000,12000])
	plt.show()

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



def ModulateRRC(state):
	argv = state['argv']
	config = state['config']
	print(f'Started Shaped QPSK Modulator process')
	print(f'Reading settings for RRC Pulse Shaping Filter')
	#PulseFilter = pulse_filter.GetGaussFilterConfig(state)
	#PulseFilter = pulse_filter.InitGaussFilter(PulseFilter)
	PulseFilter = pulse_filter.GetRRCFilterConfig(state)
	PulseFilter = pulse_filter.InitRRCFilter(PulseFilter)
	PulseFilter['SymbolMap'] = pulse_filter.GetSymbolMapConfig(state)
	NCO = nco.GetNCOConfig(config, 1, "TX NCO ")
	NCO = nco.InitNCO(NCO)

	PulseFilter = pulse_filter.GenPulseFilterPatterns(PulseFilter)
	print(max(PulseFilter['FilterPatterns']))
	plt.plot(PulseFilter['FilterPatterns'])
	plt.show()

	#BitStream = pulse_filter.ExpandSampleStream(state['InputData'], PulseFilter)
	BitStream = pulse_filter.BytesToSymbols(state['InputData'], PulseFilter)
	plt.figure()
	plt.plot(BitStream)
	plt.title('BitStream')
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
	plt.title('Modulating Waveform')
	plt.show()

	Baseband = np.zeros(len(ModulatingWaveform))
	i = 0
	for Amplitude in ModulatingWaveform:
		NCO = nco.UpdateNCO(NCO)
		Baseband[i] = Amplitude * NCO['Sine'] // PulseFilter['amplitude']
		i += 1

	plt.figure()
	plt.plot(Baseband)
	plt.title('Baseband')
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

		report_file.write(fo.GenInt16ArrayC('RRCFilter', PulseFilter['Taps'] * np.rint(65536 / np.sum(PulseFilter['Taps'])), PulseFilter['Oversample']))
		print(np.sum(PulseFilter['Taps']))


		report_file.write('\n')
		report_file.write(fo.GenInt16ArrayC(f'QPSKFilterPatterns', PulseFilter['FilterPatterns'], PulseFilter['Oversample']))
		report_file.write('\n\n')

		report_file.close()


	return