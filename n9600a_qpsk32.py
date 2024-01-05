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
	this['OutputSampleRate'] = this['InputSampleRate'] // this['intermediate decimation']
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
			this['NCO'] = nco.UpdateNCO(this['NCO'])
			this['SineOutput'][index] = this['NCO']['Sine']
			this['CosineOutput'][index] = this['NCO']['Cosine']
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

	print(f'Started QPSK Demodulation process')
	print(f'Reading settings for Filter Decimator')

	FilterDecimator = input_filter.GetInputFilterConfig(state)
	FilterDecimator = demod.InitFilterDecimator(FilterDecimator)

	DemodulatorNumber = 1
	id_string = "QPSK Demodulator "
	print(f'Reading settings for {id_string}{DemodulatorNumber}')
	QPSKDemodulator = GetQPSKDemodConfig(config, DemodulatorNumber, id_string)
	QPSKDemodulator['InputSampleRate'] = FilterDecimator['OutputSampleRate']
	QPSKDemodulator = InitQPSKDemod(QPSKDemodulator)

	ReceivePulseFilter = pulse_filter.GetRRCFilterConfig(state)
	ReceivePulseFilter['sample rate'] = QPSKDemodulator['OutputSampleRate']
	ReceivePulseFilter = pulse_filter.InitRRCFilter(ReceivePulseFilter)

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

	DataSlicer['InputSampleRate'] = QPSKDemodulator['OutputSampleRate'] // QPSKDemodulator['post decimation']
	print('Data Slicer Sample Rate: ', DataSlicer['InputSampleRate'])
	DataSlicer = demod.InitDataSlicer(DataSlicer)

	print(f'Initializing Descrambler and AX25Decoder')
	AX25Decoder = demod.InitAX25Decoder()
	Descrambler = {}
	Descrambler['Polynomial'] = int('0x21001',16) # G3RUH poly
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
	QPSKDemodulator['InputBuffer'] = FilterDecimator['FilterBuffer']
	QPSKDemodulator = DemodulateQPSK(QPSKDemodulator)
	print(f'Done.')

	ReceivePulseFilter['Taps'] = np.rint(ReceivePulseFilter['Taps'] * 8191)

	FilteredIOutput = np.convolve(QPSKDemodulator['I_LPFOutput'], ReceivePulseFilter['Taps'], 'valid') // 65536
	FilteredQOutput = np.convolve(QPSKDemodulator['Q_LPFOutput'], ReceivePulseFilter['Taps'], 'valid') // 65536

	FilteredIOutput = FilteredIOutput[0::QPSKDemodulator['post decimation']]
	FilteredQOutput = FilteredQOutput[0::QPSKDemodulator['post decimation']]

	print(f'\nSlicing, differential decoding, and AX25 decoding data. ')

	total_packets = 0

	DataSlicer['IInput'] = FilteredIOutput
	DataSlicer['QInput'] = FilteredQOutput

	DataSlicer = demod.SliceIQData(DataSlicer)

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
			filename = f'Packet-{total_packets}_CRC-{format(CRC,"#06x")}'
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
		scipy.io.wavfile.write(dirname+"DemodSignal.wav", FilterDecimator['OutputSampleRate'], QPSKDemodulator['Result'].astype(np.int16))
		scipy.io.wavfile.write(dirname+"LoopFilter.wav", FilterDecimator['OutputSampleRate'],QPSKDemodulator['LoopFilterOutput'].astype(np.int16))
		scipy.io.wavfile.write(dirname+"I_LPF.wav", FilterDecimator['OutputSampleRate'], QPSKDemodulator['I_LPFOutput'].astype(np.int16))
		scipy.io.wavfile.write(dirname+"Q_LPF.wav", FilterDecimator['OutputSampleRate'], QPSKDemodulator['Q_LPFOutput'].astype(np.int16))

		scipy.io.wavfile.write(dirname+"SamplePulse.wav", FilterDecimator['OutputSampleRate'], QPSKDemodulator['SamplePulse'].astype(np.int16))

		scipy.io.wavfile.write(dirname+"I_Mixer.wav", FilterDecimator['OutputSampleRate'], QPSKDemodulator['I_Mixer'].astype(np.int16))
		scipy.io.wavfile.write(dirname+"Q_Mixer.wav", FilterDecimator['OutputSampleRate'], QPSKDemodulator['Q_Mixer'].astype(np.int16))

	if state['plots'] == True:
		plt.figure()
		plt.subplot(221)
		plt.plot(FilterDecimator['FilterBuffer'])
		plt.plot(QPSKDemodulator['I_LPFOutput'])
		plt.plot(FilterDecimator['EnvelopeBuffer'])
		plt.title('I LPF Output')
		plt.legend(['Filtered Input','I_LPF Output','Envelope'])
		plt.subplot(222)
		plt.plot(QPSKDemodulator['LoopMixer'])
		plt.plot(QPSKDemodulator['LoopIntegral'])
		plt.plot(QPSKDemodulator['LoopFilterOutput'])
		plt.legend(['LoopMixer','LoopIntegral', 'LoopFilter'])
		plt.title('Loop Filter Output')
		plt.subplot(223)
		plt.plot(FilteredIOutput)
		plt.plot(FilteredQOutput)
		plt.title('Filtered I/Q Output')
		plt.subplot(224)
		plt.plot(QPSKDemodulator['NCOControlOutput'])
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
