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

def GetDPSKDemodulatorConfig(config, num):
	dpsk_demodulator = {}
	try:
		dpsk_demodulator['Enabled'] = config[f'DPSK Demodulator {num}'].getboolean('enabled')
	except:
		print(f'{argv[1]} [DPSK Demodulator {num}] \'enabled\' is missing or invalid')
		sys.exit(-2)
	try:
		dpsk_demodulator['AutoCorrelatorLag'] = int(config[f'DPSK Demodulator {num}']['autocorrelator lag'])
	except:
		print(f'{argv[1]} [DPSK Demodulator {num}] \'autocorrelator lag\' is missing or invalid')
		sys.exit(-2)
	try:
		dpsk_demodulator['CorrelatorShift'] = int(config[f'DPSK Demodulator {num}']['correlator shift'])
	except:
		print(f'{argv[1]} [DPSK Demodulator {num}] \'correlator shift\' is missing or invalid')
		sys.exit(-2)
	try:
		dpsk_demodulator['OutputFilter'] = strings.StringToIntArray(config[f'DPSK Demodulator {num}']['output filter taps'])
	except:
		print(f'{argv[1]} [DPSK Demodulator {num}] \'output filter taps\' is missing or invalid')
		sys.exit(-2)

	try:
		dpsk_demodulator['OutputFilterShift'] = int(config[f'DPSK Demodulator {num}']['output filter shift'])
	except:
		print(f'{argv[1]} [DPSK Demodulator {num}] \'output filter shift\' is missing or invalid')
		sys.exit(-2)
	return dpsk_demodulator

def FullProcess(state):
	argv = state['argv']
	config = state['config']
	
	print(f'Started DPSK2 process')
	print(f'Reading settings for Filter Decimator')
	FilterDecimator = input_filter.GetInputFilterConfig(state)
	FilterDecimator = demod.InitFilterDecimator(FilterDecimator)
	DPSKDemodulator = []
	DataSlicer = []
	DemodulatorCount = 0
	for DemodulatorNumber in range(4):
		DPSKDemodulator.append({})
		DataSlicer.append({})
		if config.has_section(f'DPSK Demodulator {DemodulatorNumber}'):
			print(f'Reading settings for DPSK Demodulator {DemodulatorNumber}')
			DemodulatorCount += 1
			DPSKDemodulator[DemodulatorNumber] = GetDPSKDemodulatorConfig(config, DemodulatorNumber)
			DPSKDemodulator[DemodulatorNumber]['InputSampleRate'] = FilterDecimator['OutputSampleRate']
			try:
				DataSlicer[DemodulatorNumber]['BitRate'] = int(config['Data Slicer 1']['slicer bit rate'])
			except:
				print(f'{argv[1]} [Data Slicer 1] \'slicer bit rate\' is missing or invalid')
				sys.exit(-2)
			try:
				DataSlicer[DemodulatorNumber]['PLLStep'] = int(config['Data Slicer 1']['pll clock step'])
			except:
				print(f'{argv[1]} [Data Slicer 1] \'pll clock step\' is missing or invalid')
				sys.exit(-2)
			try:
				DataSlicer[DemodulatorNumber]['PLLFeedbackGain'] = float(config['Data Slicer 1']['pll feedback gain'])
			except:
				print(f'{argv[1]} [Data Slicer 1] \'pll feedback gain\' is missing or invalid')
				sys.exit(-2)

			try:
				DataSlicer[DemodulatorNumber]['LoopFilter'] = strings.StringToIntArray(config['Data Slicer 1']['loop filter taps'])
			except:
				print(f'{argv[1]} [Data Slicer 1] \'loop filter taps\' is missing or invalid')
				sys.exit(-2)

			DataSlicer[DemodulatorNumber]['InputSampleRate'] = FilterDecimator['OutputSampleRate']
			DPSKDemodulator[DemodulatorNumber] = demod.InitDPSKDemod(DPSKDemodulator[DemodulatorNumber])
			DataSlicer[DemodulatorNumber] = demod.InitDataSlicer2(DataSlicer[DemodulatorNumber])

	DifferentialDecoderA = [{}]
	AX25Decoder = [{}]
	
	for index in range(1,DemodulatorCount+1):
		print(f'Initializing DifferentialDecoder {index} and AX25Decoder {index}')
		DifferentialDecoderA.append({})
		DifferentialDecoderA[index] = demod.InitDifferentialDecoder()
		AX25Decoder.append({})
		AX25Decoder[index] = demod.InitAX25Decoder()

	try:
		samplerate, audio = scipy.io.wavfile.read(argv[2])
		# Take two bits of resolution away
		audio = audio >> (16 - FilterDecimator['InputBitCount'])
	except:
		print('Unable to open wave file.')
		sys.exit(-2)

	print("Opened file. \r\nSample rate:", samplerate, "\r\nLength:", len(audio))


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

	total_packets = 0
	data = np.array([])
	filtered_signal_buffer = np.zeros(round(len(audio) / FilterDecimator['DecimationRate']) + 1)
	demod_sig_buffer1 = np.zeros(round(len(audio) / FilterDecimator['DecimationRate']) + 1)

	print(f'\nFiltering and decimating audio. ')
	FilterDecimator['FilterBuffer'] = audio
	FilterDecimator = demod.FilterDecimate(FilterDecimator)
	print(f'Done.')

	scipy.io.wavfile.write(dirname+"FilteredSignal.wav", FilterDecimator['OutputSampleRate'], FilterDecimator['FilterBuffer'].astype(np.int16))

	print(f'\nDemodulating audio. ')
	DPSKDemodulator[1]['InputBuffer'] = FilterDecimator['FilterBuffer']
	DPSKDemodulator[1] = demod.DemodulateDPSK3(DPSKDemodulator[1])
	print(f'Done.')

	print(f'\nSlicing, differential decoding, and AX25 decoding data. ')
	loop_count = len(DPSKDemodulator[1]['Result'])
	PhaseAccumulator = np.zeros(loop_count)
	PLLControl = np.zeros(loop_count)

	for index in range(loop_count):
		DataSlicer[1]['NewSample'] = DPSKDemodulator[1]['Result'][index]
		DataSlicer[1] = demod.ProgSliceData2(DataSlicer[1])
		#print(f'{DataSlicer[1]["PLLControl"]}')
		PhaseAccumulator[index] = DataSlicer[1]['PLLClock']
		PLLControl[index] = DataSlicer[1]['PLLControl']
		for data_bit in DataSlicer[1]['Result']:
			DifferentialDecoderA[1]['NewBit'] = data_bit
			DifferentialDecoderA[1] = demod.ProgDifferentialDecode(DifferentialDecoderA[1])
			AX25Decoder[1]['NewBit'] = DifferentialDecoderA[1]['Result']
			# if losing_index != 1:
			AX25Decoder[1] = demod.ProgDecodeAX25(AX25Decoder[1])
			if AX25Decoder[1]['OutputTrigger'] == True:
				AX25Decoder[1]['OutputTrigger'] = False
				# Check for unioqueness
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
	scipy.io.wavfile.write(dirname+"PhaseAccumulator.wav", FilterDecimator['OutputSampleRate'], PhaseAccumulator.astype(np.int16))
	scipy.io.wavfile.write(dirname+"PLLControl.wav", FilterDecimator['OutputSampleRate'], PLLControl.astype(np.int16))
	# scipy.io.wavfile.write(dirname+"DemodSignal2.wav", FilterDecimator['OutputSampleRate'], demod_sig_buffer2.astype(np.int16))
	#scipy.io.wavfile.write(dirname+"FilteredSignal.wav", FilterDecimator['OutputSampleRate'], filtered_signal_buffer.astype(np.int16))

	# Generate and save report file
	report_file_name = f'run{run_number}_report.txt'
	try:
		report_file = open(dirname + report_file_name, 'w+')
	except:
		print('Unable to create report file.')
	with report_file:
		report_file.write('# Command line: ')
		for argument in argv:
			report_file.write(f'{argument} ')
		report_file.write('\n#\n########## Begin Transcribed .ini file: ##########\n')
		try:
			ini_file = open(argv[1])
		except:
			report_file.write('Unable to open .ini file.')
		with ini_file:
			for character in ini_file:
				report_file.write(character)

		report_file.write('\n\n########## End Transcribed .ini file: ##########\n')

		report_file.write('\n\n# Demodulator performance:\n')
		report_file.write('\n')
		report_file.write(f'# Total packets: {total_packets}')
		report_file.write('\n')
		report_file.close()

	print('total packets: ', total_packets)
	print('made new directory: ', dirname)
	print('done')
	return