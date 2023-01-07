import sys
import configparser
import struct
import scipy.io.wavfile
import numpy as np
import os
import n9600a_progdemod as demod
import format_output as fo

def StringToIntArray(input_string):
	working_string = ''
	result = np.array([])
	for character in input_string:
		if (character >= '0' and character <= '9') or ((character == '-') or (character == '.')):
			working_string += character
		else:
			if working_string != '':
				result = np.append(result, np.array([int(working_string)]))
				working_string = ''
	return result

def GetInputFilterConfig(config):
	FilterDecimator = {}
	try:
	 	FilterDecimator['Filter'] = StringToIntArray(config['Input Filter']['taps'])
	except:
		print(f'{sys.argv[1]} [Input Filter] \'taps\' is missing or invalid')
		sys.exit(-2)
	# Get input filter AGC option:
	try:
		FilterDecimator['InputAGCEnabled'] = config['Input Filter'].getboolean('agc enabled')
	except:
		print(f'{sys.argv[1]} [Input Filter] \'agc enabled\' is missing or invalid')
		sys.exit(-2)

	# Get Input Filter AGC Attack Rate:
	try:
		FilterDecimator['InputAGCAttackRate'] = int(config['Input Filter']['agc attack rate'])
	except:
		print(f'{sys.argv[1]} [Input Filter] \'agc attack rate\' is missing or invalid')
		sys.exit(-2)

	# Get Input Filter AGC Sustain Period:
	try:
		FilterDecimator['InputAGCSustainPeriod'] = int(config['Input Filter']['agc sustain period'])
	except:
		print(f'{sys.argv[1]} [Input Filter] \'agc sustain period\' is missing or invalid')
		sys.exit(-2)

	# Get Input Filter AGC Decay Rate:
	try:
		FilterDecimator['InputAGCDecayRate'] = int(config['Input Filter']['agc decay rate'])
	except:
		print(f'{sys.argv[1]} [Input Filter] \'agc decay rate\' is missing or invalid')
		sys.exit(-2)

	# Get Input Filter ADC bit count:
	try:
		FilterDecimator['InputBitCount'] = int(config['Input Filter']['input bit count'])
	except:
		print(f'{sys.argv[1]} [Input Filter] \'input bit count\' is missing or invalid')
		sys.exit(-2)
	# Get input filter Sample Rate:
	try:
		FilterDecimator['InputSampleRate'] = int(config['Input Filter']['sample rate'])
	except:
		print(f'{sys.argv[1]} [Input Filter] \'sample rate\' is missing or invalid')
		sys.exit(-2)

	# Get input filter Decimation Rate:
	try:
		FilterDecimator['DecimationRate'] = int(config['Input Filter']['decimation'])
	except:
		print(f'{sys.argv[1]} [Input Filter] \'decimation\' is missing or invalid')
		sys.exit(-2)
	return FilterDecimator

if len(sys.argv) < 3:
	print("Not enough arguments. Usage: py -3 n9600a_rx.py <ini file> <wav file>")
	sys.exit(-1)

def GetAFSKDemodulatorConfig(config, num):
	# Read settings for AFSK Demodulator
	afsk_demodulator = {}
	try:
		afsk_demodulator['Enabled'] = config[f'AFSK Demodulator {num}'].getboolean('enabled')
	except:
		print(f'{sys.argv[1]} [AFSK Demodulator {num}] \'enabled\' is missing or invalid')
		sys.exit(-2)

	try:
		afsk_demodulator['CorrelatorTapCount'] = int(config[f'AFSK Demodulator {num}']['correlator tap count'])
	except:
		print(f'{sys.argv[1]} [AFSK Demodulator {num}] \'correlator tap count\' is missing or invalid')
		sys.exit(-2)

	try:
		afsk_demodulator['MarkAmplitude'] = float(config[f'AFSK Demodulator {num}']['mark amplitude'])
	except:
		print(f'{sys.argv[1]} [AFSK Demodulator {num}] \'mark amplitude\' is missing or invalid')
		sys.exit(-2)

	try:
		afsk_demodulator['SpaceRatio'] = float(config[f'AFSK Demodulator {num}']['space gain'])
	except:
		print(f'{sys.argv[1]} [AFSK Demodulator {num}] \'space gain\' is missing or invalid')
		sys.exit(-2)

	try:
		afsk_demodulator['CorrelatorShift'] = int(config[f'AFSK Demodulator {num}']['correlator shift'])
	except:
		print(f'{sys.argv[1]} [AFSK Demodulator {num}] \'correlator shift\' is missing or invalid')
		sys.exit(-2)

	try:
		afsk_demodulator['MarkFreq'] = int(config[f'AFSK Demodulator {num}']['mark frequency'])
	except:
		print(f'{sys.argv[1]} [AFSK Demodulator {num}] \'mark frequency\' is missing or invalid')
		sys.exit(-2)

	try:
		afsk_demodulator['SpaceFreq'] = int(config[f'AFSK Demodulator {num}']['space frequency'])
	except:
		print(f'{sys.argv[1]} [AFSK Demodulator {num}] \'space frequency\' is missing or invalid')
		sys.exit(-2)

	try:
		afsk_demodulator['SquareSumBitCount'] = int(config[f'AFSK Demodulator {num}']['square sum bit count'])
	except:
		print(f'{sys.argv[1]} [AFSK Demodulator {num}] \'square sum bit count\' is missing or invalid')
		sys.exit(-2)

	try:
		afsk_demodulator['OutputFilter'] = StringToIntArray(config[f'AFSK Demodulator {num}']['output filter taps'])
	except:
		print(f'{sys.argv[1]} [AFSK Demodulator {num}] \'output filter taps\' is missing or invalid')
		sys.exit(-2)

	try:
		afsk_demodulator['OutputFilterShift'] = int(config[f'AFSK Demodulator {num}']['output filter shift'])
	except:
		print(f'{sys.argv[1]} [AFSK Demodulator {num}] \'output filter shift\' is missing or invalid')
		sys.exit(-2)

	try:
		afsk_demodulator['ChopAudio'] = config[f'AFSK Demodulator {num}'].getboolean('chop audio enable')
	except:
		print(f'{sys.argv[1]} [AFSK Demodulator {num}] \'chop audio enable\' is missing or invalid')
		sys.exit(-2)

	try:
		afsk_demodulator['SqrtBitCount'] = int(config[f'AFSK Demodulator {num}']['sqrt bit count'])
	except:
		print(f'{sys.argv[1]} [AFSK Demodulator {num}] \'sqrt bit count\' is missing or invalid')
		sys.exit(-2)
	return afsk_demodulator

def GetGFSKDemodulatorConfig(config, num):
	gfsk_demodulator = {}
	try:
		gfsk_demodulator['Enabled'] = config[f'GFSK Demodulator {num}'].getboolean('enabled')
	except:
		print(f'{sys.argv[1]} [GFSK Demodulator {num}] \'enabled\' is missing or invalid')
		sys.exit(-2)
	return gfsk_demodulator

def GetDPSKDemodulatorConfig(config, num):
	dpsk_demodulator = {}
	try:
		dpsk_demodulator['Enabled'] = config[f'DPSK Demodulator {num}'].getboolean('enabled')
	except:
		print(f'{sys.argv[1]} [DPSK Demodulator {num}] \'enabled\' is missing or invalid')
		sys.exit(-2)
	try:
		dpsk_demodulator['AutoCorrelatorLag'] = int(config[f'DPSK Demodulator {num}']['autocorrelator lag'])
	except:
		print(f'{sys.argv[1]} [DPSK Demodulator {num}] \'autocorrelator lag\' is missing or invalid')
		sys.exit(-2)
	try:
		dpsk_demodulator['CorrelatorShift'] = int(config[f'DPSK Demodulator {num}']['correlator shift'])
	except:
		print(f'{sys.argv[1]} [DPSK Demodulator {num}] \'correlator shift\' is missing or invalid')
		sys.exit(-2)
	try:
		dpsk_demodulator['OutputFilter'] = StringToIntArray(config[f'DPSK Demodulator {num}']['output filter taps'])
	except:
		print(f'{sys.argv[1]} [DPSK Demodulator {num}] \'output filter taps\' is missing or invalid')
		sys.exit(-2)

	try:
		dpsk_demodulator['OutputFilterShift'] = int(config[f'DPSK Demodulator {num}']['output filter shift'])
	except:
		print(f'{sys.argv[1]} [DPSK Demodulator {num}] \'output filter shift\' is missing or invalid')
		sys.exit(-2)
	return dpsk_demodulator

# read demodulator description from ini file:
config = configparser.ConfigParser()
try:
	config.read(sys.argv[1])
except:
	print(f'Unable to open ini file {sys.argv[1]}')
	sys.exit(-2)
print(f'Opened ini file {sys.argv[1]}')

print(f'Checking demodulator type')
try:
	DemodulatorType = config['General']['demodulator']
except:
	print(f'{sys.argv[1]} [General] \'demodulator\' is missing or invalid')
	sys.exit(-2)
print(f'Demodulator type is {DemodulatorType}')

print(f'Reading settings for Filter Decimator')
FilterDecimator = GetInputFilterConfig(config)
FilterDecimator = demod.InitFilterDecimator(FilterDecimator)

if DemodulatorType == 'afsk':
	# Initialize AFSK Demodulators
	AFSKDemodulator = []
	DataSlicer = []
	DemodulatorCount = 0
	for DemodulatorNumber in range(4):
		AFSKDemodulator.append({})
		DataSlicer.append({})
		if config.has_section(f'AFSK Demodulator {DemodulatorNumber}'):
			print(f'Reading settings for AFSK Demodulator {DemodulatorNumber}')
			DemodulatorCount += 1
			AFSKDemodulator[DemodulatorNumber] = GetAFSKDemodulatorConfig(config, DemodulatorNumber)
			AFSKDemodulator[DemodulatorNumber]['InputSampleRate'] = FilterDecimator['OutputSampleRate']
			try:
				DataSlicer[DemodulatorNumber]['BitRate'] = int(config['Data Slicer 1']['slicer bit rate'])
			except:
				print(f'{sys.argv[1]} [Data Slicer 1] \'slicer bit rate\' is missing or invalid')
				sys.exit(-2)
			try:
				DataSlicer[DemodulatorNumber]['Rate'] = float(config['Data Slicer 1']['slicer lock rate'])
			except:
				print(f'{sys.argv[1]} [Data Slicer 1] \'slicer lock rate\' is missing or invalid')
				sys.exit(-2)

			DataSlicer[DemodulatorNumber]['InputSampleRate'] = FilterDecimator['OutputSampleRate']
			AFSKDemodulator[DemodulatorNumber] = demod.InitAFSKDemod(AFSKDemodulator[DemodulatorNumber])
			DataSlicer[DemodulatorNumber] = demod.InitDataSlicer(DataSlicer[DemodulatorNumber])

	DifferentialDecoder = [{}]
	AX25Decoder = [{}]
	for index in range(1,DemodulatorCount+1):
		print(f'Initializing DifferentialDecoder {index} and AX25Decoder {index}')
		DifferentialDecoder.append({})
		AX25Decoder.append({})
		DifferentialDecoder[index] = demod.InitDifferentialDecoder()
		AX25Decoder[index] = demod.InitAX25Decoder()

	try:
		samplerate, audio = scipy.io.wavfile.read(sys.argv[2])
		# Take two bits of resolution away
		audio = audio >> (16 - FilterDecimator['InputBitCount'])
	except:
		print('Unable to open wave file.')
		sys.exit(-2)

	print("Opened file. \r\nSample rate:", samplerate, "\r\nLength:", len(audio))

	index1 = 0
	index2 = 0
	index3 = 0
	envelope_index = 0

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
	duplicate_packets = 0

	data = np.array([])
	filtered_signal_buffer = np.zeros(round(len(audio) / FilterDecimator['DecimationRate']) + 1)
	demod_sig_buffer1 = np.zeros(round(len(audio) / FilterDecimator['DecimationRate']) + 1)
	demod_sig_buffer2 = np.zeros(round(len(audio) / FilterDecimator['DecimationRate']) + 1)
	demod_sig_buffer3 = np.zeros(round(len(audio) / FilterDecimator['DecimationRate']) + 1)

	print(f'\nFiltering and decimating audio. ')
	FilterDecimator['FilterBuffer'] = audio
	FilterDecimator = demod.FilterDecimate(FilterDecimator)
	print(f'Done.')

	scipy.io.wavfile.write(dirname+"FilteredSignal.wav", FilterDecimator['OutputSampleRate'], FilterDecimator['FilterBuffer'].astype(np.int16))

	print(f'\nDemodulating audio. ')
	AFSKDemodulator[1]['CorrelatorBuffer'] = FilterDecimator['FilterBuffer']
	AFSKDemodulator[2]['CorrelatorBuffer'] = FilterDecimator['FilterBuffer']
	AFSKDemodulator[3]['CorrelatorBuffer'] = FilterDecimator['FilterBuffer']
	AFSKDemodulator[1] = demod.DemodulateAFSK(AFSKDemodulator[1])
	AFSKDemodulator[2] = demod.DemodulateAFSK(AFSKDemodulator[2])
	AFSKDemodulator[3] = demod.DemodulateAFSK(AFSKDemodulator[3])
	print(f'Done.')

	print(f'\nSlicing, differential decoding, and AX25 decoding data. ')
	loop_count = np.min([len(AFSKDemodulator[1]['Result']), len(AFSKDemodulator[2]['Result']), len(AFSKDemodulator[3]['Result'])])
	for index in range(loop_count):
		trials = [DataSlicer[1]['AvgPhaseError'], DataSlicer[2]['AvgPhaseError'], DataSlicer[3]['AvgPhaseError']]
		winning_value = 100
		winning_index = 0
		losing_value = 0
		losing_index = 0
		for test in range(3):
			# if trials[test] > winning_value:
			# 	winning_index = test
			# 	winning_value = trials[test]
			if trials[test] > losing_value:
				losing_index = test
				losing_value = trials[test]
		losing_index += 1
		DataSlicer[1]['NewSample'] = AFSKDemodulator[1]['Result'][index]
		DataSlicer[1] = demod.ProgSliceData(DataSlicer[1])
		for data_bit in DataSlicer[1]['Result']:
			DifferentialDecoder[1]['NewBit'] = data_bit
			DifferentialDecoder[1] = demod.ProgDifferentialDecode(DifferentialDecoder[1])
			AX25Decoder[1]['NewBit'] = DifferentialDecoder[1]['Result']
			# if losing_index != 1:
			AX25Decoder[1] = demod.ProgDecodeAX25(AX25Decoder[1])
			if AX25Decoder[1]['OutputTrigger'] == True:
				AX25Decoder[1]['OutputTrigger'] = False
				# Check for unioqueness
				if ((AX25Decoder[2]['CRCAge'] > 10) and (AX25Decoder[3]['CRCAge'] > 10)) or ((AX25Decoder[1]['CRC'][0] != AX25Decoder[2]['CRC'][0]) and (AX25Decoder[1]['CRC'] != AX25Decoder[3]['CRC'])):
					total_packets += 1
					CRC = AX25Decoder[1]['CRC'][0]
					decodernum = '1'
					filename = f'Packet-{total_packets}_CRC-{format(CRC,"#06x")}_decoder-{decodernum}_Index-{index1}'
					print(f'{dirname+filename} {DataSlicer[1]["AvgPhaseError"]} {DataSlicer[2]["AvgPhaseError"]} {DataSlicer[3]["AvgPhaseError"]} ')
					# try:
					# 	bin_file = open(dirname + filename + '.bin', '+wb')
					# except:
					# 	pass
					# with bin_file:
					# 	for byte in AX25Decoder[2]['Output']:
					# 		bin_file.write(byte.astype('uint8'))
					# 	bin_file.close()
				else:
					if (AX25Decoder[1]['CRC'][0] == AX25Decoder[2]['CRC'][0]) and AX25Decoder[2]['CRCAge'] <= 10:
						duplicate_packets += 1
						if AX25Decoder[1]['UniqueFlag'] == True:
							AX25Decoder[1]['UniquePackets'] -= 1
							AX25Decoder[1]['UniqueFlag'] = False
						if AX25Decoder[2]['UniqueFlag'] == True:
							AX25Decoder[2]['UniquePackets'] -= 1
							AX25Decoder[2]['UniqueFlag'] = False
						print('Decoder 1 Duplicate 2, bit delay: ', AX25Decoder[2]['CRCAge'])
					if (AX25Decoder[1]['CRC'][0] == AX25Decoder[3]['CRC'][0]) and AX25Decoder[3]['CRCAge'] <= 10:
						duplicate_packets += 1
						if AX25Decoder[1]['UniqueFlag'] == True:
							AX25Decoder[1]['UniquePackets'] -= 1
							AX25Decoder[1]['UniqueFlag'] = False
						if AX25Decoder[3]['UniqueFlag'] == True:
							AX25Decoder[3]['UniquePackets'] -= 1
							AX25Decoder[3]['UniqueFlag'] = False
						print('Decoder 1 Duplicate 3, bit delay: ', AX25Decoder[3]['CRCAge'])

		DataSlicer[2]['NewSample'] = AFSKDemodulator[2]['Result'][index]
		DataSlicer[2] = demod.ProgSliceData(DataSlicer[2])
		for data_bit in DataSlicer[2]['Result']:
			DifferentialDecoder[2]['NewBit'] = data_bit
			DifferentialDecoder[2] = demod.ProgDifferentialDecode(DifferentialDecoder[2])
			AX25Decoder[2]['NewBit'] = DifferentialDecoder[2]['Result']
			# if losing_index != 2:
			AX25Decoder[2] = demod.ProgDecodeAX25(AX25Decoder[2])
			if AX25Decoder[2]['OutputTrigger'] == True:
				AX25Decoder[2]['OutputTrigger'] = False
				# Check for uniqueness
				if ((AX25Decoder[1]['CRCAge'] > 10) and (AX25Decoder[3]['CRCAge'] > 10)) or ((AX25Decoder[2]['CRC'][0] != AX25Decoder[1]['CRC'][0]) and (AX25Decoder[2]['CRC'] != AX25Decoder[3]['CRC'])):
					total_packets += 1
					CRC = AX25Decoder[2]['CRC'][0]
					decodernum = '2'
					filename = f'Packet-{total_packets}_CRC-{format(CRC,"#06x")}_decoder-{decodernum}_Index-{index1}'
					print(f'{dirname+filename} {DataSlicer[1]["AvgPhaseError"]} {DataSlicer[2]["AvgPhaseError"]} {DataSlicer[3]["AvgPhaseError"]} ')
					# try:
					# 	bin_file = open(dirname + filename + '.bin', '+wb')
					# except:
					# 	pass
					# with bin_file:
					# 	for byte in AX25Decoder[2]['Output']:
					# 		bin_file.write(byte.astype('uint8'))
					# 	bin_file.close()
				else:
					if (AX25Decoder[2]['CRC'][0] == AX25Decoder[1]['CRC'][0]) and AX25Decoder[1]['CRCAge'] <= 10:
						duplicate_packets += 1
						if AX25Decoder[2]['UniqueFlag'] == True:
							AX25Decoder[2]['UniquePackets'] -= 1
							AX25Decoder[2]['UniqueFlag'] = False
						if AX25Decoder[1]['UniqueFlag'] == True:
							AX25Decoder[1]['UniquePackets'] -= 1
							AX25Decoder[1]['UniqueFlag'] = False
						print('Decoder 2 Duplicate 1, bit delay: ', AX25Decoder[1]['CRCAge'])
					if (AX25Decoder[2]['CRC'][0] == AX25Decoder[3]['CRC'][0]) and AX25Decoder[3]['CRCAge'] <= 10:
						duplicate_packets += 1
						if AX25Decoder[2]['UniqueFlag'] == True:
							AX25Decoder[2]['UniquePackets'] -= 1
							AX25Decoder[2]['UniqueFlag'] = False
						if AX25Decoder[3]['UniqueFlag'] == True:
							AX25Decoder[3]['UniquePackets'] -= 1
							AX25Decoder[3]['UniqueFlag'] = False
						print('Decoder 2 Duplicate 3, bit delay: ', AX25Decoder[3]['CRCAge'])

		DataSlicer[3]['NewSample'] = AFSKDemodulator[3]['Result'][index]
		DataSlicer[3] = demod.ProgSliceData(DataSlicer[3])
		for data_bit in DataSlicer[3]['Result']:
			DifferentialDecoder[3]['NewBit'] = data_bit
			DifferentialDecoder[3] = demod.ProgDifferentialDecode(DifferentialDecoder[3])
			AX25Decoder[3]['NewBit'] = DifferentialDecoder[3]['Result']
			# if losing_index != 3:
			AX25Decoder[3] = demod.ProgDecodeAX25(AX25Decoder[3])
			if AX25Decoder[3]['OutputTrigger'] == True:
				AX25Decoder[3]['OutputTrigger'] = False
				# Check for uniqueness
				if ((AX25Decoder[1]['CRCAge'] > 10) and (AX25Decoder[2]['CRCAge'] > 10)) or ((AX25Decoder[3]['CRC'][0] != AX25Decoder[1]['CRC'][0]) and (AX25Decoder[3]['CRC'] != AX25Decoder[2]['CRC'])):
					total_packets += 1
					CRC = AX25Decoder[2]['CRC'][0]
					decodernum = '3'
					filename = f'Packet-{total_packets}_CRC-{format(CRC,"#06x")}_decoder-{decodernum}_Index-{index1}'
					print(f'{dirname+filename} {DataSlicer[1]["AvgPhaseError"]} {DataSlicer[2]["AvgPhaseError"]} {DataSlicer[3]["AvgPhaseError"]} ')
					# try:
					# 	bin_file = open(dirname + filename + '.bin', '+wb')
					# except:
					# 	pass
					# with bin_file:
					# 	for byte in AX25Decoder[2]['Output']:
					# 		bin_file.write(byte.astype('uint8'))
					# 	bin_file.close()
				else:
					if (AX25Decoder[3]['CRC'][0] == AX25Decoder[1]['CRC'][0]) and AX25Decoder[1]['CRCAge'] <= 10:
						duplicate_packets += 1
						if AX25Decoder[3]['UniqueFlag'] == True:
							AX25Decoder[3]['UniquePackets'] -= 1
							AX25Decoder[3]['UniqueFlag'] = False
						if AX25Decoder[1]['UniqueFlag'] == True:
							AX25Decoder[1]['UniquePackets'] -= 1
							AX25Decoder[1]['UniqueFlag'] = False
						print('Decoder 3 Duplicate 1, bit delay: ', AX25Decoder[1]['CRCAge'])
					if (AX25Decoder[3]['CRC'][0] == AX25Decoder[2]['CRC'][0]) and AX25Decoder[2]['CRCAge'] <= 10:
						duplicate_packets += 1
						if AX25Decoder[3]['UniqueFlag'] == True:
							AX25Decoder[3]['UniquePackets'] -= 1
							AX25Decoder[3]['UniqueFlag'] = False
						if AX25Decoder[2]['UniqueFlag'] == True:
							AX25Decoder[2]['UniquePackets'] -= 1
							AX25Decoder[2]['UniqueFlag'] = False
						print('Decoder 3 Duplicate 2, bit delay: ', AX25Decoder[2]['CRCAge'])

	# scipy.io.wavfile.write(dirname+"DemodSignal1.wav", FilterDecimator['OutputSampleRate'], demod_sig_buffer1.astype(np.int16))
	# scipy.io.wavfile.write(dirname+"DemodSignal2.wav", FilterDecimator['OutputSampleRate'], demod_sig_buffer2.astype(np.int16))
	# scipy.io.wavfile.write(dirname+"FilteredSignal.wav", FilterDecimator['OutputSampleRate'], filtered_signal_buffer.astype(np.int16))

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

		report_file.write(fo.GenInt16ArrayC('MarkCorCos1', AFSKDemodulator[1]['MarkCOS'], 32))
		report_file.write(fo.GenInt16ArrayC('MarkCorSin1', AFSKDemodulator[1]['MarkSIN'], 32))
		report_file.write(fo.GenInt16ArrayC('SpaceCorCos1', AFSKDemodulator[1]['SpaceCOS'], 32))
		report_file.write(fo.GenInt16ArrayC('SpaceCorSin1', AFSKDemodulator[1]['SpaceSIN'], 32))
		report_file.write('\n')
		report_file.write(fo.GenInt16ArrayC('MarkCorCos2', AFSKDemodulator[2]['MarkCOS'], 32))
		report_file.write(fo.GenInt16ArrayC('MarkCorSin2', AFSKDemodulator[2]['MarkSIN'], 32))
		report_file.write(fo.GenInt16ArrayC('SpaceCorCos2', AFSKDemodulator[2]['SpaceCOS'], 32))
		report_file.write(fo.GenInt16ArrayC('SpaceCorSin2', AFSKDemodulator[2]['SpaceSIN'], 32))
		report_file.write('\n')
		report_file.write(fo.GenInt16ArrayC('MarkCorCos3', AFSKDemodulator[3]['MarkCOS'], 32))
		report_file.write(fo.GenInt16ArrayC('MarkCorSin3', AFSKDemodulator[3]['MarkSIN'], 32))
		report_file.write(fo.GenInt16ArrayC('SpaceCorCos3', AFSKDemodulator[3]['SpaceCOS'], 32))
		report_file.write(fo.GenInt16ArrayC('SpaceCorSin3', AFSKDemodulator[3]['SpaceSIN'], 32))

		report_file.write('\nMaximum correlation values:')
		tstep = 1.0 / AFSKDemodulator[1]['InputSampleRate']
		time = np.arange(0, tstep * AFSKDemodulator[1]['CorrelatorTapCount'], tstep)
		max_mark_cos_input = np.rint(24576 * (np.cos(2 * AFSKDemodulator[1]['MarkFreq'] * np.pi * time)))
		max_mark_sin_input = np.rint(24576 * (np.sin(2 * AFSKDemodulator[1]['MarkFreq'] * np.pi * time)))
		max_space_cos_input = np.rint(24576 * (np.cos(2 * AFSKDemodulator[1]['SpaceFreq'] * np.pi * time)))
		max_space_sin_input = np.rint(24576 * (np.sin(2 * AFSKDemodulator[1]['SpaceFreq'] * np.pi * time)))

		mark_sin_1 = np.convolve(max_mark_sin_input, AFSKDemodulator[1]["MarkSIN"], "valid") / pow(2, (16 + AFSKDemodulator[1]['CorrelatorShift']))
		mark_cos_1 = np.convolve(max_mark_sin_input, AFSKDemodulator[1]["MarkCOS"], "valid") / pow(2, (16 + AFSKDemodulator[1]['CorrelatorShift']))
		mark_square_sum_1 = mark_sin_1**2 + mark_cos_1**2

		space_sin_1 = np.convolve(max_space_sin_input, AFSKDemodulator[1]["SpaceSIN"], "valid") / pow(2, (16 + AFSKDemodulator[1]['CorrelatorShift']))
		space_cos_1 = np.convolve(max_space_sin_input, AFSKDemodulator[1]["SpaceCOS"], "valid") / pow(2, (16 + AFSKDemodulator[1]['CorrelatorShift']))
		space_square_sum_1 = space_sin_1**2 + space_cos_1**2

		mark_sin_2 = np.convolve(max_mark_sin_input, AFSKDemodulator[2]["MarkSIN"], "valid") / pow(2, (16 + AFSKDemodulator[2]['CorrelatorShift']))
		mark_cos_2 = np.convolve(max_mark_sin_input, AFSKDemodulator[2]["MarkCOS"], "valid") / pow(2, (16 + AFSKDemodulator[2]['CorrelatorShift']))
		mark_square_sum_2 = mark_sin_2**2 + mark_cos_2**2

		space_sin_2 = np.convolve(max_space_sin_input, AFSKDemodulator[2]["SpaceSIN"], "valid") / pow(2, (16 + AFSKDemodulator[2]['CorrelatorShift']))
		space_cos_2 = np.convolve(max_space_sin_input, AFSKDemodulator[2]["SpaceCOS"], "valid") / pow(2, (16 + AFSKDemodulator[2]['CorrelatorShift']))
		space_square_sum_2 = space_sin_2**2 + space_cos_2**2

		mark_sin_3 = np.convolve(max_mark_sin_input, AFSKDemodulator[3]["MarkSIN"], "valid") / pow(2, (16 + AFSKDemodulator[2]['CorrelatorShift']))
		mark_cos_3 = np.convolve(max_mark_sin_input, AFSKDemodulator[3]["MarkCOS"], "valid") / pow(2, (16 + AFSKDemodulator[2]['CorrelatorShift']))
		mark_square_sum_3 = mark_sin_2**2 + mark_cos_2**2

		space_sin_3 = np.convolve(max_space_sin_input, AFSKDemodulator[3]["SpaceSIN"], "valid") / pow(2, (16 + AFSKDemodulator[3]['CorrelatorShift']))
		space_cos_3 = np.convolve(max_space_sin_input, AFSKDemodulator[3]["SpaceCOS"], "valid") / pow(2, (16 + AFSKDemodulator[3]['CorrelatorShift']))
		space_square_sum_3 = space_sin_3**2 + space_cos_3**2

		report_file.write('\n')
		report_file.write(f'\nMark Sin Correlator 1: {int(np.rint(mark_sin_1[0]))}')
		report_file.write(f'\nMark Cos Correlator 1: {int(np.rint(mark_cos_1[0]))}')
		report_file.write(f'\nMark Square Sum Correlator 1: {int(np.rint(mark_square_sum_1[0]))}')
		report_file.write(f'\nMark Square Sum / Square Scale: {int(np.rint(mark_square_sum_1[0] / AFSKDemodulator[1]["SquareScale"]))}')
		# report_file.write(f'\n Log2(Mark Square Sum / Sqrt Table Size): {np.log2(mark_square_sum_1[0] // 2**AFSKDemodulator[1]["SqrtBitCount"]):.2f}')

		report_file.write('\n')
		report_file.write(f'\nSpace Sin Correlator 1: {int(np.rint(space_sin_1[0]))}')
		report_file.write(f'\nSpace Cos Correlator 1: {int(np.rint(space_cos_1[0]))}')
		report_file.write(f'\nSpace Square Sum Correlator 1: {int(np.rint(space_square_sum_1[0]))}')
		report_file.write(f'\nMark Square Sum / Square Scale: {int(np.rint(space_square_sum_1[0] / AFSKDemodulator[1]["SquareScale"]))}')
		# report_file.write(f'\n Log2(Space Square Sum / Sqrt Table Size): {np.log2(space_square_sum_1[0] // 2**AFSKDemodulator[1]["SqrtBitCount"]):.2f}')

		report_file.write('\n')
		report_file.write(f'\nMark Sin Correlator 2: {int(np.rint(mark_sin_2[0]))}')
		report_file.write(f'\nMark Cos Correlator 2: {int(np.rint(mark_cos_2[0]))}')
		report_file.write(f'\nMark Square Sum Correlator 2: {int(np.rint(mark_square_sum_2[0]))}')
		report_file.write(f'\nMark Square Sum / Square Scale: {int(np.rint(mark_square_sum_2[0] / AFSKDemodulator[2]["SquareScale"]))}')
		# report_file.write(f'\n Log2(Mark Square Sum / Sqrt Table Size): {np.log2(mark_square_sum_2[0] // 2**AFSKDemodulator[2]["SqrtBitCount"]):.2f}')

		report_file.write('\n')
		report_file.write(f'\nSpace Sin Correlator 2: {int(np.rint(space_sin_2[0]))}')
		report_file.write(f'\nSpace Cos Correlator 2: {int(np.rint(space_cos_2[0]))}')
		report_file.write(f'\nSpace Square Sum Correlator 2: {int(np.rint(space_square_sum_2[0]))}')
		report_file.write(f'\nMark Square Sum / Square Scale: {int(np.rint(space_square_sum_2[0] / AFSKDemodulator[2]["SquareScale"]))}')
		# report_file.write(f'\n Log2(Space Square Sum / Sqrt Table Size): {np.log2(space_square_sum_2[0] // 2**AFSKDemodulator[2]["SqrtBitCount"]):.2f}')


		report_file.write('\n\n')
		report_file.write(f'Square Scale 1: {AFSKDemodulator[1]["SquareScale"]}')

		report_file.write('\n')
		report_file.write(f'Square Scale 2: {AFSKDemodulator[2]["SquareScale"]}')

		report_file.write('\n\n')
		report_file.write(fo.GenInt16ArrayC(f'SquareRoot{AFSKDemodulator[1]["SqrtBitCount"]}', AFSKDemodulator[1]['SqrtTable'], 16))

		report_file.write('\n\n# Demodulator performance:\n')
		report_file.write('\n')
		report_file.write(f'# Total packets: {total_packets}')
		report_file.write('\n')
		report_file.write(f'# Duplicate packets: {duplicate_packets}')
		report_file.write('\n')
		report_file.write(f'# Demodulator 1 unique packets: {AX25Decoder[1]["UniquePackets"] // 2}')
		report_file.write('\n')
		report_file.write(f'# Demodulator 1 total packets: {AX25Decoder[1]["PacketCount"]}')
		report_file.write('\n')
		report_file.write(f'# Demodulator 2 unique packets: {AX25Decoder[2]["UniquePackets"] // 2}')
		report_file.write('\n')
		report_file.write(f'# Demodulator 2 total packets: {AX25Decoder[2]["PacketCount"]}')
		report_file.write('\n')
		report_file.write(f'# Demodulator 3 unique packets: {AX25Decoder[3]["UniquePackets"] // 2}')
		report_file.write('\n')
		report_file.write(f'# Demodulator 3 total packets: {AX25Decoder[3]["PacketCount"]}')
		report_file.write('\n')
		report_file.close()

	print('total packets: ', total_packets)
	print('duplicate_packets: ', duplicate_packets)
	print('Decoder 1 unique packets: ', AX25Decoder[1]['UniquePackets'] // 2)
	print('Decoder 1 total packets: ', AX25Decoder[1]['PacketCount'])
	print('Decoder 2 unique packets: ', AX25Decoder[2]['UniquePackets'] // 2)
	print('Decoder 2 total packets: ', AX25Decoder[2]['PacketCount'])
	print('Decoder 3 unique packets: ', AX25Decoder[3]['UniquePackets'] // 2)
	print('Decoder 3 total packets: ', AX25Decoder[3]['PacketCount'])
	print('made new directory: ', dirname)
	print('done')

elif DemodulatorType == 'dpsk':
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
				print(f'{sys.argv[1]} [Data Slicer 1] \'slicer bit rate\' is missing or invalid')
				sys.exit(-2)
			try:
				DataSlicer[DemodulatorNumber]['Rate'] = float(config['Data Slicer 1']['slicer lock rate'])
			except:
				print(f'{sys.argv[1]} [Data Slicer 1] \'slicer lock rate\' is missing or invalid')
				sys.exit(-2)

			DataSlicer[DemodulatorNumber]['InputSampleRate'] = FilterDecimator['OutputSampleRate']
			DPSKDemodulator[DemodulatorNumber] = demod.InitDPSKDemod(DPSKDemodulator[DemodulatorNumber])
			DataSlicer[DemodulatorNumber] = demod.InitDataSlicer(DataSlicer[DemodulatorNumber])

	DifferentialDecoder = [{}]
	AX25Decoder = [{}]
	for index in range(1,DemodulatorCount+1):
		print(f'Initializing DifferentialDecoder {index} and AX25Decoder {index}')
		DifferentialDecoder.append({})
		AX25Decoder.append({})
		DifferentialDecoder[index] = demod.InitDifferentialDecoder()
		AX25Decoder[index] = demod.InitAX25Decoder()

	try:
		samplerate, audio = scipy.io.wavfile.read(sys.argv[2])
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
	for index in range(loop_count):
		DataSlicer[1]['NewSample'] = DPSKDemodulator[1]['Result'][index]
		DataSlicer[1] = demod.ProgSliceData(DataSlicer[1])
		for data_bit in DataSlicer[1]['Result']:
			DifferentialDecoder[1]['NewBit'] = data_bit
			DifferentialDecoder[1] = demod.ProgDifferentialDecode(DifferentialDecoder[1])
			AX25Decoder[1]['NewBit'] = DifferentialDecoder[1]['Result']
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

		report_file.write('\n\n# Demodulator performance:\n')
		report_file.write('\n')
		report_file.write(f'# Total packets: {total_packets}')
		report_file.write('\n')
		report_file.close()

	print('total packets: ', total_packets)
	print('made new directory: ', dirname)
	print('done')
elif DemodulatorType == 'gfsk':
	GFSKDemodulator = []
	DataSlicer = []
	DemodulatorCount = 0
	for DemodulatorNumber in range(4):
		GFSKDemodulator.append({})
		DataSlicer.append({})
		if config.has_section(f'GFSK Demodulator {DemodulatorNumber}'):
			print(f'Reading settings for GFSK Demodulator {DemodulatorNumber}')
			DemodulatorCount += 1
			GFSKDemodulator[DemodulatorNumber] = GetGFSKDemodulatorConfig(config, DemodulatorNumber)
			GFSKDemodulator[DemodulatorNumber]['InputSampleRate'] = FilterDecimator['OutputSampleRate']
			try:
				DataSlicer[DemodulatorNumber]['BitRate'] = int(config[f'Data Slicer {DemodulatorNumber}']['slicer bit rate'])
			except:
				print(f'{sys.argv[1]} [Data Slicer {DemodulatorNumber}] \'slicer bit rate\' is missing or invalid')
				sys.exit(-2)
			try:
				DataSlicer[DemodulatorNumber]['Rate'] = float(config[f'Data Slicer {DemodulatorNumber}']['slicer lock rate'])
			except:
				print(f'{sys.argv[1]} [Data Slicer {DemodulatorNumber}] \'slicer lock rate\' is missing or invalid')
				sys.exit(-2)

			DataSlicer[DemodulatorNumber]['InputSampleRate'] = FilterDecimator['OutputSampleRate']
			GFSKDemodulator[DemodulatorNumber] = demod.InitGFSKDemod(GFSKDemodulator[DemodulatorNumber])
			DataSlicer[DemodulatorNumber] = demod.InitDataSlicer(DataSlicer[DemodulatorNumber])

	DifferentialDecoder = [{}]
	AX25Decoder = [{}]
	for index in range(1,DemodulatorCount+1):
		print(f'Initializing DifferentialDecoder {index} and AX25Decoder {index}')
		DifferentialDecoder.append({})
		AX25Decoder.append({})
		DifferentialDecoder[index] = demod.InitDifferentialDecoder()
		AX25Decoder[index] = demod.InitAX25Decoder()

	try:
		samplerate, audio = scipy.io.wavfile.read(sys.argv[2])
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
	GFSKDemodulator[1]['InputBuffer'] = FilterDecimator['FilterBuffer']
	GFSKDemodulator[1] = demod.DemodulateGFSK(GFSKDemodulator[1])
	print(f'Done.')

	print(f'\nSlicing, differential decoding, and AX25 decoding data. ')
	loop_count = len(GFSKDemodulator[1]['Result'])
	for index in range(loop_count):
		DataSlicer[1]['NewSample'] = GFSKDemodulator[1]['Result'][index]
		DataSlicer[1] = demod.ProgSliceData(DataSlicer[1])
		for data_bit in DataSlicer[1]['Result']:
			DifferentialDecoder[1]['NewBit'] = data_bit
			DifferentialDecoder[1] = demod.ProgDifferentialDecode(DifferentialDecoder[1])
			AX25Decoder[1]['NewBit'] = DifferentialDecoder[1]['Result']
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

	scipy.io.wavfile.write(dirname+"DemodSignal.wav", FilterDecimator['OutputSampleRate'], GFSKDemodulator[1]['Result'].astype(np.int16))
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

		report_file.write('\n\n# Demodulator performance:\n')
		report_file.write('\n')
		report_file.write(f'# Total packets: {total_packets}')
		report_file.write('\n')
		report_file.close()

	print('total packets: ', total_packets)
	print('made new directory: ', dirname)
	print('done')
