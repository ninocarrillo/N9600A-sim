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
		afsk_demodulator['OutputFilter'] = strings.StringToIntArray(config[f'AFSK Demodulator {num}']['output filter taps'])
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

	try:
		afsk_demodulator['sqrt shift'] = int(config[f'AFSK Demodulator {num}']['sqrt shift'])
	except:
		print(f'{sys.argv[1]} [AFSK Demodulator {num}] \'sqrt shift\' is missing or invalid')
		sys.exit(-2)


	try:
		afsk_demodulator['OutputFilterDecimationRate'] = int(config[f'AFSK Demodulator {num}']['output filter decimation'])
	except:
		print(f'{sys.argv[1]} [AFSK Demodulator {num}] \'output filter decimation\' is missing or invalid')
		sys.exit(-2)
	try:
		afsk_demodulator['CorrelatorDecimationRate'] = int(config[f'AFSK Demodulator {num}']['correlator decimation'])
	except:
		print(f'{sys.argv[1]} [AFSK Demodulator {num}] \'correlator decimation\' is missing or invalid')
		sys.exit(-2)
	return afsk_demodulator


def FullProcess(state):
	argv = state['argv']
	config = state['config']

	print(f'Started AFSK process')
	print(f'Reading settings for Filter Decimator')
	FilterDecimator = input_filter.GetInputFilterConfig(state)
	FilterDecimator = demod.InitFilterDecimator(FilterDecimator)

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
				print(f'{sys.argv[1]} [Data Slicer {DemodulatorNumber}] \'slicer bit rate\' is missing or invalid')
				sys.exit(-2)
			try:
				DataSlicer[DemodulatorNumber]['Rate'] = float(config['Data Slicer 1']['slicer lock rate'])
			except:
				print(f'{sys.argv[1]} [Data Slicer {DemodulatorNumber}] \'slicer lock rate\' is missing or invalid')
				sys.exit(-2)

			AFSKDemodulator[DemodulatorNumber] = demod.InitAFSKDemod(AFSKDemodulator[DemodulatorNumber])
			DataSlicer[DemodulatorNumber]['InputSampleRate'] = AFSKDemodulator[DemodulatorNumber]['OutputSampleRate']
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
	print('Trying to make a new directory.')
	while True:
		run_number = run_number + 1
		dirname = f'./run{run_number}/'
		try:
			os.mkdir(dirname)
		except:
			print(dirname + ' exists')
			continue
		print('Made ' + dirname)
		break


	total_packets = 0
	duplicate_packets = 0

	data = np.array([])
	filtered_signal_buffer = np.zeros(round(len(audio) / FilterDecimator['DecimationRate']) + 1)

	print(f'\nFiltering and decimating audio. ')
	FilterDecimator['FilterBuffer'] = audio
	FilterDecimator = demod.FilterDecimate(FilterDecimator)
	print(f'Done.')

	scipy.io.wavfile.write(dirname+"FilteredSignal.wav", FilterDecimator['OutputSampleRate'], FilterDecimator['FilterBuffer'].astype(np.int16))
	print(f'Wrote file {dirname}FilteredSignal.wav at {FilterDecimator["OutputSampleRate"]} samples per second.')

	print(f'\nDemodulating audio.')
	loop_count = 0
	for index in range(1, DemodulatorCount + 1):
		AFSKDemodulator[index]['CorrelatorBuffer'] = FilterDecimator['FilterBuffer']
		AFSKDemodulator[index] = demod.DemodulateAFSK(AFSKDemodulator[index])
		if loop_count == 0:
			try:
				loop_count = len(AFSKDemodulator[index]['Result'])
			except:
				pass
		else:
			try:
				if loop_count > len(AFSKDemodulator[index]['Result']):
					loop_count = len(AFSKDemodulator[index]['Result'])
			except:
				pass

	print(f'Done.')

	print(f'\nSlicing, differential decoding, and AX25 decoding data. ')

	for index in range(loop_count):
		for demod_index in range(1, DemodulatorCount + 1):
			if AFSKDemodulator[demod_index]['Enabled'] == True:
				DataSlicer[demod_index]['NewSample'] = AFSKDemodulator[demod_index]['Result'][index]
				DataSlicer[demod_index] = demod.ProgSliceData(DataSlicer[demod_index])
				for data_bit in DataSlicer[demod_index]['Result']:
					DifferentialDecoder[demod_index]['NewBit'] = data_bit
					DifferentialDecoder[demod_index] = demod.ProgDifferentialDecode(DifferentialDecoder[demod_index])
					AX25Decoder[demod_index]['NewBit'] = DifferentialDecoder[demod_index]['Result']
					AX25Decoder[demod_index] = demod.ProgDecodeAX25(AX25Decoder[demod_index])
					if AX25Decoder[demod_index]['OutputTrigger'] == True:
						AX25Decoder[demod_index]['OutputTrigger'] = False
						AX25Decoder[demod_index]['LastPacketSampleIndex'] = DataSlicer[demod_index]['SampleIndex']
						# Check for unioqueness
						unique = True
						others_old = True
						for index2 in range(1, DemodulatorCount + 1):
							if index2 != demod_index:
								try:
									if abs(AX25Decoder[index2]['LastPacketSampleIndex'] - AX25Decoder[demod_index]['LastPacketSampleIndex']) < 100:
										others_old = False
								except:
									pass
						# others_old = False
						if others_old == False:
							# there is another new packet. Need to check it.
							for index2 in range(1, DemodulatorCount + 1):
								if index2 != demod_index:
									if AX25Decoder[demod_index]['CRC'][0] == AX25Decoder[index2]['CRC'][0]:
										# print(AX25Decoder[demod_index]['CRC'][0])
										# print(AX25Decoder[index2]['CRC'][0])
										AX25Decoder[index2]['UniquePackets'] -= 1
										unique = False
						if unique == True:
						# if ((AX25Decoder[2]['CRCAge'] > 10) and (AX25Decoder[3]['CRCAge'] > 10)) or ((AX25Decoder[1]['CRC'][0] != AX25Decoder[2]['CRC'][0]) and (AX25Decoder[1]['CRC'] != AX25Decoder[3]['CRC'])):
							total_packets += 1
							AX25Decoder[demod_index]['UniquePackets'] += 1
							CRC = AX25Decoder[demod_index]['CRC'][0]
							filename = f'Packet-{total_packets}_CRC-{format(CRC,"#06x")}_decoder-{demod_index}'
							print(f'{dirname+filename}')
							# for byte in AX25Decoder[demod_index]['Output']:
								# print(hex(byte), end = ' ')
							# print(' ')
							# try:
								# bin_file = open(dirname + filename + '.bin', '+wb')
							# except:
								# pass
							# with bin_file:
								# for byte in AX25Decoder[demod_index]['Output']:
									# bin_file.write(byte.astype('uint8'))
								# bin_file.close()


	scipy.io.wavfile.write(dirname+"DemodSignal1.wav", FilterDecimator['OutputSampleRate'], AFSKDemodulator[1]['Result'].astype(np.int16))
	#scipy.io.wavfile.write(dirname+"DemodSignal2.wav", FilterDecimator['OutputSampleRate'], AFSKDemodulator[2]['Result'].astype(np.int16))
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


		report_file.write(f'\nFilter Decimator Output Sample Rate: {FilterDecimator["OutputSampleRate"]}.')
		print(f'\nFilter Decimator Output Sample Rate: {FilterDecimator["OutputSampleRate"]}.')

		tstep = 1.0 / AFSKDemodulator[1]['InputSampleRate']

		time = np.arange(0, tstep * AFSKDemodulator[1]['CorrelatorTapCount'], tstep)

		index = 0
		while index < DemodulatorCount:
			index += 1

			report_file.write(f'\nDemodulator {index} Output Sample Rate: {AFSKDemodulator[index]["OutputSampleRate"]}.')
			print(f'\nDemodulator {index} Output Sample Rate: {AFSKDemodulator[index]["OutputSampleRate"]}.')
			report_file.write(fo.GenInt16ArrayC(f'MarkCorCos{index}', AFSKDemodulator[index]['MarkCOS'], 8))
			report_file.write(fo.GenInt16ArrayC(f'MarkCorSin{index}', AFSKDemodulator[index]['MarkSIN'], 8))
			report_file.write(fo.GenInt16ArrayC(f'SpaceCorCos{index}', AFSKDemodulator[index]['SpaceCOS'], 8))
			report_file.write(fo.GenInt16ArrayC(f'SpaceCorSin{index}', AFSKDemodulator[index]['SpaceSIN'], 8))
			report_file.write('\n')

			max_mark_cos_input = np.rint(24576 * (np.cos(2 * AFSKDemodulator[index]['MarkFreq'] * np.pi * time)))
			max_mark_sin_input = np.rint(24576 * (np.sin(2 * AFSKDemodulator[index]['MarkFreq'] * np.pi * time)))
			max_space_cos_input = np.rint(24576 * (np.cos(2 * AFSKDemodulator[index]['SpaceFreq'] * np.pi * time)))
			max_space_sin_input = np.rint(24576 * (np.sin(2 * AFSKDemodulator[index]['SpaceFreq'] * np.pi * time)))

			mark_sin = np.convolve(max_mark_sin_input, AFSKDemodulator[index]["MarkSIN"], "valid") / pow(2, (16 + AFSKDemodulator[index]['CorrelatorShift']))
			mark_cos = np.convolve(max_mark_sin_input, AFSKDemodulator[index]["MarkCOS"], "valid") / pow(2, (16 + AFSKDemodulator[index]['CorrelatorShift']))
			mark_square_sum = mark_sin**2 + mark_cos**2

			space_sin = np.convolve(max_space_sin_input, AFSKDemodulator[index]["SpaceSIN"], "valid") / pow(2, (16 + AFSKDemodulator[index]['CorrelatorShift']))
			space_cos = np.convolve(max_space_sin_input, AFSKDemodulator[index]["SpaceCOS"], "valid") / pow(2, (16 + AFSKDemodulator[index]['CorrelatorShift']))
			space_square_sum = space_sin**2 + space_cos**2

			report_file.write(f'\nMaximum correlation values Demodulator {index}:')
			report_file.write('\n')
			report_file.write(f'\nMark Sin Correlator: {int(np.rint(mark_sin[0]))}')
			report_file.write(f'\nMark Cos Correlator: {int(np.rint(mark_cos[0]))}')
			report_file.write(f'\nMark Square Sum Correlator: {int(np.rint(mark_square_sum[0]))}')
			report_file.write(f'\nMark Square Sum / Square Scale: {int(np.rint(mark_square_sum[0] / AFSKDemodulator[index]["SquareScale"]))}')
			# report_file.write(f'\n Log2(Mark Square Sum / Sqrt Table Size): {np.log2(mark_square_sum_1[0] // 2**AFSKDemodulator[1]["SqrtBitCount"]):.2f}')

			report_file.write('\n')
			report_file.write(f'\nSpace Sin Correlator: {int(np.rint(space_sin[0]))}')
			report_file.write(f'\nSpace Cos Correlator: {int(np.rint(space_cos[0]))}')
			report_file.write(f'\nSpace Square Sum Correlator: {int(np.rint(space_square_sum[0]))}')
			report_file.write(f'\nMark Square Sum / Square Scale: {int(np.rint(space_square_sum[0] / AFSKDemodulator[index]["SquareScale"]))}')
			# report_file.write(f'\n Log2(Space Square Sum / Sqrt Table Size): {np.log2(space_square_sum_1[0] // 2**AFSKDemodulator[1]["SqrtBitCount"]):.2f}')


			report_file.write(f'\nSquare Scale {index}: {AFSKDemodulator[index]["SquareScale"]}')

			report_file.write(fo.GenInt16ArrayC(f'\nSquareRoot{AFSKDemodulator[index]["SqrtBitCount"]}', AFSKDemodulator[1]['SqrtTable'], 8))


			report_file.write(f'\nDemodulator {index} total packets: {AX25Decoder[index]["PacketCount"]}')

			# print(f'Decoder {index} unique packets: ', AX25Decoder[index]['UniquePackets'])
			print(f'Decoder {index} total packets: ', AX25Decoder[index]['PacketCount'])

		report_file.write('\n\n# Demodulator performance:\n')
		report_file.write('\n')
		report_file.write(f'# Total packets: {total_packets}')
		report_file.write('\n')

	print(f'\nTotal packets: {total_packets}.')
	print(f'\nDone.')
	return
