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

def ModulateRRC(state):
	argv = state['argv']
	config = state['config']
	print(f'Started Shaped 4FSK Modulator process')
	print(f'Reading settings for RRC Pulse Shaping Filter')
	PulseFilter = pulse_filter.GetRRCFilterConfig(state)
	PulseFilter = pulse_filter.InitRRCFilter(PulseFilter)
	PulseFilter['SymbolMap'] = pulse_filter.GetSymbolMapConfig(state)
	BitStream = pulse_filter.ExpandSampleStream(state['InputData'], PulseFilter)

	waveform = np.convolve(PulseFilter['Taps'], BitStream)
	waveform_2 = np.convolve(PulseFilter['Taps'], waveform)
	PulseFilter['RC'] = np.convolve(PulseFilter['Taps'], PulseFilter['Taps'], 'same')
	plt.figure()
	plt.subplot(221)
	plt.plot(PulseFilter['Time'], PulseFilter['Taps'], 'b')
	plt.plot(PulseFilter['Time'], PulseFilter['RC'], 'r')
	plt.xticks(PulseFilter['SymbolTicks'])
	plt.xticks(color='w')
	plt.grid(True)

	plt.subplot(222)
	plt.plot(waveform, 'b')
	plt.plot(waveform_2, 'r')

	eye_data = pulse_filter.GenEyeData2(waveform_2, PulseFilter['Oversample'], PulseFilter['Oversample'] // 2)
	plt.subplot(223)
	plt.plot(eye_data)

	fft_n = len(waveform)
	x = np.linspace(0.0, fft_n * PulseFilter['TimeStep'], fft_n, endpoint = False)
	x_fft = fftfreq(fft_n, PulseFilter['TimeStep'])[:fft_n//2]
	waveform_fft = fft(waveform)
	plt.subplot(224)
	plt.plot(x_fft, 10*np.log(2.0/fft_n * np.abs(waveform_fft[0:fft_n//2])))
	plt.xlim(0,10000)
	plt.ylim(-100,10)
	plt.show()
	# GFSKDemodulator = []
	# DataSlicer = []
	# DemodulatorCount = 0
	# for DemodulatorNumber in range(4):
		# GFSKDemodulator.append({})
		# DataSlicer.append({})
		# if config.has_section(f'GFSK Demodulator {DemodulatorNumber}'):
			# print(f'Reading settings for GFSK Demodulator {DemodulatorNumber}')
			# DemodulatorCount += 1
			# GFSKDemodulator[DemodulatorNumber] = GetGFSKDemodulatorConfig(config, DemodulatorNumber)
			# GFSKDemodulator[DemodulatorNumber]['InputSampleRate'] = FilterDecimator['OutputSampleRate']
			# try:
				# DataSlicer[DemodulatorNumber]['BitRate'] = int(config[f'Data Slicer {DemodulatorNumber}']['slicer bit rate'])
			# except:
				# print(f'{sys.argv[1]} [Data Slicer {DemodulatorNumber}] \'slicer bit rate\' is missing or invalid')
				# sys.exit(-2)
			# try:
				# DataSlicer[DemodulatorNumber]['Rate'] = float(config[f'Data Slicer {DemodulatorNumber}']['slicer lock rate'])
			# except:
				# print(f'{sys.argv[1]} [Data Slicer {DemodulatorNumber}] \'slicer lock rate\' is missing or invalid')
				# sys.exit(-2)

			# DataSlicer[DemodulatorNumber]['InputSampleRate'] = FilterDecimator['OutputSampleRate']
			# GFSKDemodulator[DemodulatorNumber] = demod.InitGFSKDemod(GFSKDemodulator[DemodulatorNumber])
			# DataSlicer[DemodulatorNumber] = demod.InitDataSlicer(DataSlicer[DemodulatorNumber])

	# DifferentialDecoder = [{}]
	# AX25Decoder = [{}]
	# Descrambler = [{}]
	# for index in range(1,DemodulatorCount+1):
		# print(f'Initializing DifferentialDecoder {index}, AX25Decoder {index}, Descrambler {index}')
		# DifferentialDecoder.append({})
		# AX25Decoder.append({})
		# Descrambler.append({})
		# DifferentialDecoder[index] = demod.InitDifferentialDecoder()
		# AX25Decoder[index] = demod.InitAX25Decoder()
		# Descrambler[index]['Polynomial'] = int('0x63003',16) # G3RUH poly * differential decoding
		# Descrambler[index] = demod.InitDescrambler(Descrambler[index])

	# try:
		# samplerate, audio = scipy.io.wavfile.read(sys.argv[2])
		# # Take two bits of resolution away
		# audio = audio >> (16 - FilterDecimator['InputBitCount'])
	# except:
		# print('Unable to open wave file.')
		# sys.exit(-2)

	# print("Opened file. \r\nSample rate:", samplerate, "\r\nLength:", len(audio))


	# #generate a new directory for the reports
	# run_number = 0
	# print('trying to make a new directory')
	# while True:
		# run_number = run_number + 1
		# dirname = f'./run{run_number}/'
		# try:
			# os.mkdir(dirname)
		# except:
			# print(dirname + ' exists')
			# continue
		# break

	# total_packets = 0
	# data = np.array([])
	# filtered_signal_buffer = np.zeros(round(len(audio) / FilterDecimator['DecimationRate']) + 1)
	# demod_sig_buffer1 = np.zeros(round(len(audio) / FilterDecimator['DecimationRate']) + 1)

	# print(f'\nFiltering and decimating audio. ')
	# FilterDecimator['FilterBuffer'] = audio
	# FilterDecimator = demod.FilterDecimate(FilterDecimator)
	# print(f'Done.')

	# scipy.io.wavfile.write(dirname+"FilteredSignal.wav", FilterDecimator['OutputSampleRate'], FilterDecimator['FilterBuffer'].astype(np.int16))

	# print(f'\nDemodulating audio. ')
	# GFSKDemodulator[1]['InputBuffer'] = FilterDecimator['FilterBuffer']
	# GFSKDemodulator[1] = demod.DemodulateGFSK(GFSKDemodulator[1])
	# print(GFSKDemodulator[1]['Result'])
	# print(f'Done.')

	# print(f'\nSlicing, differential decoding, and AX25 decoding data. ')
	# loop_count = len(GFSKDemodulator[1]['Result'])
	# for index in range(loop_count):
		# DataSlicer[1]['NewSample'] = GFSKDemodulator[1]['Result'][index]
		# DataSlicer[1] = demod.ProgSliceData(DataSlicer[1])
		# for data_bit in DataSlicer[1]['Result']:
			# Descrambler[1]['NewBit'] = data_bit
			# Descrambler[1] = demod.ProgUnscramble(Descrambler[1])
			# AX25Decoder[1]['NewBit'] = Descrambler[1]['Result']
			# AX25Decoder[1] = demod.ProgDecodeAX25(AX25Decoder[1])
			# if AX25Decoder[1]['OutputTrigger'] == True:
				# AX25Decoder[1]['OutputTrigger'] = False
				# total_packets += 1
				# CRC = AX25Decoder[1]['CRC'][0]
				# decodernum = '1'
				# filename = f'Packet-{total_packets}_CRC-{format(CRC,"#06x")}_decoder-{decodernum}_Index-{index}'
				# print(f'{dirname+filename}')
				# # try:
				# # 	bin_file = open(dirname + filename + '.bin', '+wb')
				# # except:
				# # 	pass
				# # with bin_file:
				# # 	for byte in AX25Decoder[2]['Output']:
				# # 		bin_file.write(byte.astype('uint8'))
				# # 	bin_file.close()

	# scipy.io.wavfile.write(dirname+"DemodSignal.wav", FilterDecimator['OutputSampleRate'], GFSKDemodulator[1]['Result'].astype(np.int16))
	# # scipy.io.wavfile.write(dirname+"DemodSignal2.wav", FilterDecimator['OutputSampleRate'], demod_sig_buffer2.astype(np.int16))
	# #scipy.io.wavfile.write(dirname+"FilteredSignal.wav", FilterDecimator['OutputSampleRate'], filtered_signal_buffer.astype(np.int16))

	# # Generate and save report file
	# report_file_name = f'run{run_number}_report.txt'
	# try:
		# report_file = open(dirname + report_file_name, 'w+')
	# except:
		# print('Unable to create report file.')
	# with report_file:
		# report_file.write('# Command line: ')
		# for argument in sys.argv:
			# report_file.write(f'{argument} ')
		# report_file.write('\n#\n########## Begin Transcribed .ini file: ##########\n')
		# try:
			# ini_file = open(sys.argv[1])
		# except:
			# report_file.write('Unable to open .ini file.')
		# with ini_file:
			# for character in ini_file:
				# report_file.write(character)

		# report_file.write('\n\n########## End Transcribed .ini file: ##########\n')

		# report_file.write('\n\n# Demodulator performance:\n')
		# report_file.write('\n')
		# report_file.write(f'# Total packets: {total_packets}')
		# report_file.write('\n')
		# report_file.close()

	# print('total packets: ', total_packets)
	# print('made new directory: ', dirname)
	# print('done')

	return

def ModulateGauss(state):
	argv = state['argv']
	config = state['config']
	print(f'Started Shaped 4FSK Modulator process')
	print(f'Reading settings for Gauss Pulse Shaping Filter')
	PulseFilter = pulse_filter.GetGaussFilterConfig(state)
	PulseFilter = pulse_filter.InitGaussFilter(PulseFilter)
	PulseFilter['SymbolMap'] = pulse_filter.GetSymbolMapConfig(state)
	BitStream = pulse_filter.ExpandSampleStream(state['InputData'], PulseFilter)

	waveform = np.convolve(PulseFilter['Taps'], BitStream)
	ReceiveFilter = np.ones(17) / 17
	waveform_2 = np.convolve(ReceiveFilter, waveform)
	PulseFilter['RC'] = np.convolve(PulseFilter['Taps'], PulseFilter['Taps'], 'same')
	plt.figure()
	plt.subplot(221)
	plt.plot(PulseFilter['Time'], PulseFilter['Taps'], 'b')
	#plt.plot(PulseFilter['Time'], PulseFilter['RC'], 'r')
	plt.xticks(PulseFilter['SymbolTicks'])
	plt.xticks(color='w')
	plt.grid(True)

	plt.subplot(222)
	plt.plot(waveform, 'b')
	plt.plot(waveform_2, 'r')


	eye_data = pulse_filter.GenEyeData2(waveform_2, PulseFilter['Oversample'], 0)
	plt.subplot(223)
	plt.plot(eye_data)


	fft_n = len(waveform)
	x = np.linspace(0.0, fft_n * PulseFilter['TimeStep'], fft_n, endpoint = False)
	x_fft = fftfreq(fft_n, PulseFilter['TimeStep'])[:fft_n//2]
	waveform_fft = fft(waveform)
	plt.subplot(224)
	plt.plot(x_fft, 10*np.log(2.0/fft_n * np.abs(waveform_fft[0:fft_n//2])))
	plt.xlim(0,10000)
	plt.ylim(-100,10)
	plt.show()
	return
