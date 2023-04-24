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
	#PulseFilter['RC'] = np.convolve(PulseFilter['Taps'], PulseFilter['Taps'], 'same')
	plt.figure()
	plt.suptitle(f"RRC 4FSK Rolloff Rate:{PulseFilter['rolloff rate']}, Span:{PulseFilter['symbol span']}, Sample Rate:{PulseFilter['sample rate']}")
	plt.subplot(221)
	plt.plot(PulseFilter['Time'], PulseFilter['Taps'], 'b')
	plt.plot(PulseFilter['Time'], PulseFilter['RC'], 'r')
	plt.xticks(PulseFilter['SymbolTicks'])
	plt.xticks(color='w')
	#plt.xlabel("Symbol Intervals")
	plt.title("Impulse Response")
	plt.legend(["RRC", "RC"])
	plt.grid(True)

	plt.subplot(222)
	plt.plot(waveform, 'b')
	plt.plot(waveform_2, 'r')
	plt.title("Modulation Waveform")
	plt.legend(["Post Transmit Filter", "Post Receive Filter"])

	eye_data = pulse_filter.GenEyeData2(waveform_2, PulseFilter['Oversample'], PulseFilter['Oversample'] // 2)
	plt.subplot(223)
	plt.plot(eye_data)
	plt.title("Eye Diagram")

	fft_n = len(waveform)
	x = np.linspace(0.0, fft_n * PulseFilter['TimeStep'], fft_n, endpoint = False)
	x_fft = fftfreq(fft_n, PulseFilter['TimeStep'])[:fft_n//2]
	waveform_fft = fft(waveform)
	fft_max = max(abs(waveform_fft))
	waveform_fft = waveform_fft / fft_max
	plt.subplot(224)
	plt.plot(x_fft, 10*np.log(np.abs(waveform_fft[0:fft_n//2])))
	plt.xlim(0,10000)
	plt.ylim(-100,10)
	plt.title("Frequency Spectrum")
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

	waveform = waveform / max(waveform)
	waveform = waveform * 32767

	waveform_2 = waveform_2 / max(waveform_2)
	waveform_2 = waveform_2 * 32767

	scipy.io.wavfile.write(dirname+"ModSignal.wav", PulseFilter['sample rate'], waveform.astype(np.int16))
	scipy.io.wavfile.write(dirname+"DemodSignal.wav", PulseFilter['sample rate'], waveform_2.astype(np.int16))
	# # scipy.io.wavfile.write(dirname+"DemodSignal2.wav", FilterDecimator['OutputSampleRate'], demod_sig_buffer2.astype(np.int16))
	# #scipy.io.wavfile.write(dirname+"FilteredSignal.wav", FilterDecimator['OutputSampleRate'], filtered_signal_buffer.astype(np.int16))

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

		report_file.write('\n\n# RRC Pulse Filter\n')
		report_file.write('\n')
		report_file.write(fo.GenInt16ArrayC(f'RRCFilter', PulseFilter['Taps'] * 65536, 8))
		report_file.write('\n')
		report_file.close()



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
	waveform_fft = fft(waveform)
	fft_max = max(abs(waveform_fft))
	waveform_fft = waveform_fft / fft_max
	plt.subplot(224)
	plt.plot(x_fft, 10*np.log(np.abs(waveform_fft[0:fft_n//2])))
	plt.xlim(0,10000)
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

	waveform = waveform / max(waveform)
	waveform = waveform * 32767
	scipy.io.wavfile.write(dirname+"ModSignal.wav", PulseFilter['sample rate'], waveform.astype(np.int16))
	return

def GaussFilterGen(state):
	argv = state['argv']
	config = state['config']
	print(f'Started GaussFilterGen')
	print(f'Reading settings for Gauss Pulse Shaping Filter')
	PulseFilter = pulse_filter.GetGaussFilterConfig(state)
	PulseFilter = pulse_filter.InitGaussFilter(PulseFilter)

	# Adjust gain of filter:
	FilterSum = np.sum(PulseFilter['Taps'])
	FilterAdj = 1 / FilterSum
	PulseFilter['TapsTrimmed'] = np.trim_zeros(PulseFilter['Taps'], trim='fb')
	PulseFilter['Taps'] = FilterAdj * PulseFilter['Taps']
	samples_per_symbol = PulseFilter['sample rate'] // PulseFilter['symbol rate']

	# i = 0
	# for tap in PulseFilter['Taps']:
	# 	PulseFilter['Taps'][i] = np.rint(tap * FilterAdj)

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

		report_file.write('\n\n#Gauss Filter Taps\n')
		report_file.write('\n')
		report_file.write(fo.GenInt16ArrayC(f'GaussFilter', np.rint(PulseFilter['TapsTrimmed'] * 65536), 8))
		report_file.write('\n\n')
		report_file.write(f'Filter Tap Sum: {np.sum(PulseFilter["TapsTrimmed"])}')

		plt.figure()
		plt.plot(PulseFilter['Time'], PulseFilter['Taps'], 'b')
		#plt.plot(PulseFilter['Time'], PulseFilter['RC'], 'r')
		plt.xticks(PulseFilter['SymbolTicks'])
		plt.xticks(color='w')
		plt.grid(True)
		plt.show()

		# generate the symbol patterns
		symbol_count = np.power(2, PulseFilter['symbol span'])
		for i0 in range(symbol_count):
			#i0 is the symbol to be factored
			y = np.zeros(samples_per_symbol * PulseFilter['symbol span'])
			factor_me = i0
			pulse_pattern = 0
			for i2 in range(PulseFilter['symbol span']):
				factor = factor_me % 2
				factor_me = factor_me // 2
				if factor == 1:
					pulse_pattern += 2**(PulseFilter['symbol span'] - i2 - 1)
					level = 1600
				else:
					level = 1800
				for i3 in range(samples_per_symbol):
					y[(i2 * samples_per_symbol) + i3] = level
			print(i0)
			#plt.figure()
			#plt.plot(y)
			#plt.show()


			z = np.rint(np.convolve(y, PulseFilter['Taps'], 'full'))
			# trim the invalid results:
			z = z[len(PulseFilter['Taps'])//2:]
			z = z[:samples_per_symbol*PulseFilter['symbol span']]

			# select the center symbol length
			x_offset = ((PulseFilter['symbol span'] * samples_per_symbol) // 2) - (samples_per_symbol // 2)
			xz = np.arange(x_offset, x_offset + samples_per_symbol)
			z = z[x_offset:x_offset+samples_per_symbol]
			plt.figure()
			plt.title(f'Pulse Pattern {pulse_pattern}')
			plt.plot(y)
			plt.plot(xz,z)
			plt.ylim(1500,1900)
			plt.show()


			report_file.write('\n')
			report_file.write(fo.GenInt16ArrayC(f'PulseTaps{pulse_pattern}', z, 8))
			report_file.write('\n\n')


		#generate a sine table:
		scale_factor = 32
		sine_table = np.zeros(PulseFilter['sample rate'] // scale_factor)
		for i in range(len(sine_table)):
			sine_table[i] = np.rint(np.sin(i * 2 * np.pi / (PulseFilter['sample rate'] // scale_factor)) * 256)

		report_file.write('\n\n#Sine Samples\n')
		report_file.write('\n')
		report_file.write(fo.GenInt16ArrayC(f'SineSamples', sine_table, 8))
		sine_table = sine_table[:(len(sine_table) // 4) + 1]
		report_file.write('\n\n#Sine Samples\n')
		report_file.write('\n')
		report_file.write(fo.GenInt16ArrayC(f'SineSamples', sine_table, 8))
		
		report_file.close()


		#generate a sine wave and save to wav file:
		
		SineSamplesPeriod = int(86400)
		ToneFreq = 1700
		TonePhase = 0
		sample_output = np.zeros(SineSamplesPeriod * 5)
		
		for sample_index in range(SineSamplesPeriod * 5):
			TonePhase += ToneFreq
			while TonePhase >= SineSamplesPeriod:
				TonePhase -= SineSamplesPeriod
			y = int(TonePhase + random.randint(0,63))
			while y >= SineSamplesPeriod:
				y -= SineSamplesPeriod
			x = SineSamplesPeriod >> 2
			
			if y < x :
				y >>= 5
				sample_output[sample_index] = sine_table[y] + 104
			elif y < (x * 2):
				y = int((SineSamplesPeriod / 2) - y) >> 5
				sample_output[sample_index] = sine_table[y] + 104
			elif y < (x * 3):
				y = int(y - (SineSamplesPeriod / 2)) >> 5
				sample_output[sample_index] = -sine_table[y] + 104
			else:
				y = y - (SineSamplesPeriod / 2)
				y = int((SineSamplesPeriod / 2) - y) >> 5
				sample_output[sample_index] = -sine_table[y] + 104
			
		plt.figure()
		plt.plot(sample_output)
		plt.show()
		
		scipy.io.wavfile.write(dirname+"ToneOutput.wav", SineSamplesPeriod, sample_output.astype(np.int16) * 64)