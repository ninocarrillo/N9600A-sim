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

def ModulateGauss(state):
	argv = state['argv']
	config = state['config']
	print(f'Started Shaped BPSK Modulator process')
	print(f'Reading settings for Gauss Pulse Shaping Filter')
	PulseFilter = pulse_filter.GetGaussFilterConfig(state)
	PulseFilter = pulse_filter.InitGaussFilter(PulseFilter)
	PulseFilter['SymbolMap'] = pulse_filter.GetSymbolMapConfig(state)
	NCO = nco.GetNCOConfig(config, 1, "TX NCO ")
	NCO = nco.InitNCO(NCO)

	BitStream = pulse_filter.ExpandSampleStream(state['InputData'], PulseFilter)
	ModulatingWaveform = np.convolve(PulseFilter['Taps'], BitStream)
	ModAmplitude = 64
	ModulatingWaveform = np.rint(ModAmplitude * ModulatingWaveform / max(ModulatingWaveform))

	Baseband = np.zeros(len(ModulatingWaveform))
	i = 0
	for Amplitude in ModulatingWaveform:
		NCO = nco.UpdateNCO(NCO)
		Baseband[i] = Amplitude * NCO['Sine']
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
		harmonic = 2
		harmonic_amplitude = 0/50
		harmonic_phase_shift = 4 * PulseFilter['sample rate'] // (8 * scale_factor)
		sine_table = np.zeros(PulseFilter['sample rate'] // scale_factor)
		for i in range(len(sine_table)):
			sine_table[i] = np.rint((1*(np.sin(i * 2 * np.pi / (PulseFilter['sample rate'] // scale_factor))) + (harmonic_amplitude * np.sin(harmonic * 2 * np.pi*(harmonic_phase_shift + i) / (PulseFilter['sample rate'] // scale_factor)))) * 256)

		plt.figure()
		plt.plot(sine_table)
		plt.show()

		report_file.write('\n\n#Sine Samples\n')
		report_file.write('\n')
		report_file.write(fo.GenInt16ArrayC(f'SineSamples', sine_table, 8))
		quarter_sine_table = sine_table[:(len(sine_table) // 4) + 1]
		half_sine_table = sine_table[:(len(sine_table) // 2) + 1]
		report_file.write('\n')
		report_file.write(fo.GenInt16ArrayC(f'QuarterSineSamples', quarter_sine_table, 8))
		report_file.write('\n')
		report_file.write(fo.GenInt16ArrayC(f'HalfSineSamples', half_sine_table, 8))

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
			y = int(TonePhase + random.randint(0,127))
			while y >= SineSamplesPeriod:
				y -= SineSamplesPeriod
			x = SineSamplesPeriod >> 2
			if y < x :
				y >>= 5
				sample_output[sample_index] = sine_table[y]
			elif y < (x * 2):
				y = int((SineSamplesPeriod / 2) - y) >> 5
				sample_output[sample_index] = sine_table[y]
			elif y < (x * 3):
				y = int(y - (SineSamplesPeriod / 2)) >> 5
				sample_output[sample_index] = -sine_table[y]
			else:
				y = y - (SineSamplesPeriod / 2)
				y = int((SineSamplesPeriod / 2) - y) >> 5
				sample_output[sample_index] = -sine_table[y]

		sample_output_2 = np.zeros(SineSamplesPeriod * 5)

		for sample_index in range(SineSamplesPeriod * 5):
			TonePhase += ToneFreq
			while TonePhase >= SineSamplesPeriod:
				TonePhase -= SineSamplesPeriod
			y = int(TonePhase + random.randint(0,127))
			while y >= SineSamplesPeriod:
				y -= SineSamplesPeriod
			x = SineSamplesPeriod >> 1
			if y < x :
				y >>= 5
				sample_output_2[sample_index] = sine_table[y]
			else:
				y = y - (SineSamplesPeriod / 2)
				y = int((SineSamplesPeriod / 2) - y) >> 5
				sample_output_2[sample_index] = -sine_table[y]

		sample_output_3 = np.zeros(SineSamplesPeriod * 5)
		ToneFreq = 1350
		ModPer = SineSamplesPeriod // ToneFreq
		ModIndex = 0
		ModAmplitude = 0

		for sample_index in range(SineSamplesPeriod * 5):
			TonePhase += ToneFreq
			while TonePhase >= SineSamplesPeriod:
				TonePhase -= SineSamplesPeriod
			y = int(TonePhase + random.randint(0,127))
			while y >= SineSamplesPeriod:
				y -= SineSamplesPeriod
			x = SineSamplesPeriod >> 2
			if y < x :
				y >>= 5
				sample_output_3[sample_index] = sine_table[y]
			elif y < (x * 2):
				y = int((SineSamplesPeriod / 2) - y) >> 5
				sample_output_3[sample_index] = sine_table[y]
			elif y < (x * 3):
				y = int(y - (SineSamplesPeriod / 2)) >> 5
				sample_output_3[sample_index] = -sine_table[y]
			else:
				y = y - (SineSamplesPeriod / 2)
				y = int((SineSamplesPeriod / 2) - y) >> 5
				sample_output_3[sample_index] = -sine_table[y]
			sample_output_3[sample_index] = sample_output_3[sample_index] * (ModAmplitude + 1)
			ModIndex += 1
			if ModIndex >= ModPer:
				ModIndex = 0
				ModAmplitude = random.randint(0,3)


		scipy.io.wavfile.write(dirname+"QuarterTableToneOutput.wav", SineSamplesPeriod, sample_output.astype(np.int16) * 64)

		scipy.io.wavfile.write(dirname+"HalfTableToneOutput.wav", SineSamplesPeriod, sample_output_2.astype(np.int16) * 64)


		scipy.io.wavfile.write(dirname+"ModQtrToneOutput.wav", SineSamplesPeriod, sample_output_3.astype(np.int16) * 8)
