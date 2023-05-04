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

		report_file.close()


	return
