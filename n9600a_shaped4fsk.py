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
	plt.suptitle(f"RRC 4FSK Rolloff Rate:{PulseFilter['rolloff rate']}, Span:{PulseFilter['symbol span']}, Sample Rate:{PulseFilter['sample rate']}")
	plt.subplot(221)
	plt.plot(PulseFilter['Time'], PulseFilter['Taps'], 'b')
	plt.plot(PulseFilter['Time'], PulseFilter['RC'], 'r')
	plt.xticks(PulseFilter['SymbolTicks'])
	plt.xticks(color='w')
	#plt.xlabel("Symbol Intervals")
	plt.title("Impulse Response")
	plt.grid(True)

	plt.subplot(222)
	plt.plot(waveform, 'b')
	plt.plot(waveform_2, 'r')
	plt.title("Modulation Waveform")

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

	scipy.io.wavfile.write(dirname+"ModSignal.wav", PulseFilter['sample rate'], waveform.astype(np.int16))
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

		report_file.write('\n\n# Filter\n')
		report_file.write('\n')
		report_file.write(f"# Filter Taps: {PulseFilter['Taps']}")
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
	return
