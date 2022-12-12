import sys
import struct
import scipy.io.wavfile
import numpy as np

def PeakDetect(signal_value, detector):
	signal_value = abs(signal_value)
	if signal_value > detector['LastValue']:
		if signal_value > detector['Envelope']:
			detector['Envelope'] = signal_value
			detector['SustainCount'] = 0
			detector['DecayRate'] = np.rint(detector['Envelope'] / detector['DecayPeriod'])
			if detector['DecayRate'] < 1:
				detector['DecayRate'] = 1
	detector['LastValue'] = signal_value
	if detector['SustainCount'] >= detector['SustainPeriod']:
		detector['Envelope'] = detector['Envelope'] - detector['DecayRate']
		if detector['Envelope'] < 0:
			detector['Envelope'] = 0
			detector['SustainCount'] = 0
	detector['SustainCount'] = detector['SustainCount'] + 1
	return detector

if len(sys.argv) < 2:
	print("Not enough arguments. Usage: py -3 afsk-1200-ax25-rx.py <wav file>")
	sys.exit(-1)

try:
	samplerate, audio = scipy.io.wavfile.read(sys.argv[1])
except:
	print('Unable to open wave file.')
	sys.exit(-2)

print("Opened file. \r\nSample rate:", samplerate, "\r\nLength:", len(audio))

# Filter and decimate the audio
input_filter = np.array([77, 53, 58, 53, 39, 18, -4, -23, -38, -49, -61, -79, -106, -142, -180, -208, -215, -190, -130, -44, 52, 134, 181, 181, 137, 70, 11, -7, 39, 155, 317, 485, 604, 629, 534, 332, 67, -191, -368, -413, -315, -117, 94, 207, 124, -207, -765, -1451, -2099, -2516, -2534, -2058, -1096, 223, 1675, 2988, 3901, 4227, 3901, 2988, 1675, 223, -1096, -2058, -2534, -2516, -2099, -1451, -765, -207, 124, 207, 94, -117, -315, -413, -368, -191, 67, 332, 534, 629, 604, 485, 317, 155, 39, -7, 11, 70, 137, 181, 181, 134, 52, -44, -130, -190, -215, -208, -180, -142, -106, -79, -61, -49, -38, -23, -4, 18, 39, 53, 58, 53, 77])
input_filter_buffer = np.zeros(len(input_filter))


#create a dictionary for the input signal peak detector
InputPeakDetector = {'AttackPeriod':1,'SustainPeriod':7200,'DecayPeriod':7200, 'AttackCount':0, 'SustainCount':0, 'DecayCount':0, 'LastValue':0, 'Envelope':0, 'DecayRate':1}

index1 = 0
index2 = 0
index3 = 0
envelope_index = 0
decimation = 1
input_filter_gain = 0
envelope = np.zeros(round(len(audio) / decimation))
filtered_signal_buffer = np.zeros(round(len(audio) / decimation))
for sample in audio:
	index1 = index1 + 1
	index2 = index2 + 1
	if index2 > len(audio) / 100:
		index2 = 0
		index3 = index3 + 1
		print(index3, InputPeakDetector['Envelope'])
	input_filter_buffer = input_filter_buffer[1:]
	input_filter_buffer = np.append(input_filter_buffer, np.array([sample]))
	if index1 == decimation:
		index1 = 0
		filtered_signal = np.rint(np.convolve(input_filter_buffer, input_filter, 'valid') / pow(2, (16 + input_filter_gain)))
		InputPeakDetector = PeakDetect(filtered_signal[0], InputPeakDetector)
		filtered_signal_buffer[envelope_index] = filtered_signal[0]
		envelope[envelope_index] = InputPeakDetector['Envelope']
		envelope_index = envelope_index + 1
		if InputPeakDetector['Envelope'] > 24576:
			input_filter_gain = input_filter_gain + 1
			print('going down', input_filter_gain)
			InputPeakDetector['Envelope'] = InputPeakDetector['Envelope'] / 2
			if input_filter_gain > 16:
				input_filter_gain = 16
		if InputPeakDetector['Envelope'] < 8192:
			input_filter_gain = input_filter_gain - 1
			InputPeakDetector['Envelope'] = InputPeakDetector['Envelope'] * 2
			print('going up', input_filter_gain)
			if input_filter_gain < -16:
				input_filter_gain = -16
		#print(filtered_signal[0], InputPeakDetector['Envelope'])
scipy.io.wavfile.write("PeakDetect.wav", round(samplerate / decimation), envelope.astype(np.int16))
scipy.io.wavfile.write("FilteredSignal.wav", round(samplerate / decimation), filtered_signal_buffer.astype(np.int16))
print('done')
