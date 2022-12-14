import sys
import struct
import scipy.io.wavfile
import numpy as np
import os


def SliceData(samples, rate, oversample):
	LockRate = rate
	length = len(samples) / oversample;
	length = np.rint(length).astype(np.int32)
	sliced_data = np.zeros(length + 10000, np.int16)
	bit_index = 0
	LastSample = 0
	PLLClock = 0.0
	PLLStep = 1000000.0
	PLLPeriod = oversample * PLLStep
	for sample in samples:
		PLLClock += PLLStep
		if (PLLClock > ((PLLPeriod / 2) - 1)):
			PLLClock -= PLLPeriod
			if (sample > 0):
				sliced_data[bit_index] = 1
			bit_index += 1
		if LastSample > 0:
			if sample <= 0:
				# Zero Crossing
				PLLClock *= LockRate
		else:
			if sample > 0:
				#Zero Crossing
				PLLClock *= LockRate
		LastSample = sample
	return sliced_data

def DifferentialDecode(bitstream):
	state = 0
	output_data = np.zeros(len(bitstream))
	output_index = 0
	for bit in bitstream:
		if bit == state:
			output_data[output_index] = 1
		output_index += 1
		state = bit
	return output_data

def CheckCRC(packet, byte_count):
	packet_val = packet[byte_count - 1] * 256
	packet_val += packet[byte_count - 2]
	# print(hex(packet_val))
	packet = packet[0:byte_count - 2] # slicing exludes the end index
	fcsval = np.uint16(0xFFFF)
	CRC_poly = np.uint16(0x8408)
	one = np.uint16(1)
	for byte in packet:
		for i in range(8):
			fcsbit = np.bitwise_and(fcsval, one)
			fcsval = np.right_shift(fcsval, 1)
			if np.bitwise_xor(fcsbit, np.bitwise_and(byte,one)) != 0:
				fcsval = np.bitwise_xor(fcsval, CRC_poly)
			byte = np.right_shift(byte, 1)
	fcs_val = np.bitwise_and(np.bitwise_not(fcsval), 0xFFFF)
	if packet_val == fcs_val:
		return [fcs_val, 1]
	else:
		return [fcs_val, 0]

def DecodeAX25(bitstream, verbose):
	output_buffer = np.zeros(1024, np.uint16)
	packet_count = 0
	ones = 0
	working_byte = np.uint16(0)
	bit_index = 0
	byte_count = 0
	stream_index = 0
	for bit in bitstream:
		stream_index += 1
		if bit == 1:
			working_byte = np.bitwise_or(working_byte, 128)
			ones += 1
			bit_index += 1
			if ones > 6:
				# abort frame for invalid bit sequence
				ones = 0
				bit_index = 0
				byte_count = 0
				# print(" Frame Abort")
			elif bit_index > 7:
				# 8 valid bits received, record byte
				bit_index = 0
				# if working_byte > 31 and working_byte < 128:
				#	print(chr(working_byte), end=' ')
				# else:
				#	print(hex(working_byte), end=' ')

				output_buffer[byte_count] = working_byte
				byte_count += 1
			else:
				working_byte = np.right_shift(working_byte, 1)
		else:
			if ones < 5:
				working_byte = np.bitwise_and(working_byte, 127)
				bit_index += 1
				if bit_index > 7:
					bit_index = 0
					# if working_byte > 31 and working_byte < 128:
					#	print(chr(working_byte), end=' ')
					# else:
					#	print(hex(working_byte), end=' ')
					try:
						output_buffer[byte_count] = working_byte
					except:
						pass
					else:
						byte_count += 1
				else:
					working_byte = np.right_shift(working_byte, 1)
			elif ones == 5:
				ones == 5
				# ignore stuffed zero
			elif ones == 6:
				# Frame complete
				if byte_count > 18:
					CRC = CheckCRC(output_buffer, byte_count)
					if CRC[1] == 1:
						packet_count += 1
						if verbose == 1:
							print(stream_index, hex(CRC[0]), packet_count)
							for i in range(byte_count - 2):
								character = output_buffer[i]
								if (character > 31) and (character < 128):
									print(''.join(chr(character)), end='')
							print("\r\n", end='')
				byte_count = 0
				bit_index = 0
			else:
				# Invalid frame
				byte_count = 0
				bit_index = 0
			ones = 0
	return packet_count


def HighLowDetect(signal_value, detector):
	#if (signal_value > detector['LastValue']):
		#signal is increasing, check for new High
	if signal_value > detector['High']:
		#detector['High'] = signal_value
		#detector['HighAttackRate'] = signal_value / detector['AttackPeriod']
		#if detector['HighAttackRate'] == 0:
			##detector['HighAttackRate'] = 1
		detector['High'] = detector['High'] + detector['AttackRate']
		if detector['High'] > signal_value:
			detector['High'] = signal_value
		detector['HighSustainCount'] = 0
		#detector['HighDecayRate'] = np.rint(detector['High'] / detector['DecayPeriod'])
		#if detector['HighDecayRate'] == 0:
			#detector['HighDecayRate'] = 1
	#if (signal_value < detector['LastValue']):
		#signal is decreasing, check for new Low
	if signal_value < detector['Low']:
		#detector['Low'] = signal_value
		#detector['LowAttackRate'] = signal_value / detector['AttackPeriod']
		#if detector['LowAttackRate'] == 0:
			#detector['LowAttackRate'] = -1
		detector['Low'] = detector['Low'] - detector['AttackRate']
		if detector['Low'] < signal_value:
			detector['Low'] = signal_value
		detector['LowSustainCount'] = 0
		#detector['LowDecayRate'] = np.rint(detector['Low'] / detector['DecayPeriod'])
		#if detector['LowDecayRate'] == 0:
			#detector['LowDecayRate'] = -1
	detector['LastValue'] = signal_value
	if detector['HighSustainCount'] >= detector['SustainPeriod']:
		detector['High'] = detector['High'] - detector['DecayRate']
		if detector['High'] <= 0:
			detector['High'] = 1
			detector['HighSustainCount'] = 0
	if detector['LowSustainCount'] >= detector['SustainPeriod']:
		detector['Low'] = detector['Low'] + detector['DecayRate']
		if detector['Low'] >= 0:
			detector['Low'] = -1
			detector['LowSustainCount'] = 0
	detector['HighSustainCount'] = detector['HighSustainCount'] + 1
	detector['LowSustainCount'] = detector['LowSustainCount'] + 1
	return detector

def PeakDetect(signal_value, detector):
	signal_value = abs(signal_value)
	if signal_value > detector['Envelope']:
		detector['Envelope'] = detector['Envelope'] + detector['AttackRate']
		if detector['Envelope'] > signal_value:
			detector['Envelope'] = signal_value
		detector['SustainCount'] = 0
	if detector['SustainCount'] >= detector['SustainPeriod']:
		detector['Envelope'] = detector['Envelope'] - detector['DecayRate']
		if detector['Envelope'] < 0:
			detector['Envelope'] = 0
			detector['SustainCount'] = 0
	detector['SustainCount'] = detector['SustainCount'] + 1
	return detector

def FilterDecimate(filter):
	filter['FilterBuffer'] = filter['FilterBuffer'][1:]
	filter['FilterBuffer'] = np.append(filter['FilterBuffer'], np.array([filter['NewSample']]))
	output_buffer = np.rint(np.convolve(filter['FilterBuffer'], filter['Filter'], 'valid') / pow(2, (16 + filter['FilterShift'])))
	filter['DataBuffer'] = np.array([])
	for data in output_buffer:
		filter['DecimationCounter'] = filter['DecimationCounter'] + 1
		if filter['DecimationCounter'] >= filter['DecimationRate']:
			filter['DecimationCounter'] = 0
			filter['DataBuffer'] = np.append(filter['DataBuffer'], np.array([data]))
			filter['PeakDetector'] = PeakDetect(data, filter['PeakDetector'])
			if filter['PeakDetector']['Envelope'] > 24576:
				filter['FilterShift'] = filter['FilterShift'] + 1
				filter['PeakDetector']['Envelope'] = filter['PeakDetector']['Envelope'] / 2
				if filter['FilterShift'] > 16:
					filter['FilterShift'] = 16
			if filter['PeakDetector']['Envelope'] < 8192:
				filter['FilterShift'] = filter['FilterShift'] - 1
				filter['PeakDetector']['Envelope'] = filter['PeakDetector']['Envelope'] * 2
				if filter['FilterShift'] < -16:
					filter['FilterShift'] = -16
	return filter

def DemodulateAFSK(demodulator):
	demodulator['CorrelatorBuffer'] = demodulator['CorrelatorBuffer'][1:]
	demodulator['CorrelatorBuffer'] = np.append(demodulator['CorrelatorBuffer'], np.array([demodulator['NewSample']]))

	mark_cos_sig = np.rint((demodulator['CorrelatorShift'] * np.convolve(demodulator['CorrelatorBuffer'], demodulator['MarkCOS'], 'valid')) / pow(2, (16 + demodulator['CorrelatorShift'])))
	mark_sin_sig = np.rint((demodulator['CorrelatorShift'] * np.convolve(demodulator['CorrelatorBuffer'], demodulator['MarkSIN'], 'valid')) / pow(2, (16 + demodulator['CorrelatorShift'])))

	mark_sig = np.add(np.square(mark_cos_sig), np.square(mark_sin_sig))
	mark_sig = np.rint(mark_sig / demodulator['SquareScale'])
	mark_sig = np.clip(mark_sig, 0, demodulator['SquareClip'])

	mark_sig = np.rint(demodulator['SquareOutputScale'] * np.sqrt(demodulator['SquareCoef'] * mark_sig))

	space_cos_sig = np.rint((demodulator['CorrelatorShift'] * np.convolve(demodulator['CorrelatorBuffer'], demodulator['SpaceCOS'], 'valid')) / pow(2, (16 + demodulator['CorrelatorShift'])))
	space_sin_sig = np.rint((demodulator['CorrelatorShift'] * np.convolve(demodulator['CorrelatorBuffer'], demodulator['SpaceSIN'], 'valid')) / pow(2, (16 + demodulator['CorrelatorShift'])))

	space_sig = np.add(np.square(space_cos_sig), np.square(space_sin_sig))
	space_sig = np.rint(space_sig / demodulator['SquareScale'])
	space_sig = np.clip(space_sig, 0, demodulator['SquareClip'])

	space_sig = np.rint(demodulator['SquareOutputScale']* np.sqrt(demodulator['SquareCoef'] * space_sig))

	demodulator['OutputFilterBuffer'] = demodulator['OutputFilterBuffer'][1:]
	demodulator['OutputFilterBuffer'] = np.append(demodulator['OutputFilterBuffer'], np.array([mark_sig - space_sig]))

	demodulator['Result'] = np.rint(np.convolve(demodulator['OutputFilterBuffer'], demodulator['OutputFilter'], 'valid') / pow(2, (16 + demodulator['OutputFilterShift'])))

	return demodulator

def ProgSliceData(slicer):
	slicer['Result'] = np.array([])
	slicer['PLLClock'] += slicer['PLLStep']
	if slicer['PLLClock'] > ((slicer['PLLPeriod'] / 2.0) - 1.0):
		slicer['PLLClock'] -= slicer['PLLPeriod']
		if slicer['NewSample'] > 0.0:
			slicer['Result'] = np.array([1])
		else:
			slicer['Result'] = np.array([0])
		if slicer['LastSample'] > 0.0:
			if slicer['NewSample'] <= 0.0:
				# Zero Crossing
				slicer['PLLClock'] *= slicer['Rate']
		else:
			if slicer['NewSample'] > 0.0:
				# Zero Crossing
				slicer['PLLClock'] *= slicer['Rate']
		slicer['LastSample'] = slicer['NewSample']
	return slicer

def ProgDifferentialDecode(decoder):
	if decoder['NewBit'] == decoder['LastBit']:
		decoder['Result'] = 1
	else:
		decoder['Result'] = 0
	decoder['LastBit'] = decoder['NewBit']
	return decoder

def ProgDecodeAX25(bitstream, verbose):
	output_buffer = np.zeros(1024, np.uint16)
	packet_count = 0
	ones = 0
	working_byte = np.uint16(0)
	bit_index = 0
	byte_count = 0
	stream_index = 0
	for bit in bitstream:
		stream_index += 1
		if bit == 1:
			working_byte = np.bitwise_or(working_byte, 128)
			ones += 1
			bit_index += 1
			if ones > 6:
				# abort frame for invalid bit sequence
				ones = 0
				bit_index = 0
				byte_count = 0
				# print(" Frame Abort")
			elif bit_index > 7:
				# 8 valid bits received, record byte
				bit_index = 0
				# if working_byte > 31 and working_byte < 128:
				#	print(chr(working_byte), end=' ')
				# else:
				#	print(hex(working_byte), end=' ')

				output_buffer[byte_count] = working_byte
				byte_count += 1
			else:
				working_byte = np.right_shift(working_byte, 1)
		else:
			if ones < 5:
				working_byte = np.bitwise_and(working_byte, 127)
				bit_index += 1
				if bit_index > 7:
					bit_index = 0
					# if working_byte > 31 and working_byte < 128:
					#	print(chr(working_byte), end=' ')
					# else:
					#	print(hex(working_byte), end=' ')
					try:
						output_buffer[byte_count] = working_byte
					except:
						pass
					else:
						byte_count += 1
				else:
					working_byte = np.right_shift(working_byte, 1)
			elif ones == 5:
				ones == 5
				# ignore stuffed zero
			elif ones == 6:
				# Frame complete
				if byte_count > 18:
					CRC = CheckCRC(output_buffer, byte_count)
					if CRC[1] == 1:
						packet_count += 1
						if verbose == 1:
							print(stream_index, hex(CRC[0]), packet_count)
							for i in range(byte_count - 2):
								character = output_buffer[i]
								if (character > 31) and (character < 128):
									print(''.join(chr(character)), end='')
							print("\r\n", end='')
				byte_count = 0
				bit_index = 0
			else:
				# Invalid frame
				byte_count = 0
				bit_index = 0
			ones = 0
	return packet_count


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

output_filter = np.array([583, 201, 157, 51, -113, -316, -525, -701, -796, -767, -580, -219, 309, 975, 1727, 2494, 3201, 3771, 4142, 4270, 4142, 3771, 3201, 2494, 1727, 975, 309, -219, -580, -767, -796, -701, -525, -316, -113, 51, 157, 201, 583])
#output_filter = np.ones(5) * 5000
output_filter_buffer = np.zeros(len(output_filter))
output_filter_shift = -3

period = 3500
attack = 3
decay = 2

Input_Fs = 28800.0
decimation = 2
Fs = Input_Fs / decimation
correlator_taps = 12

#create some dictionaries for the processing objects
InputPeakDetector = {'AttackRate':5000, 'SustainPeriod':7200, 'DecayRate':3, 'SustainCount':0, 'Envelope':0}
FilterDecimator = {'Filter':input_filter, 'DecimationRate':decimation, 'FilterBuffer':input_filter_buffer, 'DataBuffer':np.array([]), 'PeakDetector':InputPeakDetector, 'FilterShift':0, 'DecimationCounter':0, 'NewSample':0}

MarkPeakDetector = {'AttackRate':attack, 'SustainPeriod':period, 'DecayRate':decay, 'SustainCount':0, 'Envelope':0}
SpacePeakDetector = {'AttackRate':attack, 'SustainPeriod':period, 'DecayRate':decay, 'SustainCount':0, 'Envelope':0}


#set up the correlators

mark_amp = 10000
f = 2200.0

tstep = 1.0 / Fs
space_phase = 0
time = np.arange(0, tstep * correlator_taps, tstep)
time = np.add(time, space_phase)
space_amp = mark_amp
space_cos = np.rint(space_amp * (np.cos(2 * f * np.pi * time)))
space_sin = np.rint(space_amp * (np.sin(2 * f * np.pi * time)))
print(space_sin)
mark_phase = 0
f = 1200.0
time = np.arange(0, tstep * correlator_taps, tstep)
time = np.add(time, mark_phase)
mark_cos = np.rint(mark_amp * (np.cos(2 * f * np.pi * time)))
mark_sin = np.rint(mark_amp * (np.sin(2 * f * np.pi * time)))
correlator_shift = 0
correlator_shift = 2.0**(-correlator_shift)
square_scale = 18.0
square_scale = 2.0**square_scale
square_output_scale = 2.0
square_coef = 4096.0
square_clip = square_coef - 1.0

correlator_buffer = np.zeros(correlator_taps)
AFSKDemodulator1 = {'MarkCOS':mark_cos, 'MarkSIN':mark_sin, 'SpaceCOS':space_cos, 'SpaceSIN':space_sin, 'OutputFilter':output_filter, 'OutputFilterBuffer':output_filter_buffer, 'NewSample':0, 'CorrelatorBuffer':correlator_buffer, 'CorrelatorShift':correlator_shift, 'SquareScale':square_scale, 'SquareClip':square_clip, 'SquareOutputScale':square_output_scale, 'SquareCoef':square_coef, 'Result':0, 'OutputFilterShift':output_filter_shift}
DataSlicer1 = {'Rate':0.7, 'PLLClock':0.0, 'PLLStep':1000000.0, 'PLLPeriod': 12.0 * 1000000, 'LastSample':0.0, 'NewSample':0.0,'Result':0.0}
DifferentialDecoder1 = {'LastBit':0, 'NewBit':0, 'Result':0}

index1 = 0
index2 = 0
index3 = 0
envelope_index = 0

space_ratio_sum = 0

input_filter_gain = 0

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

data = np.array([])
filtered_signal_buffer = np.zeros(round(len(audio) / decimation))
demod_sig_buffer = np.zeros(round(len(audio) / decimation))
chop_audio_buffer = np.array([])
chop_filtered_audio_buffer = np.array([])
chop_demodulated_audio_buffer = np.array([])
for sample in audio:
	#chop_audio_buffer = np.append(chop_audio_buffer, np.array([sample]))
	index1 = index1 + 1
	index2 = index2 + 1
	if index2 > len(audio) / 100:
		index2 = 0
		index3 = index3 + 1
		#print(index3, InputPeakDetector['Envelope'], space_sig_ratio, space_sig_gain_error)
		print(f'{index3}')
	FilterDecimator['NewSample'] = sample
	FilterDecimator = FilterDecimate(FilterDecimator)
	for filtered_signal in FilterDecimator['DataBuffer']:
		#chop_filtered_audio_buffer = np.append(chop_filtered_audio_buffer, np.array([filtered_signal]))
		AFSKDemodulator1['NewSample'] = filtered_signal
		AFSKDemodulator1 = DemodulateAFSK(AFSKDemodulator1)
		for demodulated_signal in AFSKDemodulator1['Result']:
			#chop_demodulated_audio_buffer = np.append(chop_demodulated_audio_buffer, np.array([demodulated_signal]))
			demod_sig_buffer[envelope_index] = demodulated_signal
			envelope_index = envelope_index + 1

			#slice the data
			DataSlicer1['NewSample'] = demodulated_signal
			DataSlicer1 = ProgSliceData(DataSlicer1)
			for data_bit in DataSlicer1['Result']:
				data = np.append(data, np.array([data_bit]))
				# DifferentialDecoder1['NewBit'] = data_bit
				# DifferentialDecoder1 = ProgDifferentialDecode(DifferentialDecoder1)

				#data = np.append(data, np.array([DifferentialDecoder1['Result']]))
# scipy.io.wavfile.write(dirname+"PeakDetect.wav", round(samplerate / decimation), envelope.astype(np.int16))
# scipy.io.wavfile.write(dirname+"FilteredSignal.wav", round(samplerate / decimation), filtered_signal_buffer.astype(np.int16))

# scipy.io.wavfile.write(dirname+"MarkCorrelatorSignal.wav", round(samplerate / decimation), mark_correlator_buffer.astype(np.int16))
# scipy.io.wavfile.write(dirname+"SpaceCorrelatorSignal.wav", round(samplerate / decimation), space_correlator_buffer.astype(np.int16))
scipy.io.wavfile.write(dirname+"DemodSignal.wav", round(samplerate / decimation), demod_sig_buffer.astype(np.int16))

# scipy.io.wavfile.write(dirname+"SpaceGain.wav", round(samplerate / decimation), space_gain_buffer.astype(np.int16))

# scipy.io.wavfile.write(dirname+"SpaceEnvelope.wav", round(samplerate / decimation), space_envelope_buffer.astype(np.int16))

# scipy.io.wavfile.write(dirname+"MarkEnvelope.wav", round(samplerate / decimation), mark_envelope_buffer.astype(np.int16))

#data = SliceData(demod_sig_buffer, 0.7, 12)
#print("Length of sliced data:", len(data))
# file = open('demod_output.txt', 'w')
# np.savetxt(file, data.astype(np.int16))
# file.close()
data = DifferentialDecode(data)
count = DecodeAX25(data, 1)
print("Decoded packet count:", count)

print('made new directory: ', dirname)
print('done')
