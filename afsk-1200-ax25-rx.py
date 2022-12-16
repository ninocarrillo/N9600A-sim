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
	if signal_value > detector['High']:
		detector['High'] = detector['High'] + detector['AttackRate']
		if detector['High'] > signal_value:
			detector['High'] = signal_value
		detector['Midpoint'] = ((detector['High'] - detector['Low']) // 2) + detector['Low']
		detector['HighSustainCount'] = 0
	if signal_value < detector['Low']:
		detector['Low'] = detector['Low'] - detector['AttackRate']
		if detector['Low'] < signal_value:
			detector['Low'] = signal_value
		detector['Midpoint'] = ((detector['High'] - detector['Low']) // 2) + detector['Low']
		detector['LowSustainCount'] = 0
	detector['LastValue'] = signal_value
	if detector['HighSustainCount'] >= detector['SustainPeriod']:
		detector['High'] = detector['High'] - detector['DecayRate']
		if detector['High'] <= 0:
			detector['High'] = 1
			detector['HighSustainCount'] = 0
		detector['Midpoint'] = ((detector['High'] - detector['Low']) // 2) + detector['Low']
	if detector['LowSustainCount'] >= detector['SustainPeriod']:
		detector['Low'] = detector['Low'] + detector['DecayRate']
		if detector['Low'] >= 0:
			detector['Low'] = -1
			detector['LowSustainCount'] = 0
		detector['Midpoint'] = ((detector['High'] - detector['Low']) // 2) + detector['Low']
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
			filter['GainChange'] = 0
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
	# if mark_sig > 4095:
	# 	print('mark clip')
	mark_sig = np.clip(mark_sig, 0, demodulator['SquareClip'])

	mark_sig = np.rint(demodulator['SquareOutputScale'] * np.sqrt(demodulator['SquareCoef'] * mark_sig))

	space_cos_sig = np.rint((demodulator['CorrelatorShift'] * np.convolve(demodulator['CorrelatorBuffer'], demodulator['SpaceCOS'], 'valid')) / pow(2, (16 + demodulator['CorrelatorShift'])))
	space_sin_sig = np.rint((demodulator['CorrelatorShift'] * np.convolve(demodulator['CorrelatorBuffer'], demodulator['SpaceSIN'], 'valid')) / pow(2, (16 + demodulator['CorrelatorShift'])))

	space_sig = np.add(np.square(space_cos_sig), np.square(space_sin_sig))
	space_sig = np.rint(space_sig / demodulator['SquareScale'])
	# if space_sig > 4095:
	# 	print('space clip')
	space_sig = np.clip(space_sig, 0, demodulator['SquareClip'])

	space_sig = np.rint(demodulator['SquareOutputScale']* np.sqrt(demodulator['SquareCoef'] * space_sig))
	space_sig = np.rint(space_sig * demodulator['SpaceRatio'])

	demodulator['OutputFilterBuffer'] = demodulator['OutputFilterBuffer'][1:]
	demodulator['OutputFilterBuffer'] = np.append(demodulator['OutputFilterBuffer'], np.array([mark_sig - space_sig]))

	demodulator['Result'] = np.rint(np.convolve(demodulator['OutputFilterBuffer'], demodulator['OutputFilter'], 'valid') / pow(2, (16 + demodulator['OutputFilterShift'])))
	demodulator['EnvelopeDetector'] = HighLowDetect(demodulator['Result'], demodulator['EnvelopeDetector'])
	demodulator['Result'] = demodulator['Result'] - (demodulator['EnvelopeDetector']['Midpoint'] * 0.5)
	return demodulator


def ProgDifferentialDecode(decoder):
	if decoder['NewBit'] == decoder['LastBit']:
		decoder['Result'] = 1
	else:
		decoder['Result'] = 0
	decoder['LastBit'] = decoder['NewBit']
	return decoder


def ProgDecodeAX25(decoder):
	if decoder['NewBit'] == 1:
		decoder['WorkingByte'] = np.bitwise_or(decoder['WorkingByte'], 128)
		decoder['Ones'] += 1
		decoder['BitIndex'] += 1
		if decoder['Ones'] > 6:
			# abort frame for invalid bit sequence
			decoder['Ones'] = 0
			decoder['BitIndex'] = 0
			decoder['ByteCount'] = 0
			decoder['Result'] = np.array([]).astype('uint16')
		elif decoder['BitIndex'] >7:
			# 8 valid bits received, record byte
			decoder['BitIndex'] = 0
			decoder['Result'] = np.append(decoder['Result'], np.array([decoder['WorkingByte']]).astype('uint16'))
			decoder['ByteCount'] += 1
		else:
			decoder['WorkingByte'] = np.right_shift(decoder['WorkingByte'], 1)
	else:
		if decoder['Ones'] < 5:
			decoder['WorkingByte'] = np.bitwise_and(decoder['WorkingByte'], 127)
			decoder['BitIndex'] += 1
			if decoder['BitIndex'] > 7:
				decoder['BitIndex'] = 0
				decoder['Result'] = np.append(decoder['Result'], np.array([decoder['WorkingByte']]).astype('uint16'))
				decoder['ByteCount'] += 1
			else:
				decoder['WorkingByte'] = np.right_shift(decoder['WorkingByte'], 1)
		elif decoder['Ones'] == 5:
			pass
			# ignore stuffed zero
		elif decoder['Ones'] == 6:
			# Frame complete
			if decoder['ByteCount'] > 18:
				decoder['CRC'] = CheckCRC(decoder['Result'].astype('uint16'), len(decoder['Result']))
				if decoder['CRC'][1] == 1:
					decoder['PacketCount'] += 1
					decoder['Output'] = decoder['Result']
					decoder['OutputTrigger'] = True
					if decoder['Verbose'] == 1:
						print(hex(decoder['CRC'][0]), decoder['PacketCount'])
						for character in decoder['Result'][:-2]:
							try:
								print(''.join(chr(character)), end='')
							except:
								pass
			decoder['ByteCount'] = 0
			decoder['Result'] = np.array([]).astype('uint16')
			decoder['BitIndex'] = 0
		else:
			decoder['ByteCount'] = 0
			decoder['BitIndex'] = 0
			decoder['Result'] = np.array([]).astype('uint16')
		decoder['Ones'] = 0
	return decoder


def ProgSliceData(slicer):
	# slicer['EnvelopeDetector'] = HighLowDetect(slicer['NewSample'], slicer['EnvelopeDetector'])
	# slicer['Midpoint'] = np.rint(slicer['EnvelopeDetector']['Midpoint'] * 0.66)
	slicer['Midpoint'] = 0
	slicer['Result'] = np.array([])
	slicer['PLLClock'] += slicer['PLLStep']
	if slicer['PLLClock'] > ((slicer['PLLPeriod'] / 2.0) - 1.0):
		slicer['PLLClock'] -= slicer['PLLPeriod']
		if slicer['NewSample'] > slicer['Midpoint']:
			slicer['Result'] = np.array([1])
		else:
			slicer['Result'] = np.array([0])
	if slicer['LastSample'] > slicer['Midpoint']:
		if slicer['NewSample'] <= slicer['Midpoint']:
			# Zero Crossing
			slicer['PLLClock'] *= slicer['Rate']
	else:
		if slicer['NewSample'] > slicer['Midpoint']:
			# Zero Crossing
			slicer['PLLClock'] *= slicer['Rate']
	slicer['LastSample'] = slicer['NewSample']
	return slicer

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



# /*
#
# FIR filter designed with
# http://t-filter.appspot.com
#
# sampling frequency: 28800 Hz
#
# fixed point precision: 16 bits
#
# * 0 Hz - 300 Hz
#   gain = 0
#   desired attenuation = -70 dB
#   actual attenuation = n/a
#
# * 1000 Hz - 2400 Hz
#   gain = 1
#   desired ripple = 2 dB
#   actual ripple = n/a
#
# * 3600 Hz - 14400 Hz
#   gain = 0
#   desired attenuation = -70 dB
#   actual attenuation = n/a
#
# */


#input_filter = np.array([-20, -49, -93, -143, -184, -199, -171, -95, 20, 149, 262, 330, 340, 303, 247, 207, 213, 266, 338, 382, 347, 206, -22, -276, -471, -533, -441, -244, -62, -41, -294, -844, -1585, -2297, -2698, -2541, -1709, -272, 1507, 3229, 4476, 4931, 4476, 3229, 1507, -272, -1709, -2541, -2698, -2297, -1585, -844, -294, -41, -62, -244, -441, -533, -471, -276, -22, 206, 347, 382, 338, 266, 213, 207, 247, 303, 340, 330, 262, 149, 20, -95, -171, -199, -184, -143, -93, -49, -20])
input_filter_buffer = np.zeros(len(input_filter))

output_filter = np.array([583, 201, 157, 51, -113, -316, -525, -701, -796, -767, -580, -219, 309, 975, 1727, 2494, 3201, 3771, 4142, 4270, 4142, 3771, 3201, 2494, 1727, 975, 309, -219, -580, -767, -796, -701, -525, -316, -113, 51, 157, 201, 583])
#output_filter = np.ones(5) * 5000
output_filter_buffer = np.zeros(len(output_filter))
output_filter_shift = -3

period = 3500
attack = 3
decay = 2

Input_Fs = 28800
decimation = 2
Fs = Input_Fs // decimation
correlator_taps = 12

#create some dictionaries for the processing objects
InputPeakDetector = {'AttackRate':50, 'SustainPeriod':3000, 'DecayRate':50, 'SustainCount':0, 'Envelope':0}
FilterDecimator = {'Filter':input_filter, 'DecimationRate':decimation, 'FilterBuffer':input_filter_buffer, 'DataBuffer':np.array([]), 'PeakDetector':InputPeakDetector, 'FilterShift':0, 'DecimationCounter':0, 'NewSample':0}

#generate the correlator window fuunction
N = correlator_taps
L = N + 1
window =np.zeros(N)
for n in range(N):
	# window[n] = 1.0 - np.abs((n-(N/2)) / (L/2))

	# this is the Hann window
	window[n] = 0.5*(1-np.cos(2*np.pi*n/N))
	window[n] = 1

#set up the correlators

mark_amp = 15000
f = 2200.0

tstep = 1.0 / Fs
space_phase = 0
time = np.arange(0, tstep * correlator_taps, tstep)
time = np.add(time, space_phase)
space_amp = mark_amp
space_cos = np.rint(space_amp * (np.cos(2 * f * np.pi * time)))
space_sin = np.rint(space_amp * (np.sin(2 * f * np.pi * time)))
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

space_cos = np.rint(np.multiply(space_cos,window))
space_sin = np.rint(np.multiply(space_sin,window))
mark_cos = np.rint(np.multiply(mark_cos,window))
mark_sin = np.rint(np.multiply(mark_sin,window))
print(mark_sin)
print(space_sin)
print(mark_cos)
print(space_cos)
print(window)

correlator_buffer = np.zeros(correlator_taps)

SlicerEnvelope = {'AttackRate':1, 'DecayRate':1, 'SustainPeriod':450, 'High':0, 'Low':0, 'HighSustainCount':0, 'LowSustainCount':0, 'Midpoint':0}
AFSKDemodulator1 = {'MarkCOS':mark_cos, 'MarkSIN':mark_sin, 'SpaceCOS':space_cos, 'SpaceSIN':space_sin, 'SpaceRatio':1.0, 'OutputFilter':output_filter, 'OutputFilterBuffer':output_filter_buffer, 'NewSample':0, 'CorrelatorBuffer':correlator_buffer, 'CorrelatorShift':correlator_shift, 'SquareScale':square_scale, 'SquareClip':square_clip, 'SquareOutputScale':square_output_scale, 'SquareCoef':square_coef, 'Result':0, 'OutputFilterShift':output_filter_shift, 'EnvelopeDetector':SlicerEnvelope}
AFSKDemodulator2 = {'MarkCOS':mark_cos, 'MarkSIN':mark_sin, 'SpaceCOS':space_cos, 'SpaceSIN':space_sin, 'SpaceRatio':1.0, 'OutputFilter':output_filter, 'OutputFilterBuffer':output_filter_buffer, 'NewSample':0, 'CorrelatorBuffer':correlator_buffer, 'CorrelatorShift':correlator_shift, 'SquareScale':square_scale, 'SquareClip':square_clip, 'SquareOutputScale':square_output_scale, 'SquareCoef':square_coef, 'Result':0, 'OutputFilterShift':output_filter_shift, 'EnvelopeDetector':SlicerEnvelope}
DataSlicer1 = {'Rate':0.7, 'PLLClock':0.0, 'PLLStep':1000000.0, 'PLLPeriod': 12.0 * 1000000, 'LastSample':0.0, 'NewSample':0.0,'Result':0.0, 'Midpoint':0, 'EnvelopeDetector':SlicerEnvelope}
DataSlicer2 = {'Rate':0.7, 'PLLClock':0.0, 'PLLStep':1000000.0, 'PLLPeriod': 12.0 * 1000000, 'LastSample':0.0, 'NewSample':0.0,'Result':0.0, 'Midpoint':0, 'EnvelopeDetector':SlicerEnvelope}
DifferentialDecoder1 = {'LastBit':0, 'NewBit':0, 'Result':0}
DifferentialDecoder2 = {'LastBit':0, 'NewBit':0, 'Result':0}
AX25Decoder1 = {'NewBit':0, 'BitIndex':0, 'Ones':0, 'ByteCount':0, 'WorkingByte':np.uint16(0), 'Result':np.array([]).astype('uint16'), 'CRC':np.array([]), 'PacketCount':0, 'Verbose':0, 'OutputTrigger':False}
AX25Decoder2 = {'NewBit':0, 'BitIndex':0, 'Ones':0, 'ByteCount':0, 'WorkingByte':np.uint16(0), 'Result':np.array([]).astype('uint16'), 'CRC':np.array([]), 'PacketCount':0, 'Verbose':0, 'OutputTrigger':False}


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
	chop_audio_buffer = np.append(chop_audio_buffer, np.array([sample]))
	index1 = index1 + 1
	index2 = index2 + 1
	if index2 > len(audio) / 100:
		index2 = 0
		index3 = index3 + 1
		midpoint = DataSlicer1['Midpoint']
		#print(index3, InputPeakDetector['Envelope'], space_sig_ratio, space_sig_gain_error)
		print(f'{index3}')
	FilterDecimator['NewSample'] = sample
	FilterDecimator = FilterDecimate(FilterDecimator)

	for filtered_signal in FilterDecimator['DataBuffer']:
		chop_filtered_audio_buffer = np.append(chop_filtered_audio_buffer, np.array([filtered_signal]))
		AFSKDemodulator1['NewSample'] = filtered_signal
		AFSKDemodulator1 = DemodulateAFSK(AFSKDemodulator1)
		for demodulated_signal in AFSKDemodulator1['Result']:
			chop_demodulated_audio_buffer = np.append(chop_demodulated_audio_buffer, np.array([demodulated_signal]))
			demod_sig_buffer[envelope_index] = demodulated_signal
			envelope_index = envelope_index + 1

			#slice the data
			DataSlicer1['NewSample'] = demodulated_signal
			DataSlicer1 = ProgSliceData(DataSlicer1)
			for data_bit in DataSlicer1['Result']:
				# data = np.append(data, np.array([data_bit]))
				DifferentialDecoder1['NewBit'] = data_bit
				DifferentialDecoder1 = ProgDifferentialDecode(DifferentialDecoder1)

				#data = np.append(data, np.array([DifferentialDecoder1['Result']]))
				AX25Decoder1['NewBit'] = DifferentialDecoder1['Result']
				AX25Decoder1 = ProgDecodeAX25(AX25Decoder1)
				if AX25Decoder1['OutputTrigger'] == True:
					AX25Decoder1['OutputTrigger'] = False

					packet = AX25Decoder1['PacketCount']
					CRC = AX25Decoder1['CRC'][0]
					filename = f'Packet-{packet}_CRC-{format(CRC,"#06x")}_Index-{index1}'
					print(dirname+filename)
					scipy.io.wavfile.write(dirname+filename+'-audio.wav', Input_Fs, chop_audio_buffer.astype(np.int16))
					chop_audio_buffer = np.array([])
					scipy.io.wavfile.write(dirname+filename+'-demod.wav', Fs, chop_demodulated_audio_buffer.astype(np.int16))
					chop_demodulated_audio_buffer = np.array([])
					scipy.io.wavfile.write(dirname+filename+'-filtered.wav', Fs, chop_filtered_audio_buffer.astype(np.int16))
					chop_filtered_audio_buffer = np.array([])



scipy.io.wavfile.write(dirname+"DemodSignal.wav", round(samplerate / decimation), demod_sig_buffer.astype(np.int16))


print('made new directory: ', dirname)
print('done')
