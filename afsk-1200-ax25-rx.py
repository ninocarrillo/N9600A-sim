import sys
import struct
import scipy.io.wavfile
import numpy as np
import os
import n9600a_progdemod as demod

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
output_filter = np.ones(3) * 5000
output_filter_buffer = np.zeros(len(output_filter))
output_filter_buffer2 = np.zeros(len(output_filter))
output_filter_shift = -3

period = 3500
attack = 3
decay = 2

Input_Fs = 28800
decimation = 2
Fs = Input_Fs // decimation
correlator_taps = 12
samplesperbit = Fs // 1200

#create some dictionaries for the processing objects
InputPeakDetector = {'AttackRate':50, 'SustainPeriod':2000, 'DecayRate':50, 'SustainCount':0, 'Envelope':0}
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

mark_amp = 20000
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
correlator_buffer2 = np.zeros(correlator_taps)

SlicerEnvelope = {'AttackRate':1, 'DecayRate':1, 'SustainPeriod':300, 'High':0, 'Low':0, 'HighSustainCount':0, 'LowSustainCount':0, 'Midpoint':0}
AFSKDemodulator1 = {'MarkCOS':mark_cos, 'MarkSIN':mark_sin, 'SpaceCOS':space_cos, 'SpaceSIN':space_sin, 'SpaceRatio':1.0, 'OutputFilter':output_filter, 'OutputFilterBuffer':output_filter_buffer, 'NewSample':0, 'CorrelatorBuffer':correlator_buffer, 'CorrelatorShift':correlator_shift, 'SquareScale':square_scale, 'SquareClip':square_clip, 'SquareOutputScale':square_output_scale, 'SquareCoef':square_coef, 'Result':0, 'OutputFilterShift':output_filter_shift, 'EnvelopeDetector':SlicerEnvelope}
AFSKDemodulator2 = {'MarkCOS':mark_cos, 'MarkSIN':mark_sin, 'SpaceCOS':space_cos, 'SpaceSIN':space_sin, 'SpaceRatio':1.39, 'OutputFilter':output_filter, 'OutputFilterBuffer':output_filter_buffer, 'NewSample':0, 'CorrelatorBuffer':correlator_buffer, 'CorrelatorShift':correlator_shift, 'SquareScale':square_scale, 'SquareClip':square_clip, 'SquareOutputScale':square_output_scale, 'SquareCoef':square_coef, 'Result':0, 'OutputFilterShift':output_filter_shift, 'EnvelopeDetector':SlicerEnvelope}
DataSlicer1 = {'Rate':0.7, 'PLLClock':0.0, 'PLLStep':1000000.0, 'PLLPeriod': samplesperbit * 1000000, 'LastSample':0.0, 'NewSample':0.0,'Result':0.0, 'Midpoint':0, 'EnvelopeDetector':SlicerEnvelope}
DataSlicer2 = {'Rate':0.7, 'PLLClock':0.0, 'PLLStep':1000000.0, 'PLLPeriod': samplesperbit * 1000000, 'LastSample':0.0, 'NewSample':0.0,'Result':0.0, 'Midpoint':0, 'EnvelopeDetector':SlicerEnvelope}
DifferentialDecoder1 = {'LastBit':0, 'NewBit':0, 'Result':0}
DifferentialDecoder2 = {'LastBit':0, 'NewBit':0, 'Result':0}
AX25Decoder1 = {'NewBit':0, 'BitIndex':0, 'Ones':0, 'ByteCount':0, 'WorkingByte':np.uint16(0), 'Result':np.array([]).astype('uint16'), 'CRC':np.array([0,0]), 'PacketCount':0, 'Verbose':0, 'OutputTrigger':False, 'CRCAge':1000000, 'UniquePackets':0}
AX25Decoder2 = {'NewBit':0, 'BitIndex':0, 'Ones':0, 'ByteCount':0, 'WorkingByte':np.uint16(0), 'Result':np.array([]).astype('uint16'), 'CRC':np.array([0,0]), 'PacketCount':0, 'Verbose':0, 'OutputTrigger':False, 'CRCAge':1000000, 'UniquePackets':0}

print(AFSKDemodulator1)

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


total_packets = 0
duplicate_packets = 0

data = np.array([])
filtered_signal_buffer = np.zeros(round(len(audio) / decimation))
demod_sig_buffer = np.zeros(round(len(audio) / decimation))
mark_sig_buffer = np.zeros(round(len(audio) / decimation))
space_sig_buffer = np.zeros(round(len(audio) / decimation))
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
	FilterDecimator = demod.FilterDecimate(FilterDecimator)

	for filtered_signal in FilterDecimator['DataBuffer']:
		chop_filtered_audio_buffer = np.append(chop_filtered_audio_buffer, np.array([filtered_signal]))
		AFSKDemodulator1['NewSample'] = filtered_signal
		filtered_signal_buffer[envelope_index] = filtered_signal
		AFSKDemodulator1 = demod.DemodulateAFSK(AFSKDemodulator1)
		for demodulated_signal in AFSKDemodulator1['Result']:
			chop_demodulated_audio_buffer = np.append(chop_demodulated_audio_buffer, np.array([demodulated_signal]))
			demod_sig_buffer[envelope_index] = demodulated_signal
			mark_sig_buffer[envelope_index] = AFSKDemodulator1['MarkSig']
			space_sig_buffer[envelope_index] = AFSKDemodulator1['SpaceSig']
			envelope_index = envelope_index + 1

			#slice the data
			DataSlicer1['NewSample'] = demodulated_signal
			DataSlicer1 = demod.ProgSliceData(DataSlicer1)
			for data_bit in DataSlicer1['Result']:
				# data = np.append(data, np.array([data_bit]))
				DifferentialDecoder1['NewBit'] = data_bit
				DifferentialDecoder1 = demod.ProgDifferentialDecode(DifferentialDecoder1)

				#data = np.append(data, np.array([DifferentialDecoder1['Result']]))
				AX25Decoder1['NewBit'] = DifferentialDecoder1['Result']
				AX25Decoder1 = demod.ProgDecodeAX25(AX25Decoder1)
				if AX25Decoder1['OutputTrigger'] == True:
					AX25Decoder1['OutputTrigger'] = False
					# Check for unioqueness
					if AX25Decoder2['CRCAge'] > 10 or AX25Decoder1['CRC'][0] != AX25Decoder2['CRC'][0]:
						total_packets += 1
						CRC = AX25Decoder1['CRC'][0]
						decodernum = '1'
						filename = f'Packet-{total_packets}_CRC-{format(CRC,"#06x")}_decoder-{decodernum}_Index-{index1}'
						print(dirname+filename)
						scipy.io.wavfile.write(dirname+filename+'-audio.wav', Input_Fs, chop_audio_buffer.astype(np.int16))
						chop_audio_buffer = np.array([])
						scipy.io.wavfile.write(dirname+filename+'-demod.wav', Fs, chop_demodulated_audio_buffer.astype(np.int16))
						chop_demodulated_audio_buffer = np.array([])
						scipy.io.wavfile.write(dirname+filename+'-filtered.wav', Fs, chop_filtered_audio_buffer.astype(np.int16))
						chop_filtered_audio_buffer = np.array([])
					else:
						if AX25Decoder1['CRC'][0] == AX25Decoder2['CRC'][0]:
							duplicate_packets += 1
							AX25Decoder1['UniquePackets'] -= 1
							AX25Decoder2['UniquePackets'] -= 1
							print('Decoder 1 Duplicate, bit delay: ', AX25Decoder2['CRCAge'])



		# Second AFSKDemodulator
		AFSKDemodulator2['NewSample'] = filtered_signal
		AFSKDemodulator2 = demod.DemodulateAFSK(AFSKDemodulator2)
		for demodulated_signal in AFSKDemodulator2['Result']:
			DataSlicer2['NewSample'] = demodulated_signal
			DataSlicer2 = demod.ProgSliceData(DataSlicer2)
			for data_bit in DataSlicer2['Result']:
				DifferentialDecoder2['NewBit'] = data_bit
				DifferentialDecoder2 = demod.ProgDifferentialDecode(DifferentialDecoder2)
				AX25Decoder2['NewBit'] = DifferentialDecoder2['Result']
				AX25Decoder2 = demod.ProgDecodeAX25(AX25Decoder2)
				if AX25Decoder2['OutputTrigger'] == True:
					AX25Decoder2['OutputTrigger'] = False
					# Check for uniqueness
					if AX25Decoder1['CRCAge'] > 10 or AX25Decoder1['CRC'][0] != AX25Decoder2['CRC'][0]:
						total_packets += 1
						CRC = AX25Decoder1['CRC'][0]
						decodernum = '2'
						filename = f'Packet-{total_packets}_CRC-{format(CRC,"#06x")}_decoder-{decodernum}_Index-{index1}'
						print(dirname+filename)
						scipy.io.wavfile.write(dirname+filename+'-audio.wav', Input_Fs, chop_audio_buffer.astype(np.int16))
						chop_audio_buffer = np.array([])
						scipy.io.wavfile.write(dirname+filename+'-demod.wav', Fs, chop_demodulated_audio_buffer.astype(np.int16))
						chop_demodulated_audio_buffer = np.array([])
						scipy.io.wavfile.write(dirname+filename+'-filtered.wav', Fs, chop_filtered_audio_buffer.astype(np.int16))
						chop_filtered_audio_buffer = np.array([])
					else:
						if AX25Decoder1['CRC'][0] == AX25Decoder2['CRC'][0]:
							duplicate_packets += 1
							AX25Decoder1['UniquePackets'] -= 1
							AX25Decoder2['UniquePackets'] -= 1
							print('Decoder 2 Duplicate, bit delay: ', AX25Decoder1['CRCAge'])

print(AFSKDemodulator1)
scipy.io.wavfile.write(dirname+"DemodSignal.wav", round(samplerate / decimation), demod_sig_buffer.astype(np.int16))
scipy.io.wavfile.write(dirname+"FilteredSignal.wav", round(samplerate / decimation), filtered_signal_buffer.astype(np.int16))
scipy.io.wavfile.write(dirname+"MarkSignal.wav", round(samplerate / decimation), mark_sig_buffer.astype(np.int16))
scipy.io.wavfile.write(dirname+"SpaceSignal.wav", round(samplerate / decimation), space_sig_buffer.astype(np.int16))


print('total packets: ', total_packets)
print('duplicate_packets: ', duplicate_packets)
print('Decoder 1 unique packets: ', AX25Decoder1['UniquePackets'])
print('Decoder 1 total packets: ', AX25Decoder1['PacketCount'])
print('Decoder 2 unique packets: ', AX25Decoder2['UniquePackets'])
print('Decoder 2 total packets: ', AX25Decoder1['PacketCount'])
print('made new directory: ', dirname)
print('done')
