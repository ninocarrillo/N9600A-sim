import sys
import configparser
import struct
import scipy.io.wavfile
import numpy as np
import os
import n9600a_progdemod as demod

def StringToIntArray(input_string):
	working_string = ''
	result = np.array([])
	for character in input_string:
		if (character >= '0' and character <= '9') or ((character == '-') or (character == '.')):
			working_string += character
		else:
			if working_string != '':
				result = np.append(result, np.array([int(working_string)]))
				working_string = ''
	return result

if len(sys.argv) < 3:
	print("Not enough arguments. Usage: py -3 n9600a_rx.py <ini file> <wav file>")
	sys.exit(-1)

# read demodulator description from ini file:
config = configparser.ConfigParser()
try:
	config.read(sys.argv[1])
except:
	print(f'Unable to open ini file {sys.argv[1]}')
	sys.exit(-2)

# Get input filter taps
FilterDecimator = {}
print(f'Opened ini file {sys.argv[1]}')
try:
 	FilterDecimator['Filter'] = StringToIntArray(config['Input Filter']['taps'])
except:
	print(f'{sys.argv[1]} [Input Filter] \'taps\' is missing or invalid')
	sys.exit(-2)

# Get input filter Sample Rate:
try:
	FilterDecimator['InputSampleRate'] = int(config['Input Filter']['sample rate'])
except:
	print(f'{sys.argv[1]} [Input Filter] \'sample rate\' is missing or invalid')
	sys.exit(-2)

# Get input filter Decimation Rate:
try:
	FilterDecimator['DecimationRate'] = int(config['Input Filter']['decimation'])
except:
	print(f'{sys.argv[1]} [Input Filter] \'decimation\' is missing or invalid')
	sys.exit(-2)

# Get input filter AGC option:
try:
	FilterDecimator['InputAGCEnabled'] = config['Input Filter'].getboolean('agc enabled')
except:
	print(f'{sys.argv[1]} [Input Filter] \'agc enabled\' is missing or invalid')
	sys.exit(-2)

# Get Input Filter AGC Attack Rate:
try:
	FilterDecimator['InputAGCAttackRate'] = int(config['Input Filter']['agc attack rate'])
except:
	print(f'{sys.argv[1]} [Input Filter] \'agc attack rate\' is missing or invalid')
	sys.exit(-2)

# Get Input Filter AGC Sustain Period:
try:
	FilterDecimator['InputAGCSustainPeriod'] = int(config['Input Filter']['agc sustain period'])
except:
	print(f'{sys.argv[1]} [Input Filter] \'agc sustain period\' is missing or invalid')
	sys.exit(-2)

# Get Input Filter AGC Decay Rate:
try:
	FilterDecimator['InputAGCDecayRate'] = int(config['Input Filter']['agc decay rate'])
except:
	print(f'{sys.argv[1]} [Input Filter] \'agc decay rate\' is missing or invalid')
	sys.exit(-2)

# Get Input Filter ADC bit count:
try:
	FilterDecimator['InputBitCount'] = int(config['Input Filter']['input bit count'])
except:
	print(f'{sys.argv[1]} [Input Filter] \'input bit count\' is missing or invalid')
	sys.exit(-2)

# Create the Input Filter dictionaries
FilterDecimator = demod.InitFilterDecimator(FilterDecimator)

# Read settings for AFSK Demodulator 1
AFSKDemodulator1 = {}
try:
	AFSKDemodulator1['Enabled'] = config['AFSK Demodulator 1'].getboolean('enabled')
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 1] \'enabled\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator1['CorrelatorTapCount'] = int(config['AFSK Demodulator 1']['correlator tap count'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 1] \'correlator tap count\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator1['MarkAmplitude'] = float(config['AFSK Demodulator 1']['mark amplitude'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 1] \'mark amplitude\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator1['SpaceAmplitude'] = float(config['AFSK Demodulator 1']['space amplitude'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 1] \'space amplitude\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator1['SpaceRatio'] = float(config['AFSK Demodulator 1']['space gain'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 1] \'space gain\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator1['OffsetRemovalEnabled'] = config['AFSK Demodulator 1'].getboolean('offset removal enabled')
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 1] \'offset removal enabled\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator1['MarkFreq'] = int(config['AFSK Demodulator 1']['mark frequency'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 1] \'mark frequency\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator1['SpaceFreq'] = int(config['AFSK Demodulator 1']['space frequency'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 1] \'space frequency\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator1['OutputFilter'] = StringToIntArray(config['AFSK Demodulator 1']['output filter taps'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 1] \'output filter taps\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator1['OutputFilterShift'] = int(config['AFSK Demodulator 1']['output filter shift'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 1] \'output filter shift\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator1['EnvelopeAttackRate'] = int(config['AFSK Demodulator 1']['envelope attack rate'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 1] \'envelope attack rate\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator1['OffsetRemovalRate'] = float(config['AFSK Demodulator 1']['offset removal rate'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 1] \'offset removal rate\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator1['EnvelopeSustainPeriod'] = int(config['AFSK Demodulator 1']['envelope sustain period'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 1] \'envelope sustain period\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator1['EnvelopeDecayRate'] = int(config['AFSK Demodulator 1']['envelope decay rate'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 1] \'envelope decay rate\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator1['ChopAudio'] = config['AFSK Demodulator 1'].getboolean('chop audio enable')
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 1] \'chop audio enable\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator1['SqrtBitCount'] = int(config['AFSK Demodulator 1']['sqrt bit count'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 1] \'sqrt bit count\' is missing or invalid')
	sys.exit(-2)

DataSlicer1 = {}
try:
	DataSlicer1['BitRate'] = int(config['AFSK Demodulator 1']['slicer bit rate'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 1] \'slicer bit rate\' is missing or invalid')
	sys.exit(-2)
try:
	DataSlicer1['Rate'] = float(config['AFSK Demodulator 1']['slicer lock rate'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 1] \'slicer lock rate\' is missing or invalid')
	sys.exit(-2)

# Initialize AFSK Demodulator 1
AFSKDemodulator1['InputSampleRate'] = FilterDecimator['OutputSampleRate']
DataSlicer1['InputSampleRate'] = FilterDecimator['OutputSampleRate']
AFSKDemodulator1 = demod.InitAFSKDemod(AFSKDemodulator1)
DataSlicer1 = demod.InitDataSlicer(DataSlicer1)

# Read settings for AFSK Demodulator 1
AFSKDemodulator2 = {}
try:
	AFSKDemodulator2['Enabled'] = config['AFSK Demodulator 2'].getboolean('enabled')
except:
	print(f'{sys.argv[1]} [AFSKDemodulator2] \'enabled\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator2['CorrelatorTapCount'] = int(config['AFSK Demodulator 2']['correlator tap count'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 2] \'correlator tap count\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator2['MarkAmplitude'] = float(config['AFSK Demodulator 2']['mark amplitude'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 2] \'mark amplitude\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator2['SpaceAmplitude'] = float(config['AFSK Demodulator 2']['space amplitude'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 2] \'space amplitude\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator2['SpaceRatio'] = float(config['AFSK Demodulator 2']['space gain'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 2] \'space gain\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator2['OffsetRemovalEnabled'] = config['AFSK Demodulator 2'].getboolean('offset removal enabled')
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 2] \'offset removal enabled\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator2['MarkFreq'] = int(config['AFSK Demodulator 2']['mark frequency'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 2] \'mark frequency\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator2['SpaceFreq'] = int(config['AFSK Demodulator 2']['space frequency'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 2] \'space frequency\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator2['OutputFilter'] = StringToIntArray(config['AFSK Demodulator 2']['output filter taps'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 2] \'output filter taps\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator2['OutputFilterShift'] = int(config['AFSK Demodulator 2']['output filter shift'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 2] \'output filter shift\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator2['OffsetRemovalRate'] = float(config['AFSK Demodulator 2']['offset removal rate'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 2] \'offset removal rate\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator2['EnvelopeAttackRate'] = int(config['AFSK Demodulator 2']['envelope attack rate'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 2] \'envelope attack rate\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator2['EnvelopeSustainPeriod'] = int(config['AFSK Demodulator 2']['envelope sustain period'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 2] \'envelope sustain period\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator2['EnvelopeDecayRate'] = int(config['AFSK Demodulator 2']['envelope decay rate'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 2] \'envelope decay rate\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator2['ChopAudio'] = config['AFSK Demodulator 2'].getboolean('chop audio enable')
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 2] \'chop audio enable\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator2['SqrtBitCount'] = int(config['AFSK Demodulator 2']['sqrt bit count'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 2] \'sqrt bit count\' is missing or invalid')
	sys.exit(-2)

DataSlicer2 = {}
try:
	DataSlicer2['BitRate'] = int(config['AFSK Demodulator 2']['slicer bit rate'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 2] \'slicer bit rate\' is missing or invalid')
	sys.exit(-2)
try:
	DataSlicer2['Rate'] = float(config['AFSK Demodulator 2']['slicer lock rate'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 2] \'slicer lock rate\' is missing or invalid')
	sys.exit(-2)

# Initialize AFSK Demodulator 2
AFSKDemodulator2['InputSampleRate'] = FilterDecimator['OutputSampleRate']
DataSlicer2['InputSampleRate'] = FilterDecimator['OutputSampleRate']
AFSKDemodulator2 = demod.InitAFSKDemod(AFSKDemodulator2)
DataSlicer2 = demod.InitDataSlicer(DataSlicer2)

try:
	samplerate, audio = scipy.io.wavfile.read(sys.argv[2])
	# Take two bits of resolution away
	audio = audio >> (16 - FilterDecimator['InputBitCount'])
except:
	print('Unable to open wave file.')
	sys.exit(-2)

print("Opened file. \r\nSample rate:", samplerate, "\r\nLength:", len(audio))

# space_cos = np.rint(np.multiply(space_cos,window))
# space_sin = np.rint(np.multiply(space_sin,window))
# mark_cos = np.rint(np.multiply(mark_cos,window))
# mark_sin = np.rint(np.multiply(mark_sin,window))
# print(mark_sin)
# print(space_sin)
# print(mark_cos)
# print(space_cos)
# print(window)

DifferentialDecoder1 = demod.InitDifferentialDecoder()
DifferentialDecoder2 = demod.InitDifferentialDecoder()

# DifferentialDecoder1 = {'LastBit':0, 'NewBit':0, 'Result':0}
# DifferentialDecoder2 = {'LastBit':0, 'NewBit':0, 'Result':0}
AX25Decoder1 = demod.InitAX25Decoder()
AX25Decoder2 = demod.InitAX25Decoder()
# AX25Decoder1 = {'NewBit':0, 'BitIndex':0, 'Ones':0, 'ByteCount':0, 'WorkingByte':np.uint16(0), 'Result':np.array([]).astype('uint16'), 'CRC':np.array([0,0]), 'PacketCount':0, 'Verbose':0, 'OutputTrigger':False, 'CRCAge':1000000, 'UniquePackets':0}
# AX25Decoder2 = {'NewBit':0, 'BitIndex':0, 'Ones':0, 'ByteCount':0, 'WorkingByte':np.uint16(0), 'Result':np.array([]).astype('uint16'), 'CRC':np.array([0,0]), 'PacketCount':0, 'Verbose':0, 'OutputTrigger':False, 'CRCAge':1000000, 'UniquePackets':0}

# print(AFSKDemodulator1['SqrtTable'])
# print(AX25Decoder2)
# sys.exit(-3)
index1 = 0
index2 = 0
index3 = 0
envelope_index = 0

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
filtered_signal_buffer = np.zeros(round(len(audio) / FilterDecimator['DecimationRate']) + 1)
demod_sig_buffer1 = np.zeros(round(len(audio) / FilterDecimator['DecimationRate']) + 1)
demod_sig_buffer2 = np.zeros(round(len(audio) / FilterDecimator['DecimationRate']) + 1)
# mark_sig_buffer1 = np.zeros(round(len(audio) / FilterDecimator['DecimationRate']) + 1)
# space_sig_buffer1 = np.zeros(round(len(audio) / FilterDecimator['DecimationRate']) + 1)
# mark_sig_buffer2 = np.zeros(round(len(audio) / FilterDecimator['DecimationRate']) + 1)
# space_sig_buffer2 = np.zeros(round(len(audio) / FilterDecimator['DecimationRate']) + 1)

chop_audio_buffer = np.array([])
chop_filtered_audio_buffer = np.array([])
chop_demodulated_audio_buffer1 = np.array([])
chop_demodulated_audio_buffer2 = np.array([])
for sample in audio:
	if AFSKDemodulator1['ChopAudio'] == True:
		chop_audio_buffer = np.append(chop_audio_buffer, np.array([sample]))
	index1 = index1 + 1
	index2 = index2 + 1
	if index2 > len(audio) / 100:
		index2 = 0
		index3 = index3 + 1
		midpoint = DataSlicer1['Midpoint']
		#print(index3, InputPeakDetector['Envelope'], space_sig_ratio, space_sig_gain_error)
		print(f'{index3}')
		print(AFSKDemodulator1['MarkClip'], AFSKDemodulator1['SpaceClip'], AFSKDemodulator2['MarkClip'], AFSKDemodulator2['SpaceClip'])
		AFSKDemodulator1['MarkClip'] = False
		AFSKDemodulator1['SpaceClip'] = False
		AFSKDemodulator2['MarkClip'] = False
		AFSKDemodulator2['SpaceClip'] = False
	FilterDecimator['NewSample'] = sample
	FilterDecimator = demod.FilterDecimate(FilterDecimator)
	for filtered_signal in FilterDecimator['DataBuffer']:
		filtered_signal_buffer[envelope_index] = filtered_signal
		if AFSKDemodulator1['ChopAudio'] == True:
			chop_filtered_audio_buffer = np.append(chop_filtered_audio_buffer, np.array([filtered_signal]))
		AFSKDemodulator1['NewSample'] = filtered_signal
		AFSKDemodulator1 = demod.DemodulateAFSK(AFSKDemodulator1)
		for demodulated_signal in AFSKDemodulator1['Result']:
			if AFSKDemodulator1['ChopAudio'] == True:
				chop_demodulated_audio_buffer1 = np.append(chop_demodulated_audio_buffer1, np.array([demodulated_signal]))
			demod_sig_buffer1[envelope_index] = demodulated_signal
			# mark_sig_buffer1[envelope_index] = AFSKDemodulator1['MarkSig']
			# space_sig_buffer1[envelope_index] = AFSKDemodulator1['SpaceSig']
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
						try:
							bin_file = open(dirname + filename + '.bin', '+wb')
						except:
							pass
						with bin_file:
							for byte in AX25Decoder2['Result']:
								bin_file.write(byte.astype('uint8'))
							bin_file.close()
						if AFSKDemodulator1['ChopAudio'] == True:
							scipy.io.wavfile.write(dirname+filename+'-audio.wav', FilterDecimator['InputSampleRate'], chop_audio_buffer.astype(np.int16))
							chop_audio_buffer = np.array([])
							scipy.io.wavfile.write(dirname+filename+'-filtered.wav', FilterDecimator['OutputSampleRate'], chop_filtered_audio_buffer.astype(np.int16))
							chop_filtered_audio_buffer = np.array([])
							scipy.io.wavfile.write(dirname+filename+'-demod.wav', FilterDecimator['OutputSampleRate'], chop_demodulated_audio_buffer1.astype(np.int16))
							chop_demodulated_audio_buffer1 = np.array([])
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
			if AFSKDemodulator2['ChopAudio'] == True:
				chop_demodulated_audio_buffer2 = np.append(chop_demodulated_audio_buffer2, np.array([demodulated_signal]))
			demod_sig_buffer2[envelope_index] = demodulated_signal
			# mark_sig_buffer2[envelope_index] = AFSKDemodulator2['MarkSig']
			# space_sig_buffer2[envelope_index] = AFSKDemodulator2['SpaceSig']
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
						CRC = AX25Decoder2['CRC'][0]
						decodernum = '2'
						filename = f'Packet-{total_packets}_CRC-{format(CRC,"#06x")}_decoder-{decodernum}_Index-{index1}'
						print(dirname+filename)
						try:
							bin_file = open(dirname + filename + '.bin', '+wb')
						except:
							pass
						with bin_file:
							for byte in AX25Decoder2['Result']:
								bin_file.write(byte.astype('uint8'))
							bin_file.close()
						if AFSKDemodulator2['ChopAudio'] == True:
							scipy.io.wavfile.write(dirname+filename+'-audio.wav', FilterDecimator['InputSampleRate'], chop_audio_buffer.astype(np.int16))
							chop_audio_buffer = np.array([])
							scipy.io.wavfile.write(dirname+filename+'-filtered.wav', FilterDecimator['OutputSampleRate'], chop_filtered_audio_buffer.astype(np.int16))
							chop_filtered_audio_buffer = np.array([])
							scipy.io.wavfile.write(dirname+filename+'-demod.wav', FilterDecimator['OutputSampleRate'], chop_demodulated_audio_buffer2.astype(np.int16))
							chop_demodulated_audio_buffer2 = np.array([])
					else:
						if AX25Decoder1['CRC'][0] == AX25Decoder2['CRC'][0]:
							duplicate_packets += 1
							AX25Decoder1['UniquePackets'] -= 1
							AX25Decoder2['UniquePackets'] -= 1
							print('Decoder 2 Duplicate, bit delay: ', AX25Decoder1['CRCAge'])

scipy.io.wavfile.write(dirname+"DemodSignal1.wav", FilterDecimator['OutputSampleRate'], demod_sig_buffer1.astype(np.int16))
scipy.io.wavfile.write(dirname+"DemodSignal2.wav", FilterDecimator['OutputSampleRate'], demod_sig_buffer2.astype(np.int16))
scipy.io.wavfile.write(dirname+"FilteredSignal.wav", FilterDecimator['OutputSampleRate'], filtered_signal_buffer.astype(np.int16))
# scipy.io.wavfile.write(dirname+"MarkSignal1.wav", FilterDecimator['OutputSampleRate'], mark_sig_buffer1.astype(np.int16))
# scipy.io.wavfile.write(dirname+"SpaceSignal1.wav", FilterDecimator['OutputSampleRate'], space_sig_buffer1.astype(np.int16))
# scipy.io.wavfile.write(dirname+"MarkSignal2.wav", FilterDecimator['OutputSampleRate'], mark_sig_buffer2.astype(np.int16))
# scipy.io.wavfile.write(dirname+"SpaceSignal2.wav", FilterDecimator['OutputSampleRate'], space_sig_buffer2.astype(np.int16))

# print(AFSKDemodulator1)

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
	report_file.write('\n\n# Demodulator performance:\n')
	report_file.write('\n')
	report_file.write(f'# Total packets: {total_packets}')
	report_file.write('\n')
	report_file.write(f'# Duplicate packets: {duplicate_packets}')
	report_file.write('\n')
	report_file.write(f'# Demodulator 1 unique packets: {AX25Decoder1["UniquePackets"]}')
	report_file.write('\n')
	report_file.write(f'# Demodulator 1 total packets: {AX25Decoder1["PacketCount"]}')
	report_file.write('\n')
	report_file.write(f'# Demodulator 2 unique packets: {AX25Decoder2["UniquePackets"]}')
	report_file.write('\n')
	report_file.write(f'# Demodulator 2 total packets: {AX25Decoder2["PacketCount"]}')
	report_file.write('\n')
	report_file.close()

print('total packets: ', total_packets)
print('duplicate_packets: ', duplicate_packets)
print('Decoder 1 unique packets: ', AX25Decoder1['UniquePackets'])
print('Decoder 1 total packets: ', AX25Decoder1['PacketCount'])
print('Decoder 2 unique packets: ', AX25Decoder2['UniquePackets'])
print('Decoder 2 total packets: ', AX25Decoder2['PacketCount'])
print('made new directory: ', dirname)
print('done')
