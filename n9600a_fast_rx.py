import sys
import configparser
import struct
import scipy.io.wavfile
import numpy as np
import os
import n9600a_progdemod as demod
import format_output as fo

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
	AFSKDemodulator1['SpaceRatio'] = float(config['AFSK Demodulator 1']['space gain'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 1] \'space gain\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator1['CorrelatorShift'] = int(config['AFSK Demodulator 1']['correlator shift'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 1] \'correlator shift\' is missing or invalid')
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
	AFSKDemodulator1['SquareSumBitCount'] = int(config['AFSK Demodulator 1']['square sum bit count'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 1] \'square sum bit count\' is missing or invalid')
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

# Read settings for AFSK Demodulator 2
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
	AFSKDemodulator2['SpaceRatio'] = float(config['AFSK Demodulator 2']['space gain'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 2] \'space gain\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator2['CorrelatorShift'] = int(config['AFSK Demodulator 2']['correlator shift'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 2] \'correlator shift\' is missing or invalid')
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
	AFSKDemodulator2['SquareSumBitCount'] = int(config['AFSK Demodulator 2']['square sum bit count'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 2] \'square sum bit count\' is missing or invalid')
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


# Read settings for AFSK Demodulator 3
AFSKDemodulator3 = {}
try:
	AFSKDemodulator3['Enabled'] = config['AFSK Demodulator 3'].getboolean('enabled')
except:
	print(f'{sys.argv[1]} [AFSKDemodulator3] \'enabled\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator3['CorrelatorTapCount'] = int(config['AFSK Demodulator 3']['correlator tap count'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 3] \'correlator tap count\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator3['MarkAmplitude'] = float(config['AFSK Demodulator 3']['mark amplitude'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 3] \'mark amplitude\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator3['SpaceRatio'] = float(config['AFSK Demodulator 3']['space gain'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 3] \'space gain\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator3['CorrelatorShift'] = int(config['AFSK Demodulator 3']['correlator shift'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 3] \'correlator shift\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator3['MarkFreq'] = int(config['AFSK Demodulator 3']['mark frequency'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 3] \'mark frequency\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator3['SpaceFreq'] = int(config['AFSK Demodulator 3']['space frequency'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 3] \'space frequency\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator3['SquareSumBitCount'] = int(config['AFSK Demodulator 3']['square sum bit count'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 3] \'square sum bit count\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator3['OutputFilter'] = StringToIntArray(config['AFSK Demodulator 3']['output filter taps'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 3] \'output filter taps\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator3['OutputFilterShift'] = int(config['AFSK Demodulator 3']['output filter shift'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 3] \'output filter shift\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator3['ChopAudio'] = config['AFSK Demodulator 3'].getboolean('chop audio enable')
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 3] \'chop audio enable\' is missing or invalid')
	sys.exit(-2)

try:
	AFSKDemodulator3['SqrtBitCount'] = int(config['AFSK Demodulator 3']['sqrt bit count'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 3] \'sqrt bit count\' is missing or invalid')
	sys.exit(-2)

DataSlicer3 = {}
try:
	DataSlicer3['BitRate'] = int(config['AFSK Demodulator 3']['slicer bit rate'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 3] \'slicer bit rate\' is missing or invalid')
	sys.exit(-2)
try:
	DataSlicer3['Rate'] = float(config['AFSK Demodulator 3']['slicer lock rate'])
except:
	print(f'{sys.argv[1]} [AFSK Demodulator 3] \'slicer lock rate\' is missing or invalid')
	sys.exit(-2)

# Initialize AFSK Demodulator 2
AFSKDemodulator3['InputSampleRate'] = FilterDecimator['OutputSampleRate']
DataSlicer3['InputSampleRate'] = FilterDecimator['OutputSampleRate']
AFSKDemodulator3 = demod.InitAFSKDemod(AFSKDemodulator3)
DataSlicer3 = demod.InitDataSlicer(DataSlicer3)

try:
	samplerate, audio = scipy.io.wavfile.read(sys.argv[2])
	# Take two bits of resolution away
	audio = audio >> (16 - FilterDecimator['InputBitCount'])
except:
	print('Unable to open wave file.')
	sys.exit(-2)

print("Opened file. \r\nSample rate:", samplerate, "\r\nLength:", len(audio))

DifferentialDecoder1 = demod.InitDifferentialDecoder()
DifferentialDecoder2 = demod.InitDifferentialDecoder()
DifferentialDecoder3 = demod.InitDifferentialDecoder()
AX25Decoder1 = demod.InitAX25Decoder()
AX25Decoder2 = demod.InitAX25Decoder()
AX25Decoder3 = demod.InitAX25Decoder()

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
demod_sig_buffer3 = np.zeros(round(len(audio) / FilterDecimator['DecimationRate']) + 1)

print(f'\nFiltering and decimating audio. ')
FilterDecimator['FilterBuffer'] = audio
FilterDecimator = demod.FilterDecimate(FilterDecimator)
print(f'Done.')

scipy.io.wavfile.write(dirname+"FilteredSignal.wav", FilterDecimator['OutputSampleRate'], FilterDecimator['FilterBuffer'].astype(np.int16))

print(f'\nDemodulating audio. ')
AFSKDemodulator1['CorrelatorBuffer'] = FilterDecimator['FilterBuffer']
AFSKDemodulator2['CorrelatorBuffer'] = FilterDecimator['FilterBuffer']
AFSKDemodulator3['CorrelatorBuffer'] = FilterDecimator['FilterBuffer']
AFSKDemodulator1 = demod.DemodulateAFSK(AFSKDemodulator1)
AFSKDemodulator2 = demod.DemodulateAFSK(AFSKDemodulator2)
AFSKDemodulator3 = demod.DemodulateAFSK(AFSKDemodulator3)
print(f'Done.')

print(f'\nSlicing, differential decoding, and AX25 decoding data. ')
loop_count = np.min([len(AFSKDemodulator1['Result']), len(AFSKDemodulator2['Result']), len(AFSKDemodulator3['Result'])])
for index in range(loop_count):
	trials = [DataSlicer1['AvgPhaseError'], DataSlicer2['AvgPhaseError'], DataSlicer3['AvgPhaseError']]
	winning_value = 100
	winning_index = 0
	losing_value = 0
	losing_index = 0
	for test in range(3):
		# if trials[test] > winning_value:
		# 	winning_index = test
		# 	winning_value = trials[test]
		if trials[test] > losing_value:
			losing_index = test
			losing_value = trials[test]
	losing_index += 1
	DataSlicer1['NewSample'] = AFSKDemodulator1['Result'][index]
	DataSlicer1 = demod.ProgSliceData(DataSlicer1)
	for data_bit in DataSlicer1['Result']:
		DifferentialDecoder1['NewBit'] = data_bit
		DifferentialDecoder1 = demod.ProgDifferentialDecode(DifferentialDecoder1)
		AX25Decoder1['NewBit'] = DifferentialDecoder1['Result']
		# if losing_index != 1:
		AX25Decoder1 = demod.ProgDecodeAX25(AX25Decoder1)
		if AX25Decoder1['OutputTrigger'] == True:
			AX25Decoder1['OutputTrigger'] = False
			# Check for unioqueness
			if ((AX25Decoder2['CRCAge'] > 10) and (AX25Decoder3['CRCAge'] > 10)) or ((AX25Decoder1['CRC'][0] != AX25Decoder2['CRC'][0]) and (AX25Decoder1['CRC'] != AX25Decoder3['CRC'])):
				total_packets += 1
				CRC = AX25Decoder1['CRC'][0]
				decodernum = '1'
				filename = f'Packet-{total_packets}_CRC-{format(CRC,"#06x")}_decoder-{decodernum}_Index-{index1}'
				print(f'{dirname+filename} {DataSlicer1["AvgPhaseError"]} {DataSlicer2["AvgPhaseError"]} {DataSlicer3["AvgPhaseError"]} ')
				# try:
				# 	bin_file = open(dirname + filename + '.bin', '+wb')
				# except:
				# 	pass
				# with bin_file:
				# 	for byte in AX25Decoder2['Output']:
				# 		bin_file.write(byte.astype('uint8'))
				# 	bin_file.close()
			else:
				if (AX25Decoder1['CRC'][0] == AX25Decoder2['CRC'][0]) and AX25Decoder2['CRCAge'] <= 10:
					duplicate_packets += 1
					if AX25Decoder1['UniqueFlag'] == True:
						AX25Decoder1['UniquePackets'] -= 1
						AX25Decoder1['UniqueFlag'] = False
					if AX25Decoder2['UniqueFlag'] == True:
						AX25Decoder2['UniquePackets'] -= 1
						AX25Decoder2['UniqueFlag'] = False
					print('Decoder 1 Duplicate 2, bit delay: ', AX25Decoder2['CRCAge'])
				if (AX25Decoder1['CRC'][0] == AX25Decoder3['CRC'][0]) and AX25Decoder3['CRCAge'] <= 10:
					duplicate_packets += 1
					if AX25Decoder1['UniqueFlag'] == True:
						AX25Decoder1['UniquePackets'] -= 1
						AX25Decoder1['UniqueFlag'] = False
					if AX25Decoder3['UniqueFlag'] == True:
						AX25Decoder3['UniquePackets'] -= 1
						AX25Decoder3['UniqueFlag'] = False
					print('Decoder 1 Duplicate 3, bit delay: ', AX25Decoder3['CRCAge'])

	DataSlicer2['NewSample'] = AFSKDemodulator2['Result'][index]
	DataSlicer2 = demod.ProgSliceData(DataSlicer2)
	for data_bit in DataSlicer2['Result']:
		DifferentialDecoder2['NewBit'] = data_bit
		DifferentialDecoder2 = demod.ProgDifferentialDecode(DifferentialDecoder2)
		AX25Decoder2['NewBit'] = DifferentialDecoder2['Result']
		# if losing_index != 2:
		AX25Decoder2 = demod.ProgDecodeAX25(AX25Decoder2)
		if AX25Decoder2['OutputTrigger'] == True:
			AX25Decoder2['OutputTrigger'] = False
			# Check for uniqueness
			if ((AX25Decoder1['CRCAge'] > 10) and (AX25Decoder3['CRCAge'] > 10)) or ((AX25Decoder2['CRC'][0] != AX25Decoder1['CRC'][0]) and (AX25Decoder2['CRC'] != AX25Decoder3['CRC'])):
				total_packets += 1
				CRC = AX25Decoder2['CRC'][0]
				decodernum = '2'
				filename = f'Packet-{total_packets}_CRC-{format(CRC,"#06x")}_decoder-{decodernum}_Index-{index1}'
				print(f'{dirname+filename} {DataSlicer1["AvgPhaseError"]} {DataSlicer2["AvgPhaseError"]} {DataSlicer3["AvgPhaseError"]} ')
				# try:
				# 	bin_file = open(dirname + filename + '.bin', '+wb')
				# except:
				# 	pass
				# with bin_file:
				# 	for byte in AX25Decoder2['Output']:
				# 		bin_file.write(byte.astype('uint8'))
				# 	bin_file.close()
			else:
				if (AX25Decoder2['CRC'][0] == AX25Decoder1['CRC'][0]) and AX25Decoder1['CRCAge'] <= 10:
					duplicate_packets += 1
					if AX25Decoder2['UniqueFlag'] == True:
						AX25Decoder2['UniquePackets'] -= 1
						AX25Decoder2['UniqueFlag'] = False
					if AX25Decoder1['UniqueFlag'] == True:
						AX25Decoder1['UniquePackets'] -= 1
						AX25Decoder1['UniqueFlag'] = False
					print('Decoder 2 Duplicate 1, bit delay: ', AX25Decoder1['CRCAge'])
				if (AX25Decoder2['CRC'][0] == AX25Decoder3['CRC'][0]) and AX25Decoder3['CRCAge'] <= 10:
					duplicate_packets += 1
					if AX25Decoder2['UniqueFlag'] == True:
						AX25Decoder2['UniquePackets'] -= 1
						AX25Decoder2['UniqueFlag'] = False
					if AX25Decoder3['UniqueFlag'] == True:
						AX25Decoder3['UniquePackets'] -= 1
						AX25Decoder3['UniqueFlag'] = False
					print('Decoder 2 Duplicate 3, bit delay: ', AX25Decoder3['CRCAge'])

	DataSlicer3['NewSample'] = AFSKDemodulator3['Result'][index]
	DataSlicer3 = demod.ProgSliceData(DataSlicer3)
	for data_bit in DataSlicer3['Result']:
		DifferentialDecoder3['NewBit'] = data_bit
		DifferentialDecoder3 = demod.ProgDifferentialDecode(DifferentialDecoder3)
		AX25Decoder3['NewBit'] = DifferentialDecoder3['Result']
		# if losing_index != 3:
		AX25Decoder3 = demod.ProgDecodeAX25(AX25Decoder3)
		if AX25Decoder3['OutputTrigger'] == True:
			AX25Decoder3['OutputTrigger'] = False
			# Check for uniqueness
			if ((AX25Decoder1['CRCAge'] > 10) and (AX25Decoder2['CRCAge'] > 10)) or ((AX25Decoder3['CRC'][0] != AX25Decoder1['CRC'][0]) and (AX25Decoder3['CRC'] != AX25Decoder2['CRC'])):
				total_packets += 1
				CRC = AX25Decoder2['CRC'][0]
				decodernum = '3'
				filename = f'Packet-{total_packets}_CRC-{format(CRC,"#06x")}_decoder-{decodernum}_Index-{index1}'
				print(f'{dirname+filename} {DataSlicer1["AvgPhaseError"]} {DataSlicer2["AvgPhaseError"]} {DataSlicer3["AvgPhaseError"]} ')
				# try:
				# 	bin_file = open(dirname + filename + '.bin', '+wb')
				# except:
				# 	pass
				# with bin_file:
				# 	for byte in AX25Decoder2['Output']:
				# 		bin_file.write(byte.astype('uint8'))
				# 	bin_file.close()
			else:
				if (AX25Decoder3['CRC'][0] == AX25Decoder1['CRC'][0]) and AX25Decoder1['CRCAge'] <= 10:
					duplicate_packets += 1
					if AX25Decoder3['UniqueFlag'] == True:
						AX25Decoder3['UniquePackets'] -= 1
						AX25Decoder3['UniqueFlag'] = False
					if AX25Decoder1['UniqueFlag'] == True:
						AX25Decoder1['UniquePackets'] -= 1
						AX25Decoder1['UniqueFlag'] = False
					print('Decoder 3 Duplicate 1, bit delay: ', AX25Decoder1['CRCAge'])
				if (AX25Decoder3['CRC'][0] == AX25Decoder2['CRC'][0]) and AX25Decoder2['CRCAge'] <= 10:
					duplicate_packets += 1
					if AX25Decoder3['UniqueFlag'] == True:
						AX25Decoder3['UniquePackets'] -= 1
						AX25Decoder3['UniqueFlag'] = False
					if AX25Decoder2['UniqueFlag'] == True:
						AX25Decoder2['UniquePackets'] -= 1
						AX25Decoder2['UniqueFlag'] = False
					print('Decoder 3 Duplicate 2, bit delay: ', AX25Decoder2['CRCAge'])

# scipy.io.wavfile.write(dirname+"DemodSignal1.wav", FilterDecimator['OutputSampleRate'], demod_sig_buffer1.astype(np.int16))
# scipy.io.wavfile.write(dirname+"DemodSignal2.wav", FilterDecimator['OutputSampleRate'], demod_sig_buffer2.astype(np.int16))
# scipy.io.wavfile.write(dirname+"FilteredSignal.wav", FilterDecimator['OutputSampleRate'], filtered_signal_buffer.astype(np.int16))

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

	report_file.write(fo.GenInt16ArrayC('MarkCorCos1', AFSKDemodulator1['MarkCOS'], 32))
	report_file.write(fo.GenInt16ArrayC('MarkCorSin1', AFSKDemodulator1['MarkSIN'], 32))
	report_file.write(fo.GenInt16ArrayC('SpaceCorCos1', AFSKDemodulator1['SpaceCOS'], 32))
	report_file.write(fo.GenInt16ArrayC('SpaceCorSin1', AFSKDemodulator1['SpaceSIN'], 32))
	report_file.write('\n')
	report_file.write(fo.GenInt16ArrayC('MarkCorCos2', AFSKDemodulator2['MarkCOS'], 32))
	report_file.write(fo.GenInt16ArrayC('MarkCorSin2', AFSKDemodulator2['MarkSIN'], 32))
	report_file.write(fo.GenInt16ArrayC('SpaceCorCos2', AFSKDemodulator2['SpaceCOS'], 32))
	report_file.write(fo.GenInt16ArrayC('SpaceCorSin2', AFSKDemodulator2['SpaceSIN'], 32))
	report_file.write('\n')
	report_file.write(fo.GenInt16ArrayC('MarkCorCos3', AFSKDemodulator3['MarkCOS'], 32))
	report_file.write(fo.GenInt16ArrayC('MarkCorSin3', AFSKDemodulator3['MarkSIN'], 32))
	report_file.write(fo.GenInt16ArrayC('SpaceCorCos3', AFSKDemodulator3['SpaceCOS'], 32))
	report_file.write(fo.GenInt16ArrayC('SpaceCorSin3', AFSKDemodulator3['SpaceSIN'], 32))

	report_file.write('\nMaximum correlation values:')
	tstep = 1.0 / AFSKDemodulator1['InputSampleRate']
	time = np.arange(0, tstep * AFSKDemodulator1['CorrelatorTapCount'], tstep)
	max_mark_cos_input = np.rint(24576 * (np.cos(2 * AFSKDemodulator1['MarkFreq'] * np.pi * time)))
	max_mark_sin_input = np.rint(24576 * (np.sin(2 * AFSKDemodulator1['MarkFreq'] * np.pi * time)))
	max_space_cos_input = np.rint(24576 * (np.cos(2 * AFSKDemodulator1['SpaceFreq'] * np.pi * time)))
	max_space_sin_input = np.rint(24576 * (np.sin(2 * AFSKDemodulator1['SpaceFreq'] * np.pi * time)))

	mark_sin_1 = np.convolve(max_mark_sin_input, AFSKDemodulator1["MarkSIN"], "valid") / pow(2, (16 + AFSKDemodulator1['CorrelatorShift']))
	mark_cos_1 = np.convolve(max_mark_sin_input, AFSKDemodulator1["MarkCOS"], "valid") / pow(2, (16 + AFSKDemodulator1['CorrelatorShift']))
	mark_square_sum_1 = mark_sin_1**2 + mark_cos_1**2

	space_sin_1 = np.convolve(max_space_sin_input, AFSKDemodulator1["SpaceSIN"], "valid") / pow(2, (16 + AFSKDemodulator1['CorrelatorShift']))
	space_cos_1 = np.convolve(max_space_sin_input, AFSKDemodulator1["SpaceCOS"], "valid") / pow(2, (16 + AFSKDemodulator1['CorrelatorShift']))
	space_square_sum_1 = space_sin_1**2 + space_cos_1**2

	mark_sin_2 = np.convolve(max_mark_sin_input, AFSKDemodulator2["MarkSIN"], "valid") / pow(2, (16 + AFSKDemodulator2['CorrelatorShift']))
	mark_cos_2 = np.convolve(max_mark_sin_input, AFSKDemodulator2["MarkCOS"], "valid") / pow(2, (16 + AFSKDemodulator2['CorrelatorShift']))
	mark_square_sum_2 = mark_sin_2**2 + mark_cos_2**2

	space_sin_2 = np.convolve(max_space_sin_input, AFSKDemodulator2["SpaceSIN"], "valid") / pow(2, (16 + AFSKDemodulator2['CorrelatorShift']))
	space_cos_2 = np.convolve(max_space_sin_input, AFSKDemodulator2["SpaceCOS"], "valid") / pow(2, (16 + AFSKDemodulator2['CorrelatorShift']))
	space_square_sum_2 = space_sin_2**2 + space_cos_2**2

	mark_sin_3 = np.convolve(max_mark_sin_input, AFSKDemodulator3["MarkSIN"], "valid") / pow(2, (16 + AFSKDemodulator2['CorrelatorShift']))
	mark_cos_3 = np.convolve(max_mark_sin_input, AFSKDemodulator3["MarkCOS"], "valid") / pow(2, (16 + AFSKDemodulator2['CorrelatorShift']))
	mark_square_sum_3 = mark_sin_2**2 + mark_cos_2**2

	space_sin_3 = np.convolve(max_space_sin_input, AFSKDemodulator3["SpaceSIN"], "valid") / pow(2, (16 + AFSKDemodulator3['CorrelatorShift']))
	space_cos_3 = np.convolve(max_space_sin_input, AFSKDemodulator3["SpaceCOS"], "valid") / pow(2, (16 + AFSKDemodulator3['CorrelatorShift']))
	space_square_sum_3 = space_sin_3**2 + space_cos_3**2

	report_file.write('\n')
	report_file.write(f'\nMark Sin Correlator 1: {int(np.rint(mark_sin_1[0]))}')
	report_file.write(f'\nMark Cos Correlator 1: {int(np.rint(mark_cos_1[0]))}')
	report_file.write(f'\nMark Square Sum Correlator 1: {int(np.rint(mark_square_sum_1[0]))}')
	report_file.write(f'\nMark Square Sum / Square Scale: {int(np.rint(mark_square_sum_1[0] / AFSKDemodulator1["SquareScale"]))}')
	# report_file.write(f'\n Log2(Mark Square Sum / Sqrt Table Size): {np.log2(mark_square_sum_1[0] // 2**AFSKDemodulator1["SqrtBitCount"]):.2f}')

	report_file.write('\n')
	report_file.write(f'\nSpace Sin Correlator 1: {int(np.rint(space_sin_1[0]))}')
	report_file.write(f'\nSpace Cos Correlator 1: {int(np.rint(space_cos_1[0]))}')
	report_file.write(f'\nSpace Square Sum Correlator 1: {int(np.rint(space_square_sum_1[0]))}')
	report_file.write(f'\nMark Square Sum / Square Scale: {int(np.rint(space_square_sum_1[0] / AFSKDemodulator1["SquareScale"]))}')
	# report_file.write(f'\n Log2(Space Square Sum / Sqrt Table Size): {np.log2(space_square_sum_1[0] // 2**AFSKDemodulator1["SqrtBitCount"]):.2f}')

	report_file.write('\n')
	report_file.write(f'\nMark Sin Correlator 2: {int(np.rint(mark_sin_2[0]))}')
	report_file.write(f'\nMark Cos Correlator 2: {int(np.rint(mark_cos_2[0]))}')
	report_file.write(f'\nMark Square Sum Correlator 2: {int(np.rint(mark_square_sum_2[0]))}')
	report_file.write(f'\nMark Square Sum / Square Scale: {int(np.rint(mark_square_sum_2[0] / AFSKDemodulator2["SquareScale"]))}')
	# report_file.write(f'\n Log2(Mark Square Sum / Sqrt Table Size): {np.log2(mark_square_sum_2[0] // 2**AFSKDemodulator2["SqrtBitCount"]):.2f}')

	report_file.write('\n')
	report_file.write(f'\nSpace Sin Correlator 2: {int(np.rint(space_sin_2[0]))}')
	report_file.write(f'\nSpace Cos Correlator 2: {int(np.rint(space_cos_2[0]))}')
	report_file.write(f'\nSpace Square Sum Correlator 2: {int(np.rint(space_square_sum_2[0]))}')
	report_file.write(f'\nMark Square Sum / Square Scale: {int(np.rint(space_square_sum_2[0] / AFSKDemodulator2["SquareScale"]))}')
	# report_file.write(f'\n Log2(Space Square Sum / Sqrt Table Size): {np.log2(space_square_sum_2[0] // 2**AFSKDemodulator2["SqrtBitCount"]):.2f}')


	report_file.write('\n\n')
	report_file.write(f'Square Scale 1: {AFSKDemodulator1["SquareScale"]}')

	report_file.write('\n')
	report_file.write(f'Square Scale 2: {AFSKDemodulator2["SquareScale"]}')

	report_file.write('\n\n')
	report_file.write(fo.GenInt16ArrayC(f'SquareRoot{AFSKDemodulator1["SqrtBitCount"]}', AFSKDemodulator1['SqrtTable'], 16))

	report_file.write('\n\n# Demodulator performance:\n')
	report_file.write('\n')
	report_file.write(f'# Total packets: {total_packets}')
	report_file.write('\n')
	report_file.write(f'# Duplicate packets: {duplicate_packets}')
	report_file.write('\n')
	report_file.write(f'# Demodulator 1 unique packets: {AX25Decoder1["UniquePackets"] // 2}')
	report_file.write('\n')
	report_file.write(f'# Demodulator 1 total packets: {AX25Decoder1["PacketCount"]}')
	report_file.write('\n')
	report_file.write(f'# Demodulator 2 unique packets: {AX25Decoder2["UniquePackets"] // 2}')
	report_file.write('\n')
	report_file.write(f'# Demodulator 2 total packets: {AX25Decoder2["PacketCount"]}')
	report_file.write('\n')
	report_file.write(f'# Demodulator 3 unique packets: {AX25Decoder3["UniquePackets"] // 2}')
	report_file.write('\n')
	report_file.write(f'# Demodulator 3 total packets: {AX25Decoder3["PacketCount"]}')
	report_file.write('\n')
	report_file.close()

print('total packets: ', total_packets)
print('duplicate_packets: ', duplicate_packets)
print('Decoder 1 unique packets: ', AX25Decoder1['UniquePackets'] // 2)
print('Decoder 1 total packets: ', AX25Decoder1['PacketCount'])
print('Decoder 2 unique packets: ', AX25Decoder2['UniquePackets'] // 2)
print('Decoder 2 total packets: ', AX25Decoder2['PacketCount'])
print('Decoder 3 unique packets: ', AX25Decoder3['UniquePackets'] // 2)
print('Decoder 3 total packets: ', AX25Decoder3['PacketCount'])
print('made new directory: ', dirname)
print('done')
