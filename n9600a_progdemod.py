import numpy as np
import n9600a_crc as crc

def InitAX25Decoder():
	decoder = {'NewBit':0, 'BitIndex':0, 'Ones':0, 'ByteCount':0, 'WorkingByte':np.uint16(0), 'Result':np.array([]).astype('uint16'), 'CRC':np.array([-65536,-65536]), 'PacketCount':0, 'Verbose':0, 'OutputTrigger':False, 'CRCAge':1000000, 'UniquePackets':0}
	decoder['UniqueFlag'] = False
	decoder['LastChop'] = 0

	return decoder

def InitDifferentialDecoder():
	decoder = {'LastBit':0, 'NewBit':0, 'Result':0}
	return decoder
	
def InitNCO(this):
	this.NormalizationFactor = np.rint(this.DesignSampleRate / this.WavetableSize)
	this.InPhase = 0
	this.QuadraturePhase = this.InPhase + np.rint(this.DesignSampleRate / 4)
	this.Control = 0
	this.Amplitude = pow(2,this.WavetableResolutionBits - 1) - 1
	this.WaveTable = np.zeros(this.WaveTableSize)
	for i in range(this.WavetableSize):
		this.WaveTable[i] = np.rint(this.Amplitude * sin(i * 2 * np.pi / this.WavetableSize))
	return this

def InitDescrambler():
	descrambler = {}
	descrambler['Polynomial'] = int('0x63003',16)
	descrambler['Tap'] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
	descrambler['Invert'] = True
	descrambler['TapCount'] = 0
	for i in range(32):
		if (descrambler['Polynomial'] >> i) & 1 == 1:
			descrambler['Tap'][descrambler['TapCount']] = i
			descrambler['TapCount'] += 1
	descrambler['InputMask'] = 1 << (descrambler['Tap'][descrambler['TapCount'] - 1])
	descrambler['OutputMask'] = 1 << descrambler['Tap'][1]
	descrambler['FeedbackMask'] = 1
	descrambler['BitDelay'] = descrambler['Tap'][descrambler['TapCount'] - 1] - descrambler['Tap'][1]
	descrambler['BitsInProgress'] = 0
	descrambler['Order'] = 1 << descrambler['Tap'][descrambler['TapCount'] - 1]
	descrambler['ShiftRegister'] = int('0xFFFF', 16)
	descrambler['Initialized'] = True
	return descrambler

def InitDataSlicer(data_slicer):
	data_slicer['PLLClock'] = 0
	data_slicer['PLLStep'] = 32
	data_slicer['CalculatedFeedbackRate'] = np.rint(data_slicer['Rate'] * 64)
	data_slicer['Oversample'] = data_slicer['InputSampleRate'] // data_slicer['BitRate']
	data_slicer['PLLPeriod'] = (data_slicer['InputSampleRate'] // data_slicer['BitRate']) * data_slicer['PLLStep']
	data_slicer['LastSample'] = 0
	data_slicer['NewSample'] = 0
	data_slicer['Result'] = 0
	data_slicer['Midpoint'] = 0
	data_slicer['LoosePhaseTolerance'] = 2
	data_slicer['TightPhaseTolerance'] = 1
	data_slicer['DCD'] = 0
	data_slicer['DCDLoad'] = 80
	data_slicer['Phase'] = 0
	data_slicer['PhaseError'] = 0
	data_slicer['PhaseTolerance'] = data_slicer['TightPhaseTolerance']
	data_slicer['CrossingsInSyncThreshold'] = 4
	data_slicer['CrossingsInSync'] = 0
	data_slicer['CrossingPhase'] = 0
	data_slicer['ZeroCrossing'] = False

	return data_slicer

def InitDataSlicer2(data_slicer):
	data_slicer['PLLClock'] = 0
	data_slicer['Oversample'] = data_slicer['InputSampleRate'] // data_slicer['BitRate']
	data_slicer['PLLPeriod'] = (data_slicer['InputSampleRate'] // data_slicer['BitRate']) * data_slicer['PLLStep']
	data_slicer['LastSample'] = 0
	data_slicer['NewSample'] = 0
	data_slicer['Result'] = 0
	data_slicer['Midpoint'] = 0
	data_slicer['ZeroCrossing'] = False
	data_slicer['PhaseBufferN'] = len(data_slicer['LoopFilter'])
	data_slicer['PhaseBuffer'] = np.zeros(data_slicer['PhaseBufferN'])
	data_slicer['SampleIndex'] = 0
	data_slicer['PLLControl'] = 0
	return data_slicer

def InitFilterDecimator(filter_decimator):
	filter_decimator['FilterBuffer'] = np.zeros(len(filter_decimator['Filter']))
	filter_decimator['DataBuffer'] = np.array([])
	filter_decimator['FilterShift'] = -5
	filter_decimator['DecimationCounter'] = 0
	filter_decimator['NewSample'] = 0
	filter_decimator['PeakDetector'] = {'AttackRate':filter_decimator['InputAGCAttackRate'], 'SustainPeriod': filter_decimator['InputAGCSustainPeriod'], 'DecayRate':filter_decimator['InputAGCDecayRate'], 'SustainCount':0, 'Envelope':0}
	filter_decimator['OutputSampleRate'] = filter_decimator['InputSampleRate'] // filter_decimator['DecimationRate']
	return filter_decimator

def InitGFSKDemod(demodulator):
	return demodulator

def InitDPSKDemod(demodulator):
	demodulator['CorrelatorBuffer'] = np.zeros(demodulator['AutoCorrelatorLag'])
	demodulator['OutputFilterBuffer'] = np.zeros(len(demodulator['OutputFilter']))
	demodulator['Result'] = 0
	demodulator['NewSample'] = 0
	return demodulator


def InitAFSKDemod(demodulator):
	tstep = 1.0 / demodulator['InputSampleRate']
	time = np.arange(0, tstep * demodulator['CorrelatorTapCount'], tstep)
	# demodulator['SpaceAmplitude'] = np.rint(demodulator['MarkAmplitude'] * demodulator['SpaceRatio'])
	demodulator['SpaceCOS'] = np.zeros(demodulator['CorrelatorTapCount'])
	demodulator['SpaceSIN'] = np.zeros(demodulator['CorrelatorTapCount'])
	demodulator['MarkCOS'] = np.zeros(demodulator['CorrelatorTapCount'])
	demodulator['MarkSIN'] = np.zeros(demodulator['CorrelatorTapCount'])

	demodulator['SpaceAmplitude'] = demodulator['MarkAmplitude']
	space_phase = 0
	mark_phase = 0
	freq_span = 0
	space_freq = demodulator['SpaceFreq'] - (freq_span / 2)
	mark_freq = demodulator['MarkFreq'] - (freq_span / 2)
	mark_amplitude = demodulator['MarkAmplitude']
	space_amplitude = demodulator['SpaceAmplitude']
	freq_step = freq_span / demodulator['CorrelatorTapCount']

	for index in range(demodulator['CorrelatorTapCount']):
		demodulator['SpaceCOS'][index] = np.rint(space_amplitude * np.cos(space_phase))
		demodulator['SpaceSIN'][index] = np.rint(space_amplitude * np.sin(space_phase))
		demodulator['MarkCOS'][index] = np.rint(mark_amplitude * np.cos(mark_phase))
		demodulator['MarkSIN'][index] = np.rint(mark_amplitude * np.sin(mark_phase))
		space_phase += 2 * np.pi * space_freq * tstep
		mark_phase += 2 * np.pi * mark_freq * tstep



	# demodulator['SpaceCOS'] = np.rint(demodulator['SpaceAmplitude'] * (np.cos(2 * demodulator['SpaceFreq'] * np.pi * time)))
	# demodulator['SpaceSIN'] = np.rint(demodulator['SpaceAmplitude'] * (np.sin(2 * demodulator['SpaceFreq'] * np.pi * time)))
	# demodulator['MarkCOS'] = np.rint(demodulator['MarkAmplitude'] * (np.cos(2 * demodulator['MarkFreq'] * np.pi * time)))
	# demodulator['MarkSIN'] = np.rint(demodulator['MarkAmplitude'] * (np.sin(2 * demodulator['MarkFreq'] * np.pi * time)))
	# demodulator['SquareScale'] = 2**(30 - (demodulator['SqrtBitCount'] + 2*demodulator['CorrelatorShift']))
	demodulator['SquareScale'] = 2**((demodulator['SquareSumBitCount'] - 1) - demodulator['SqrtBitCount'])
	# demodulator['SquareScale'] = 2**(31 - demodulator['SqrtBitCount'])
	# demodulator['SquareScale'] = 2**0
	demodulator['SquareCoef'] = 2**demodulator['SqrtBitCount']
	demodulator['SquareClip'] = demodulator['SquareCoef'] - 1
	demodulator['SqrtTable'] = np.zeros(demodulator['SquareCoef'])
	for N in range(demodulator['SquareCoef']):
		demodulator['SqrtTable'][N] = round(2*np.sqrt(demodulator['SquareCoef']*N))
	demodulator['CorrelatorBuffer'] = np.zeros(demodulator['CorrelatorTapCount'])
	demodulator['OutputFilterBuffer'] = np.zeros(len(demodulator['OutputFilter']))
	demodulator['Result'] = 0
	demodulator['NewSample'] = 0
	demodulator['MarkClip'] = False
	demodulator['SpaceClip'] = False
	demodulator['OutputSampleRate'] = demodulator['InputSampleRate'] // (demodulator['OutputFilterDecimationRate'] * demodulator['CorrelatorDecimationRate'])
	return demodulator

def InitAFSKSSBDemod(demodulator):
	demodulator['OutputSampleRate'] = demodulator['InputSampleRate'] // (demodulator['OutputFilterDecimationRate'] * demodulator['CorrelatorDecimationRate'])
	demodulator['CorrelatorTapCount'] = len(demodulator['MarkFilter'])
	demodulator['CorrelatorBuffer'] = np.zeros(demodulator['CorrelatorTapCount'])
	demodulator['OutputFilterBuffer'] = np.zeros(len(demodulator['OutputFilter']))
	demodulator['Result'] = 0
	return demodulator

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
	filter['FilterBuffer'] = np.rint(np.convolve(filter['FilterBuffer'], filter['Filter'], 'valid'))
	filter['FilterBuffer'] = filter['FilterBuffer'][::filter['DecimationRate']]
	index = 0
	for data in filter['FilterBuffer']:
		data = data // pow(2, (16 + filter['FilterShift']))
		filter['FilterBuffer'][index] = data
		index += 1
		if filter['InputAGCEnabled'] == True:
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
	filter['FilterBuffer'] = np.clip(filter['FilterBuffer'], -32768, 32767)
	return filter

def ProgFilterDecimate(filter):
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

def DemodulateGFSK(demodulator):
	if demodulator['Enabled'] == True:
		demodulator['Result'] = demodulator['InputBuffer']
	return demodulator

def DemodulateDPSK2(demodulator):
	if demodulator['Enabled'] == True:
		demodulator['LagBuffer'] = np.zeros(demodulator['AutoCorrelatorLag'])
		demodulator['CorrelatorBuffer'] = np.zeros(demodulator['AutoCorrelatorLag'])
		demodulator['Result'] = np.zeros(len(demodulator['InputBuffer']))
		index = 0
		for sample in demodulator['InputBuffer']:
			lag_sample = demodulator['CorrelatorBuffer'][0]
			demodulator['CorrelatorBuffer'] = demodulator['CorrelatorBuffer'][1:]
			demodulator['CorrelatorBuffer'] = np.append(demodulator['CorrelatorBuffer'], np.array([sample]))
			demodulator['LagBuffer'] = demodulator['LagBuffer'][1:]
			demodulator['LagBuffer'] = np.append(demodulator['LagBuffer'], np.array([lag_sample]))
			demodulator['Result'][index] = np.convolve(demodulator['LagBuffer'], demodulator['CorrelatorBuffer'], 'valid') // pow(2, 16 + demodulator['CorrelatorShift'])
			index += 1
		demodulator['Result'] = np.convolve(demodulator['Result'], demodulator['OutputFilter'], 'valid') // pow(2, (16 + demodulator['OutputFilterShift']))

	return demodulator

def DemodulateDPSK3(demodulator):
	if demodulator['Enabled'] == True:

		index = 0
		for sample in demodulator['InputBuffer']:
			if demodulator['InputBuffer'][index] > 0:
				demodulator['InputBuffer'][index] = 128
			else:
				demodulator['InputBuffer'][index] = -128
			index += 1
		demodulator['Result'] = np.rint(np.multiply(demodulator['InputBuffer'][:-demodulator['AutoCorrelatorLag']], demodulator['InputBuffer'][demodulator['AutoCorrelatorLag']:]))
		demodulator['Result'] = np.convolve(demodulator['Result'], demodulator['OutputFilter'], 'valid') // pow(2, (16 + demodulator['OutputFilterShift']))
	return demodulator

def DemodulateDPSK(demodulator):
	if demodulator['Enabled'] == True:
		demodulator['Result'] = np.rint(np.multiply(demodulator['InputBuffer'][:-demodulator['AutoCorrelatorLag']], demodulator['InputBuffer'][demodulator['AutoCorrelatorLag']:]) / 8192)
		demodulator['Result'] = np.convolve(demodulator['Result'], demodulator['OutputFilter'], 'valid') // pow(2, (16 + demodulator['OutputFilterShift']))
	return demodulator

def DemodulateAFSK(demodulator):
	if demodulator['Enabled'] == True:

		mark_cos_sig = np.convolve(demodulator['CorrelatorBuffer'], demodulator['MarkCOS'], 'valid') // pow(2, (16 + demodulator['CorrelatorShift']))
		mark_sin_sig = np.convolve(demodulator['CorrelatorBuffer'], demodulator['MarkSIN'], 'valid') // pow(2, (16 + demodulator['CorrelatorShift']))

		mark_sig = np.add(np.square(mark_cos_sig), np.square(mark_sin_sig))
		mark_sig += 2**(demodulator['SquareSumBitCount'] - 1)
		mark_sig %= 2**demodulator['SquareSumBitCount']
		mark_sig -= 2**(demodulator['SquareSumBitCount'] - 1)
		mark_sig = mark_sig // demodulator['SquareScale']
		mark_sig = np.clip(mark_sig, 0, demodulator['SquareClip'])

		index = 0
		for sample in mark_sig:
			# print(sample)
			mark_sig[index] = demodulator['SqrtTable'][int(sample)]
			index = index + 1

		# demodulator['MarkSig'] = mark_sig

		space_cos_sig = np.convolve(demodulator['CorrelatorBuffer'], demodulator['SpaceCOS'], 'valid') // pow(2, (16 + demodulator['CorrelatorShift']))
		space_sin_sig = np.convolve(demodulator['CorrelatorBuffer'], demodulator['SpaceSIN'], 'valid') // pow(2, (16 + demodulator['CorrelatorShift']))

		space_cos_sig = np.rint(space_cos_sig * demodulator['SpaceRatio'])
		space_sin_sig = np.rint(space_sin_sig * demodulator['SpaceRatio'])

		space_sig = np.add(np.square(space_cos_sig), np.square(space_sin_sig))
		space_sig += 2**(demodulator['SquareSumBitCount'] - 1)
		space_sig %= 2**demodulator['SquareSumBitCount']
		space_sig -= 2**(demodulator['SquareSumBitCount'] - 1)
		space_sig = space_sig // demodulator['SquareScale']
		space_sig = np.clip(space_sig, 0, demodulator['SquareClip'])

		# space_sig = demodulator['SqrtTable'][space_sig]

		index = 0
		for sample in space_sig:
			space_sig[index] = demodulator['SqrtTable'][int(sample)]
			index = index + 1
		# demodulator['SpaceSig'] = space_sig

		demodulator['OutputFilterBuffer'] = np.subtract(mark_sig, space_sig)

		demodulator['OutputFilterBuffer'] = demodulator['OutputFilterBuffer'][::demodulator['CorrelatorDecimationRate']]

		demodulator['Result'] = np.convolve(demodulator['OutputFilterBuffer'], demodulator['OutputFilter'], 'valid') // pow(2, (16 + demodulator['OutputFilterShift']))

		demodulator['Result'] = demodulator['Result'][::demodulator['OutputFilterDecimationRate']]
	return demodulator

def ProgDemodulateAFSK(demodulator):
	if demodulator['Enabled'] == True:
		demodulator['CorrelatorBuffer'] = demodulator['CorrelatorBuffer'][1:]
		demodulator['CorrelatorBuffer'] = np.append(demodulator['CorrelatorBuffer'], np.array([demodulator['NewSample']]))

		mark_cos_sig = np.rint((np.convolve(demodulator['CorrelatorBuffer'], demodulator['MarkCOS'], 'valid')) / pow(2, (16 + demodulator['CorrelatorShift'])))
		mark_sin_sig = np.rint((np.convolve(demodulator['CorrelatorBuffer'], demodulator['MarkSIN'], 'valid')) / pow(2, (16 + demodulator['CorrelatorShift'])))

		mark_sig = np.add(np.square(mark_cos_sig), np.square(mark_sin_sig))
		# mark_sig = int(np.clip(mark_sig, 0, 32767))
		mark_sig += 32768
		mark_sig %= 65536
		mark_sig -= 32768
		mark_sig = np.rint(mark_sig / demodulator['SquareScale'])
		if mark_sig > demodulator['SquareClip']:
			demodulator['MarkClip'] = True
		mark_sig = int(np.clip(mark_sig, 0, demodulator['SquareClip']))

		# mark_sig = np.rint(demodulator['SquareOutputScale'] * np.sqrt(demodulator['SquareCoef'] * mark_sig))
		mark_sig = demodulator['SqrtTable'][mark_sig]

		demodulator['MarkSig'] = mark_sig

		space_cos_sig = np.rint((np.convolve(demodulator['CorrelatorBuffer'], demodulator['SpaceCOS'], 'valid')) / pow(2, (16 + demodulator['CorrelatorShift'])))
		space_sin_sig = np.rint((np.convolve(demodulator['CorrelatorBuffer'], demodulator['SpaceSIN'], 'valid')) / pow(2, (16 + demodulator['CorrelatorShift'])))

		space_sig = np.add(np.square(space_cos_sig), np.square(space_sin_sig))
		# space_sig = int(np.clip(space_sig, 0, 32767))
		space_sig += 32768
		space_sig %= 65536
		space_sig -= 32768
		space_sig = np.rint(space_sig / demodulator['SquareScale'])
		if space_sig > demodulator['SquareClip']:
			demodulator['SpaceClip'] = True
		space_sig = int(np.clip(space_sig, 0, demodulator['SquareClip']))

		# space_sig = np.rint(demodulator['SquareOutputScale']* np.sqrt(demodulator['SquareCoef'] * space_sig))
		space_sig = demodulator['SqrtTable'][space_sig]
		#space_sig = np.rint(space_sig * demodulator['SpaceRatio'])

		demodulator['SpaceSig'] = space_sig

		demodulator['OutputFilterBuffer'] = demodulator['OutputFilterBuffer'][1:]
		demodulator['OutputFilterBuffer'] = np.append(demodulator['OutputFilterBuffer'], np.array([mark_sig - space_sig]))

		demodulator['Result'] = np.rint(np.convolve(demodulator['OutputFilterBuffer'], demodulator['OutputFilter'], 'valid') / pow(2, (16 + demodulator['OutputFilterShift'])))
		if demodulator['OffsetRemovalEnabled'] == True:
			demodulator['EnvelopeDetector'] = HighLowDetect(demodulator['Result'], demodulator['EnvelopeDetector'])
			demodulator['Result'] = demodulator['Result'] - (demodulator['EnvelopeDetector']['Midpoint'] * demodulator['OffsetRemovalRate'])
	return demodulator

def ProgUnscramble(descrambler):
	if descrambler['NewBit'] == 1:
		descrambler['ShiftRegister'] ^= descrambler['Polynomial']
	if descrambler['Invert'] == True:
		if descrambler['ShiftRegister'] & 1 == 1:
			descrambler['Result'] = 0
		else:
			descrambler['Result'] = 1
	else:
		if descrambler['ShiftRegister'] & 1 == 1:
			descrambler['Result'] = 1
		else:
			descrambler['Result'] = 0
	descrambler['ShiftRegister'] = descrambler['ShiftRegister'] >> 1
	return descrambler

def ProgDifferentialDecode(decoder):
	if decoder['NewBit'] == decoder['LastBit']:
		decoder['Result'] = 1
	else:
		decoder['Result'] = 0
	decoder['LastBit'] = decoder['NewBit']
	return decoder


def ProgDecodeAX25(decoder):
	decoder['CRCAge'] += 1
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
				decoder['CRC'] = crc.CheckCRC(decoder['Result'].astype('uint16'), len(decoder['Result']))
				if decoder['CRC'][1] == 1:
					decoder['CRCAge'] = 0
					decoder['PacketCount'] += 1
					decoder['Output'] = decoder['Result'][:-2]
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


def SliceData(slicer):
	slicer['Midpoint'] = 0
	slicer['LastSample'] = 0
	slicer['Result'] = np.zeros(int((len(slicer['Input']) / slicer['Oversample']) * 1.1))
	output_index = 0
	for slicer['NewSample'] in slicer['Input']:
		slicer['PLLClock'] += slicer['PLLStep']
		if slicer['PLLClock'] > ((slicer['PLLPeriod'] / 2.0) - 1.0):
			slicer['PLLClock'] -= slicer['PLLPeriod']
			if slicer['NewSample'] > slicer['Midpoint']:
				slicer['Result'][output_index] = 1
			else:
				slicer['Result'][output_index] = 0
			output_index += 1
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

def ProgSliceData2(slicer):
	# slicer['EnvelopeDetector'] = HighLowDetect(slicer['NewSample'], slicer['EnvelopeDetector'])
	# slicer['Midpoint'] = np.rint(slicer['EnvelopeDetector']['Midpoint'] * 0.66)
	slicer['SampleIndex'] += 1
	slicer['Midpoint'] = 0
	slicer['Result'] = np.array([])
	slicer['OutputTrigger'] = False
	slicer['PLLClock'] += (slicer['PLLStep'] + np.rint(slicer['PLLControl'] * slicer['PLLFeedbackGain']))
	if slicer['PLLClock'] > ((slicer['PLLPeriod'] // 2) - 1):
		slicer['PLLClock'] -= slicer['PLLPeriod']
		if slicer['NewSample'] > slicer['Midpoint']:
			slicer['Result'] = np.array([1])
			slicer['OutputTrigger'] = True
		else:
			slicer['Result'] = np.array([0])
			slicer['OutputTrigger'] = True
	if slicer['LastSample'] > slicer['Midpoint']:
		if slicer['NewSample'] <= slicer['Midpoint']:
			# Zero Crossing
			slicer['ZeroCrossing'] = True
			#slicer['PLLClock'] *= slicer['CalculatedFeedbackRate']
			#slicer['PLLClock'] //= 64
	else:
		if slicer['NewSample'] > slicer['Midpoint']:
			# Zero Crossing
			slicer['ZeroCrossing'] = True
			#slicer['PLLClock'] *= slicer['CalculatedFeedbackRate']
			#slicer['PLLClock'] //= 64

	if slicer['ZeroCrossing'] == True:
		slicer['ZeroCrossing'] = False

		slicer['PhaseError'] = slicer['PLLClock']

		#slicer['PLLControl'] = -slicer['PhaseError']
		slicer['PhaseBuffer'] = slicer['PhaseBuffer'][1:]
		slicer['PhaseBuffer'] = np.append(slicer['PhaseBuffer'], np.array(slicer['PhaseError']))
		slicer['PLLControl'] = -(np.convolve(slicer['PhaseBuffer'], slicer['LoopFilter'], 'valid') // pow(2, 16))
	#else:
		#slicer['PhaseBuffer'] = slicer['PhaseBuffer'][1:]
		#slicer['PhaseBuffer'] = np.append(slicer['PhaseBuffer'], np.array(slicer['PhaseError']))


	slicer['LastSample'] = slicer['NewSample']

	return slicer

def ProgSliceData(slicer):
	# slicer['EnvelopeDetector'] = HighLowDetect(slicer['NewSample'], slicer['EnvelopeDetector'])
	# slicer['Midpoint'] = np.rint(slicer['EnvelopeDetector']['Midpoint'] * 0.66)

	slicer['Midpoint'] = 0
	slicer['Result'] = np.array([])
	slicer['OutputTrigger'] = False
	slicer['PLLClock'] += slicer['PLLStep']
	if slicer['PLLClock'] > ((slicer['PLLPeriod'] // 2) - 1):
		slicer['PLLClock'] -= slicer['PLLPeriod']
		if slicer['NewSample'] > slicer['Midpoint']:
			slicer['Result'] = np.array([1])
			slicer['OutputTrigger'] = True
		else:
			slicer['Result'] = np.array([0])
			slicer['OutputTrigger'] = True
	if slicer['LastSample'] > slicer['Midpoint']:
		if slicer['NewSample'] <= slicer['Midpoint']:
			# Zero Crossing
			slicer['ZeroCrossing'] = True
			#slicer['PLLClock'] *= slicer['Rate']
			slicer['PLLClock'] *= slicer['CalculatedFeedbackRate']
			slicer['PLLClock'] //= 64
	else:
		if slicer['NewSample'] > slicer['Midpoint']:
			# Zero Crossing
			slicer['ZeroCrossing'] = True
			slicer['PLLClock'] *= slicer['CalculatedFeedbackRate']
			slicer['PLLClock'] //= 64

	if slicer['ZeroCrossing'] == True:
		slicer['ZeroCrossing'] = False
		if slicer['DCD'] < (slicer['DCDLoad'] / 2):
			slicer['PhaseTolerance'] = slicer['TightPhaseTolerance']
		slicer['PhaseError'] = abs(slicer['Phase'] - slicer['CrossingPhase'])
		if slicer['PhaseError'] < slicer['PhaseTolerance']:
			slicer['CrossingsInSync'] += 1
			if slicer['CrossingsInSync'] > slicer['CrossingsInSyncThreshold']:
				slicer['DCD'] = slicer['DCDLoad']
				slicer['PhaseTolerance'] = slicer['LoosePhaseTolerance']
		else:
			slicer['CrossingsInSync'] //= 2
		slicer['CrossingPhase'] = slicer['Phase']

	slicer['LastSample'] = slicer['NewSample']
	slicer['Phase'] += 1
	if slicer['Phase'] >= slicer['Oversample']:
		slicer['Phase'] = 0
	return slicer

def DemodulateAFSKSSB(demodulator):
	if demodulator['Enabled'] == True:
		mark_sig = abs(np.convolve(demodulator['CorrelatorBuffer'], demodulator['MarkFilter'], 'valid') // pow(2, (16 + demodulator['CorrelatorShift'])))
		space_sig = abs(np.convolve(demodulator['CorrelatorBuffer'], demodulator['SpaceFilter'], 'valid') // pow(2, (16 + demodulator['CorrelatorShift'])))
		demodulator['OutputFilterBuffer'] = np.subtract(mark_sig, space_sig)
		demodulator['OutputFilterBuffer'] = demodulator['OutputFilterBuffer'][::demodulator['CorrelatorDecimationRate']]
		demodulator['Result'] = np.convolve(demodulator['OutputFilterBuffer'], demodulator['OutputFilter'], 'valid') // pow(2, (16 + demodulator['OutputFilterShift']))
		demodulator['Result'] = demodulator['Result'][::demodulator['OutputFilterDecimationRate']]
	return demodulator
