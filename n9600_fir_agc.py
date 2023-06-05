import numpy
import n9600a_strings as strings
import matplotlib.pyplot as plot

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

def FilterDecimate(this):
	this['FilterBuffer'] = numpy.rint(numpy.convolve(this['FilterBuffer'], this['Filter'], 'valid'))
	this['FilterBuffer'] = this['FilterBuffer'][::this['DecimationRate']]
	this['EnvelopeBuffer'] = numpy.zeros(len(this['FilterBuffer']))
	index = 0
	for data in this['FilterBuffer']:

		this['EnvelopeBuffer'][index] = this['PeakDetector']['Envelope']
		data = data // pow(2, (16 + this['FilterShift']))

		scale = this['AGCScaleTable'][int(this['PeakDetector']['Envelope'] / 128)]


		this['PeakDetector'] = PeakDetect(data, this['PeakDetector'])
		this['GainChange'] = 0
		if this['PeakDetector']['Envelope'] > this['ScalerHighThreshold']:
			this['FilterShift'] = this['FilterShift'] + 1
			this['PeakDetector']['Envelope'] = this['PeakDetector']['Envelope'] / 2
			if this['FilterShift'] > 16:
				this['FilterShift'] = 16
		if this['PeakDetector']['Envelope'] < this['ScalerLowThreshold']:
			this['FilterShift'] = this['FilterShift'] - 1
			this['PeakDetector']['Envelope'] = this['PeakDetector']['Envelope'] * 2
			if this['FilterShift'] < -16:
				this['FilterShift'] = -16

		data = (data * scale) // 32768
		this['FilterBuffer'][index] = data
		index += 1

	this['FilterBuffer'] = numpy.clip(this['FilterBuffer'], -32768, 32767)
	return this

def InitFilterDecimator(config):
	this =  {}
	this['Filter'] = strings.StringToIntArray(config['taps'])
	this['DecimationRate'] = int(config['decimation_rate'])
	this['FilterBuffer'] = numpy.zeros(len(this['Filter']))
	this['DataBuffer'] = numpy.array([])
	this['FilterShift'] = -5
	this['DecimationCounter'] = 0
	this['NewSample'] = 0
	this['SampleRate'] = int(config['SampleRate']) // this['DecimationRate']
	this['ScalerHighThreshold'] = int(config['scaler_high_threshold'])
	this['ScalerLowThreshold'] = int(config['scaler_low_threshold'])
	this['PeakDetector'] = {'AttackRate':int(config['agc_attack_rate']), 'SustainPeriod': int(config['agc_sustain_period']), 'DecayRate':int(config['agc_decay_rate']), 'SustainCount':0, 'Envelope':0}
	
	
	
	# Generate 256-step AGC gain table
	# Target signal level will be from config
	this['AGCScaleTable'] = numpy.ones(256)
	for i in range(256):
		peak = 128 * i
		try:
			scale = int(config['agc_target']) / peak
		except:
			scale = 1
		if scale >= 1:
			scale = 32767
		else:
			scale = numpy.rint(scale * 32768)
		this['AGCScaleTable'][i] = scale
	return this


def doit(config):
	this = InitFilterDecimator(config)
	this['FilterBuffer'] = config['I_Audio']
	this = FilterDecimate(this)
	if bool(config['plots']) == True:
		plot.figure()
		plot.plot(this['EnvelopeBuffer'])
		plot.plot(this['FilterBuffer'])
		plot.title('Input Filter AGC')
		plot.legend(['Envelope','Filtered Signal'])
		plot.show()
	return numpy.rint(this['FilterBuffer'])
