import numpy as np

def ModulateFM(waveform, deviation, sample_rate):
	# create an FM waveform
	fm_waveform = np.zeros(len(waveform))
	t = 0
	for i in range(len(fm_waveform)):
	 	t += deviation * waveform[i] / sample_rate
	 	fm_waveform[i] = np.sin(2*np.pi*t)
	return fm_waveform
