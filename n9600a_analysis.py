import numpy as np
from scipy.fft import fft, fftfreq

def AnalyzeSpectrum(waveform, sample_rate, power_ratio):
	fft_n = len(waveform)
	time_step = 1 / sample_rate
	x = np.linspace(0.0, fft_n * time_step, fft_n, endpoint = False)
	x_fft = fftfreq(fft_n, time_step)
	waveform_fft = fft(waveform)
	fft_max = max(abs(waveform_fft))
	waveform_fft = waveform_fft / fft_max

	# create power spectrum
	waveform_psd = np.zeros(len(waveform_fft))
	for i in range(len(waveform_fft)):
		waveform_psd[i] = pow(abs(waveform_fft[i]),2)

	# calculate power bandwidth
	# first determine total
	total_power = 0
	for sample in waveform_psd:
		total_power += sample
	# now sum from center of spectrum out
	total_power = total_power / 2

	power_sum = 0
	i = -1
	while (power_sum < power_ratio) and (i < len(waveform_psd)):
		i += 1
		power_sum += waveform_psd[i] / total_power
	obw = x_fft[i] * 2

	obw_mask = np.zeros(4)
	obw_x = np.zeros(4)
	obw_mask[0] = 10*np.log10(waveform_psd[i])
	obw_mask[1] = 0
	obw_mask[2] = 0
	obw_mask[3] = 10*np.log10(waveform_psd[i])
	obw_x[0] = -x_fft[i]
	obw_x[1] = -x_fft[i]
	obw_x[2] = x_fft[i]
	obw_x[3] = x_fft[i]
	return([x_fft, 10*np.log10(waveform_psd), obw_x, obw_mask, obw])
