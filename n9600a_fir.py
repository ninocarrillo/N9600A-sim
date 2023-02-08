import struct
import numpy as np

def InitFIR(this):
	this['Buffer'] = np.zeros(len(this['Taps']))
	return this

def UpdateFIR(this, sample):
	this['Buffer'] = this['Buffer'][1:]
	this['Buffer'] = np.append(this['Buffer'], np.array([sample]))
	this['Output'] = np.rint(np.convolve(this['Buffer'], this['Taps'], 'valid') / pow(2, (16 + this['OutputShift'])))
	return this
