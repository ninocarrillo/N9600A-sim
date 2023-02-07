import sys
import configparser
import struct
import scipy.io.wavfile
import numpy as np
import os
import n9600a_progdemod as demod
import format_output as fo
import n9600a_strings as strings
import n9600a_input_filter as input_filter

def GetNCOConfig(config, num):
	this = {}
	try:
		this['Enabled'] = config[f'GFSK Demodulator {num}'].getboolean('enabled')
	except:
		print(f'{sys.argv[1]} [GFSK Demodulator {num}] \'enabled\' is missing or invalid')
		sys.exit(-2)
	return this

def Test(state):
	return
	
	