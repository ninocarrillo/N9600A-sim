import sys
import configparser
import struct
import scipy.io.wavfile
import numpy as np
import os

def GenInt16ArrayC(name, array, column_width):
	result = 'int16_t '
	result += name
	result += ' = {'
	for x in range(len(array)):
		if x % column_width == 0:
			result += ' \\\n     '
		if x == 0:
			result += f' {int(array[x])}'
		else:
			result += f', {int(array[x])}'
	result += ' };'
	return result
