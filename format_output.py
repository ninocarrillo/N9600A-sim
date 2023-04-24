import sys
import configparser
import struct
import scipy.io.wavfile
import numpy as np
import os

def GenInt16ArrayC(name, array, column_width):
	result = '\n'
	result += f'const __prog__ int16_t __attribute__((space(prog))) {name}[{len(array)}] = '
	result += '{ '
	y = len(array)
	for x in range(y):
		if x % column_width == 0:
			result += ' \\\n     '
		result += f' {int(array[x])}'
		if x < (y-1):
			result += ','
	result += ' };'
	return result
