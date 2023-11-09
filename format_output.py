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
		result += f' {int(np.rint(array[x]))}'
		if x < (y-1):
			result += ','
	result += ' };'
	return result

def GenInt16ArrayHexC(name, array, column_width):
	result = '\n'
	result += f'const __prog__ int16_t __attribute__((space(prog))) {name}[{len(array)}] = '
	result += '{ '
	y = len(array)
	for x in range(y):
		if x % column_width == 0:
			result += ' \\\n     '
		result += f' {hex(int(np.rint(array[x])))}'
		if x < (y-1):
			result += ','
	result += ' };'
	return result

def GenCSV(array):
	result = ''
	column_count = len(array)
	row_count = len(array[0])
	for row in range(row_count):
		result += '\n'
		for column in range(column_count):
			if column > 0:
				result += f', '
			result += f'{array[column][row]}'
	return result
