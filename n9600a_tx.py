import sys
import configparser
import struct
import n9600a_dpsk3 as dpsk3
import n9600a_dpsk4 as dpsk4
import n9600a_rx_afsk as afsk
import n9600a_rx_dpsk2 as dpsk2
import n9600a_rx_dpsk as dpsk
import n9600a_rx_gfsk as gfsk
import n9600a_bpsk32 as bpsk
import n9600a_shaped4fsk as shaped4fsk
import n9600a_nco as nco
import n9600a_filters as filters
import os
import numpy as np
import random

if len(sys.argv) < 2:
	print("Not enough arguments. Usage: py -3 n9600a_tx.py <ini file>")
	sys.exit(-1)

# read demodulator description from ini file:
config = configparser.ConfigParser()
try:
	config.read(sys.argv[1])
except:
	print(f'Unable to open ini file {sys.argv[1]}')
	sys.exit(-2)
print(f'Opened ini file {sys.argv[1]}')

print(f'Checking modulator type')
try:
	ModulatorType = config['General']['modulator']
except:
	print(f'{sys.argv[1]} [General] \'modulator\' is missing or invalid')
	sys.exit(-2)
print(f'Modulator type is {ModulatorType}')

state = {}
state['argv'] = sys.argv
state['config'] = config

data_count = 10000
np.random.seed(0)
state['InputData'] = np.random.randint(0,256, data_count)
for i in range(20):
	pass
	#state['InputData'][i] = 0x77
	#state['InputData'][i+20] = 0x55
# try:
# 	with open(sys.argv[2], 'rb') as f:
# 		file_size = os.path.getsize(sys.argv[2])
# 		print("File Size: ", file_size)
# 		state['InputData'] = np.zeros(file_size)
# 		index = 0
# 		while(byte := f.read(1)):
# 			#state['InputData'][index] = int.from_bytes(byte, "big")
# 			state['InputData'][index] = random.randint(0,255)
# 			index += 1
# 		f.close()
# except:
# 	print(f'Unable to open input file {sys.argv[2]}')
# 	sys.exit(-3)



if ModulatorType == 'afsk':
	afsk.Modulate(state)

elif ModulatorType == 'dpsk':
	dpsk.Modulate(state)

elif ModulatorType == 'gfsk':
	gfsk.Modulate(state)

elif ModulatorType == 'rrcfsk':
	print('starting RRCFSK modulator')
	shaped4fsk.ModulateRRC(state)

elif ModulatorType == 'gauss4fsk':
	print('starting Gauss4FSK modulator')
	shaped4fsk.ModulateGauss(state)

elif ModulatorType == 'gaussfiltergen':
	shaped4fsk.GaussFilterGen(state)

elif ModulatorType == 'bpsk':
	bpsk.Modulate(state)

elif ModulatorType == 'bpskgauss':
	bpsk.ModulateGauss(state)

elif ModulatorType == 'bpskrrc':
	bpsk.ModulateRRC(state)
