import sys
import configparser
import os
import scipy.io.wavfile
import numpy

import n9600_bit_shifter as bit_shifter
import n9600_fir_agc as fir_agc
import n9600_rrc_filter as rrc_filter

if len(sys.argv) < 3:
	print("Not enough arguments. Usage: py -3 n9600a.py <ini file> <wav file>")
	sys.exit(-1)

# read demodulator description from ini file:
config = configparser.ConfigParser()
try:
	config.read(sys.argv[1])
except:
	print(f'Unable to open ini file {sys.argv[1]}')
	sys.exit(-2)
print(f'Opened ini file {sys.argv[1]}')

print(f'Opening wav file {sys.argv[2]}')
# open wav file:
try:
	samplerate, audio = scipy.io.wavfile.read(sys.argv[2])
except:
	print('Unable to open wave file.')
	sys.exit(-2)

#generate a new directory for the reports
run_number = 0
print('Trying to make a new directory')
while True:
	run_number = run_number + 1
	dirname = f'./run{run_number}/'
	try:
		os.mkdir(dirname)
	except:
		print(dirname + ' exists')
		continue
	break
print(f'Made new directory {dirname}')

print('Initializing process object.')
lastProcessBlock = {}
lastProcessBlock['I_Audio'] = audio
lastProcessBlock['Q_Audio'] = {}
lastProcessBlock['SampleRate'] = samplerate
lastProcessBlock['block_type'] = 'initialize'
BlockNumber = 0
print(lastProcessBlock)
print(f'Processing demodulator objects:')
for ProcessBlock in config:
	print(f'- [{ProcessBlock}]')
	thisProcessBlock = {}
	thisProcessBlock['SampleRate'] = lastProcessBlock['SampleRate']
	thisProcessBlock['I_Audio'] = lastProcessBlock['I_Audio']
	thisProcessBlock['Q_Audio'] = lastProcessBlock['Q_Audio']
	thisProcessBlock['Name'] = ProcessBlock
	thisProcessBlock['Number'] = str(BlockNumber)
	for config_key in config[f"{ProcessBlock}"]:
		#print(config_key)
		thisProcessBlock[f"{config_key}"] = config[f"{ProcessBlock}"][f"{config_key}"]  

	if ProcessBlock == 'DEFAULT':
		pass
	elif thisProcessBlock['block_type'] == 'bit_shifter':
		print('**bit_shifter')
		thisProcessBlock = bit_shifter.doit(thisProcessBlock)
	elif thisProcessBlock['block_type'] == 'wavfile_writer':
		print('**wavfile_writer')
		scipy.io.wavfile.write(dirname+lastProcessBlock['Number']+f"_{lastProcessBlock['Name']}"+".wav", lastProcessBlock['SampleRate'], thisProcessBlock['I_Audio'].astype(numpy.int16))
	elif thisProcessBlock['block_type'] == 'rrc_filter':
		print('**rrc_filter')
		thisProcessBlock = rrc_filter.doit(thisProcessBlock)
	elif thisProcessBlock['block_type'] == 'fir_agc':
		print('**fir_agc')
		thisProcessBlock = fir_agc.doit(thisProcessBlock)

	lastProcessBlock = thisProcessBlock
	BlockNumber += 1

exit()