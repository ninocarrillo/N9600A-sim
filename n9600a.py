import sys
import configparser
import os
import scipy.io.wavfile
import numpy

import n9600_bit_shifter as bit_shifter
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
lastProcessBlock['OutputSampleRate'] = samplerate
lastProcessBlock['Name'] = 'initialize'
BlockNumber = 0
print(lastProcessBlock)
print(f'Processing demodulator objects:')
for ProcessBlock in config:
    print(f'- [{ProcessBlock}]')
    thisProcessBlock = {}
    thisProcessBlock['OutputSampleRate'] = lastProcessBlock['OutputSampleRate']
    thisProcessBlock['Name'] = ProcessBlock
    thisProcessBlock['Number'] = str(BlockNumber)
    for config_key in config[f"{ProcessBlock}"]:
        #print(config_key)
        thisProcessBlock[f"{config_key}"] = config[f"{ProcessBlock}"][f"{config_key}"]  
    if ProcessBlock == 'DEFAULT':
        pass
    elif ProcessBlock == 'bit_shifter':
        audio = bit_shifter.doit(thisProcessBlock, audio)
    elif ProcessBlock == 'write_wavfile':
        thisProcessBlock['OutputSampleRate'] = lastProcessBlock['OutputSampleRate']
        scipy.io.wavfile.write(dirname+lastProcessBlock['Number']+f"_{lastProcessBlock['Name']}"+".wav", thisProcessBlock['OutputSampleRate'], audio.astype(numpy.int16))
    elif ProcessBlock == 'rrc_filter':
        thisProcessBlock['OutputSampleRate'] = lastProcessBlock['OutputSampleRate']
        audio = rrc_filter.doit(thisProcessBlock, audio)


    lastProcessBlock = thisProcessBlock
    BlockNumber += 1

    
exit()