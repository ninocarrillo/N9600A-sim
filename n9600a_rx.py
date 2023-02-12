import sys
import configparser
import struct
import n9600a_dpsk3 as dpsk3
import n9600a_rx_afsk as afsk
import n9600a_rx_dpsk2 as dpsk2
import n9600a_rx_dpsk as dpsk
import n9600a_rx_gfsk as gfsk
import n9600a_nco as nco
import n9600a_filters as filters

if len(sys.argv) < 3:
	print("Not enough arguments. Usage: py -3 n9600a_rx.py <ini file> <wav file>")
	sys.exit(-1)

# read demodulator description from ini file:
config = configparser.ConfigParser()
try:
	config.read(sys.argv[1])
except:
	print(f'Unable to open ini file {sys.argv[1]}')
	sys.exit(-2)
print(f'Opened ini file {sys.argv[1]}')

print(f'Checking demodulator type')
try:
	DemodulatorType = config['General']['demodulator']
except:
	print(f'{sys.argv[1]} [General] \'demodulator\' is missing or invalid')
	sys.exit(-2)
print(f'Demodulator type is {DemodulatorType}')

state = {}
state['argv'] = sys.argv
state['config'] = config

if DemodulatorType == 'afsk':
	afsk.FullProcess(state)

elif DemodulatorType == 'dpsk':
	dpsk.FullProcess(state)

elif DemodulatorType == 'dpsk2':
	dpsk2.FullProcess(state)

elif DemodulatorType == 'gfsk':
	gfsk.FullProcess(state)

elif DemodulatorType == 'ncotest':
	nco.Test(state)

elif DemodulatorType == 'dpsk3':
	dpsk3.FullProcess(state)

elif DemodulatorType == 'iirtest':
	filters.IIRTest(state)
