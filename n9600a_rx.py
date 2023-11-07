import sys
import configparser
import struct
import n9600a_dpsk3 as dpsk3
import n9600a_dpsk4 as dpsk4
import n9600a_rx_afsk as afsk
import n9600a_rx_dpsk2 as dpsk2
import n9600a_rx_dpsk as dpsk
import n9600a_rx_gfsk as gfsk
import n9600a_rx_gfsk2 as gfsk2
import n9600a_rx_rrcfsk as rrcfsk
import n9600a_nco as nco
import n9600a_filters as filters
import n9600a_bpsk as bpsk
import n9600a_qpsk as qpsk
import n9600a_qpsk32 as qpsk32
import n9600a_bpsk32 as bpsk32

if len(sys.argv) < 3:
	print("Not enough arguments. Usage: py -3 n9600a_rx.py <ini file> <wav file>")
	sys.exit(-1)
state = {}

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


key_string = "reports"
try:
	state[f'{key_string}'] = config[f'General'].getboolean(f'{key_string}')
except:
	print(f'{sys.argv[1]} [General] \'{key_string}\' is missing or invalid')
	sys.exit(-2)

key_string = "plots"
try:
	state[f'{key_string}'] = config[f'General'].getboolean(f'{key_string}')
except:
	print(f'{sys.argv[1]} [General] \'{key_string}\' is missing or invalid')
	sys.exit(-2)


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

elif DemodulatorType == 'gfsk2':
	gfsk2.FullProcess(state)

elif DemodulatorType == 'ncotest':
	nco.Test(state)

elif DemodulatorType == 'dpsk3':
	dpsk3.FullProcess(state)

elif DemodulatorType == 'dpsk4':
	dpsk4.FullProcess(state)

elif DemodulatorType == 'iirtest':
	filters.IIRTest(state)

elif DemodulatorType == 'rrcfsk':
	rrcfsk.FullProcess(state)

elif DemodulatorType== 'bpsk':
	bpsk.FullProcess(state)

elif DemodulatorType == 'qpsk':
	qpsk.FullProcess(state)

elif DemodulatorType == 'qpsk32':
	qpsk32.FullProcess(state)

elif DemodulatorType == 'bpsk32':
	bpsk32.FullProcess(state)
