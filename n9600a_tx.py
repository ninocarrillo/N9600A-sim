import sys
import configparser
import struct
import n9600a_dpsk3 as dpsk3
import n9600a_dpsk4 as dpsk4
import n9600a_rx_afsk as afsk
import n9600a_rx_dpsk2 as dpsk2
import n9600a_rx_dpsk as dpsk
import n9600a_rx_gfsk as gfsk
import n9600a_shaped4fsk as shaped4fsk
import n9600a_nco as nco
import n9600a_filters as filters

if len(sys.argv) < 3:
	print("Not enough arguments. Usage: py -3 n9600a_tx.py <ini file> <input bin file> <output wav file>")
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

if ModulatorType == 'afsk':
	afsk.Modulate(state)

elif ModulatorType == 'dpsk':
	dpsk.Modulate(state)
	
elif ModulatorType == 'gfsk':
	gfsk.Modulate(state)

elif ModulatorType == 'rrc4fsk':
	print('starting RRC4FSK modulator')
	shaped4fsk.ModulateRRC(state)
