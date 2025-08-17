import sys
import configparser
import struct
import n9600a_nco as nco

if len(sys.argv) < 2:
	print("Not enough arguments. Usage: py -3 n9600a_test.py <ini file>")
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

state['argv'] = sys.argv
state['config'] = config

nco.Test(state)