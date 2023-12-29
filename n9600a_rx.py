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
import numpy as np
import format_output as fo

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

key_string = "iterators"
try:
	state[f'{key_string}'] = int(config[f'General'][f'{key_string}'])
except:
	print(f'{sys.argv[1]} [General] \'{key_string}\' is missing or invalid')
	sys.exit(-2)


state['argv'] = sys.argv
state['config'] = config

total_cycles = 1

if state['iterators'] > 0:
	print('Iterators Active.')
	IteratorStructure = []
	for iterator_number in range(state['iterators'] + 1):
		print(f'Reading settings for Iterator {iterator_number}')
		IteratorStructure.append({})
		if config.has_section(f'Iterator {iterator_number}'):
			key_string = "section"
			try:
				IteratorStructure[iterator_number][f'{key_string}'] = config[f'Iterator {iterator_number}'][f'{key_string}']
			except:
					print(f'{sys.argv[1]} [f\'Iterator {iterator_number}\'] \'{key_string}\' is missing or invalid')
					sys.exit(-2)
			key_string = "key name"
			try:
				IteratorStructure[iterator_number][f'{key_string}'] = config[f'Iterator {iterator_number}'][f'{key_string}']
			except:
					print(f'{sys.argv[1]} [f\'Iterator {iterator_number}\'] \'{key_string}\' is missing or invalid')
					sys.exit(-2)
			key_string = "low value"
			try:
				IteratorStructure[iterator_number][f'{key_string}'] = float(config[f'Iterator {iterator_number}'][f'{key_string}'])
			except:
					print(f'{sys.argv[1]} [f\'Iterator {iterator_number}\'] \'{key_string}\' is missing or invalid')
					sys.exit(-2)
			key_string = "high value"
			try:
				IteratorStructure[iterator_number][f'{key_string}'] = float(config[f'Iterator {iterator_number}'][f'{key_string}'])
			except:
					print(f'{sys.argv[1]} [f\'Iterator {iterator_number}\'] \'{key_string}\' is missing or invalid')
					sys.exit(-2)
			key_string = "step value"
			try:
				IteratorStructure[iterator_number][f'{key_string}'] = float(config[f'Iterator {iterator_number}'][f'{key_string}'])
			except:
					print(f'{sys.argv[1]} [f\'Iterator {iterator_number}\'] \'{key_string}\' is missing or invalid')
					sys.exit(-2)
			IteratorStructure[iterator_number]['current value'] = IteratorStructure[iterator_number]['low value']
			IteratorStructure[iterator_number]['cycles'] = np.ceil((IteratorStructure[iterator_number]['high value'] - IteratorStructure[iterator_number]['low value']) / IteratorStructure[iterator_number]['step value']) + 1
			total_cycles = total_cycles * IteratorStructure[iterator_number]['cycles']


	print(f'Total cycles for all iterations: {total_cycles}')
	# Create empty results matrix. packet count will be in the highest row index.
	IteratorResults = np.zeros((state['iterators']+1,int(total_cycles)))

	iterating = True
	cycle_number = 0
	while iterating:
		for iterator_number in range(1, state['iterators'] + 1):
			IteratorResults[iterator_number - 1][cycle_number] = IteratorStructure[iterator_number]['current value']
			config[f"{IteratorStructure[iterator_number]['section']}"][f"{IteratorStructure[iterator_number]['key name']}"] = str(IteratorStructure[iterator_number]['current value'])
			print(f"{IteratorStructure[iterator_number]['section']} {IteratorStructure[iterator_number]['key name']} {IteratorStructure[iterator_number]['current value']}")
		state['config'] = config

		if DemodulatorType == 'afsk':
			afsk.FullProcess(state)

		elif DemodulatorType == 'dpsk':
			dpsk.FullProcess(state)

		elif DemodulatorType == 'dpsk2':
			dpsk2.FullProcess(state)

		elif DemodulatorType == 'gfsk':
			IteratorResults[state['iterators']][cycle_number] = gfsk.FullProcess(state)

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
			IteratorResults[state['iterators']][cycle_number] = rrcfsk.FullProcess(state)

		elif DemodulatorType== 'bpsk':
			IteratorResults[state['iterators']][cycle_number] = bpsk.FullProcess(state)

		elif DemodulatorType == 'qpsk':
			IteratorResults[state['iterators']][cycle_number] = qpsk.FullProcess(state)

		elif DemodulatorType == 'qpsk32':
			IteratorResults[state['iterators']][cycle_number] = qpsk32.FullProcess(state)

		elif DemodulatorType == 'bpsk32':
			IteratorResults[state['iterators']][cycle_number] = bpsk32.FullProcess(state)

		#now update the appropriate iterators
		cycle_number += 1
		iterating = False
		for update_index in range(1, state['iterators'] + 1):
			IteratorStructure[update_index]['current value'] += IteratorStructure[update_index]['step value']
			if IteratorStructure[update_index]['current value'] > IteratorStructure[update_index]['high value']:
				IteratorStructure[update_index]['current value'] = IteratorStructure[update_index]['low value']
				# this iterator rolled over, update the next one (if there is a next one)

			else:
				# this iterator did not roll over, tell routine to continue and break out of this 'for'
				iterating = True
				break

	for print_index in range(1, state['iterators'] +1):
		if print_index > 1:
			print(', ', end='')
		print(IteratorStructure[print_index]['key name'], end='')
	print(', decoded packets', end='')
	result = fo.GenCSV(IteratorResults)
	print(result)


else:
	# No iterators, just run once.
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
