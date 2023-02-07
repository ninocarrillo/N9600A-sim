import numpy as np

def StringToIntArray(input_string):
	working_string = ''
	result = np.array([])
	for character in input_string:
		if (character >= '0' and character <= '9') or ((character == '-') or (character == '.')):
			working_string += character
		else:
			if working_string != '':
				result = np.append(result, np.array([int(working_string)]))
				working_string = ''
	return result