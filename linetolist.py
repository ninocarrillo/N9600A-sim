import sys
import numpy as np

if len(sys.argv) < 2:
    print("Not enough arguments. Usage: py -3 linetolist.py <text file>")
    sys.exit(-1)

data = np.array([])
try:
    with open(sys.argv[1]) as file:
        for line in file:
            #print(line.rstrip())
            data = np.append(data, np.array(int(line.rstrip())))
except:
    print(f'Unable to process file {sys.argv[1]}')
    sys.exit(-2)

filename = 'output_'
try:
    file = open(filename + sys.argv[1], 'w+')
    file.write('[ ')
    first = True
    for element in data:
        if first == True:
            file.write(str(element))
        else:
            file.write(', ' + str(element))
        first = False
    file.write(' ]')
    file.close()
except:
    print(f'Unable to write file {filename + sys.argv[1]}')
    sys.exit(-2)
    
print(data)