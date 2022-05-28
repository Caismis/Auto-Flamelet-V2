from flamelet_main import flamelet_main
from multipdf import multipdf
from input_reader import input_reader
from time import time
import datetime
import numpy as np
import sys

from write_data import write_data, write_remap
if __name__ == "__main__":
    print('--------------------------------------------------------')
    print('Auto Flamelet 0.01')
    print('========================================================')
    input_data = input_reader()
    if not input_data:
        print('Exiting...')
        sys.exit()
    else:
        print('Input data confirmed...')
    print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print('________________________________________________________')
    print('_________________Raw Data Generation____________________')
    print('1. Create steady flamelets set, data saved as *.npy.')
    print('2. Create progress/variable flamelets set, data saved as *.npy.')
    print('________________________________________________________')
    print('_________________Raw Data Intergration__________________')
    print('3. Resume from *.npy, do presumed pdf, SLFM.')
    print('4. Resume from *.npy, do presumed pdf, FPV.')
    print('________________________________________________________')
    print('_________________One Click______________________________')
    print('5. Generate SLFM.')
    print('6. Generate FPV.')
    print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++')

    usrinput = input('Enter the number...\n')
    start_time = time()
    if usrinput == '1':
        print('=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-')
        print('Create SLFM...')
        print('=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-')
        flamelet_main(input_data, False, True)
    elif usrinput == '2':
        print('=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-')
        print('Create FPV...')
        print('=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-')
        flamelet_main(input_data, True, True)
    elif usrinput == '3':
        misc = np.load(file='misc.npy', allow_pickle=True)
        print('=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-')
        print('Presumed-PDF postprocessing...')
        print('=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-')
        uniset = multipdf(input_data, switch=False)
        print('---PDF postprocessing complete---')
        write_data(input_data, uniset, misc)
    elif usrinput == '4':
        misc = np.load(file='misc.npy', allow_pickle=True)
        print('=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-')
        print('Presumed-PDF postprocessing...')
        print('=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-')
        uniset = multipdf(input_data, switch=True)
        write_remap(input_data, uniset, misc)
    elif usrinput == '5':
        print('=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-')
        print('Create SLFM...')
        print('=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-')
        flamelet_main(input_data, False, True)
        misc = np.load(file='misc.npy', allow_pickle=True)
        print('=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-')
        print('Presumed-PDF postprocessing...')
        print('=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-')
        uniset = multipdf(input_data, switch=False)
        print('---PDF postprocessing complete---')
        write_data(input_data, uniset, misc)
    elif usrinput == '6':
        print('=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-')
        print('Create FPV...')
        print('=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-')
        flamelet_main(input_data, True, True)
        misc = np.load(file='misc.npy', allow_pickle=True)
        print('=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-')
        print('Presumed-PDF postprocessing...')
        print('=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-')
        uniset = multipdf(input_data, switch=True)
        write_remap(input_data, uniset, misc)
    end_time  = time()
    duration = end_time - start_time
    print('Running Time:', str(datetime.timedelta(seconds=duration)))
