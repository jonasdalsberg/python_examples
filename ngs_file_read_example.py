#!/usr/bin/env python3

### Code to read and demultiplex NGS data based on individual barcodes
### Secondary purpose was to optimise the code with respect to the runtime
### For a course on high performance computing at DTU
### Code written by Jonas Dalsberg Jørgensen in May 2022

### PSEUDOCODE ###
# load libraries
# define function to read gzipped files
# define function to extract reads for a barcode
    # function should also write the output to a file
    # file name = a combination of input file + barcode
# open ngs file + barcode file
# execute the function for each barcode (in parallel)
    # parallelism will be written into the function

## Loading libraries
import sys, time, subprocess, gzip
from joblib import Parallel, delayed
start_time2 = time.time()

## Defining functions
# Define function for opening potentially gzipped files
def openfile(filename):
    """ Opens both gzipped and non-compressed files """
    try:
        if filename.endswith('.gz'):
            fh = gzip.open(filename, mode='rb')
        else:
            fh = open(filename, 'rb')
    except IOError as err:
        print("Could not open file {}. Reason:".format(filename), err)
    return fh

def barcode_extract(barcodes):
    """ Returns a list of barcodes based on a barcode file """
    barcode_list = list()
    for line in barcodes:
        barcode = line[-9:-1].decode()
        barcode_list.append(barcode)
    return barcode_list[1:]

def demultiplex(infile, barcodes):
    """ Demultiplexes an NGS file based on a list of barcodes. """
    start_time = time.time()
    outfile_list = list()
    for barcode in barcodes:
        outfile_list.append(open("demultiplexed_files/demulti_{}_{}".format(barcode, sys.argv[1][:-3]), "w"))
    for line in infile:
        if line.startswith(b'@'):
            barcode_detected = line[-9:-1].decode()
            if barcode_detected in barcodes:
                for i in range(len(barcodes)):
                    if barcode_detected == barcodes[i]:
                        outfile = outfile_list[i]
                        correct_barcode = True
                        break
            else:
                correct_barcode = False
        if correct_barcode == True:
            line = str(line.decode())
            outfile.write(line)
    outfile.close()

## Opening the files
print("\nOpening files {} and {}. Please wait...".format(sys.argv[1], sys.argv[2]))
start_time = time.time()
if len(sys.argv) > 1:
    try:
        infile = openfile(sys.argv[1])
    except IOError as err:
        print("The program should be run like this:")
        print("python3 program.py fastqfile barcodefile")
        sys.exit(0)
    try:
        barcodes = open(sys.argv[2], "rb")
    except IOError as err:
        print("Could not open file {}. Reason:".format(sys.argv[2]), err)
        print("The program should be run like this:")
        print("python3 program.py fastqfile barcodefile")
        sys.exit(0)
time_open = time.time() - start_time
print("Files {} and {} successfully opened.\n".format(sys.argv[1], sys.argv[2]))

## Main program
# Extracting barcodes
print("Extracting barcodes... Please wait.")
start_time = time.time()
barcodes = barcode_extract(barcodes)
time_bar = time.time() - start_time
print("Extraction successful.")

print("Demultiplexing {}. Please wait...".format(sys.argv[1]))
start_time = time.time()
demultiplex(infile, barcodes)
#demultiplex(infile, barcodes[0])
#for barcode in barcodes:
#    demultiplex(infile, barcode)
#    infile.seek(0)
#result = Parallel(n_jobs=8)(delayed(demultiplex)(infile, barcode) for barcode in barcodes)
time_demulti = time.time() - start_time
print("Demultiplexing complete.")

## Closing the files
infile.close()

## RESULTS
print("\nTIME:")
print("Time spent opening files:\t\t{:.4f} seconds.".format(time_open))
print("Time spent extracting barcodes:\t\t{:.4f} seconds".format(time_bar))
print("Time spent demultiplexing:\t\t{:.4f} seconds.".format(time_demulti))
time_total = time.time() - start_time2
print("\nTime spent in total:\t\t\t{:.4f} seconds.".format(time_total))
