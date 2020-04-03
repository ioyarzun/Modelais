#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

def write_welcoming(input_files): #executing_message
    sys.stderr.write("\n\n_______________________________________\n\n")
    sys.stderr.write("----  MODELAIS has been executed  ----\n")
    sys.stderr.write("_______________________________________\n\n\n")
    sys.stderr.write("Input correctly parsed.\nFiles used as input:\n")
    for file in input_files:
        sys.stderr.write("   - %s \n" % file)

def beginning(random_choice_id): #starting_pair
    sys.stderr.write("Randomly selected starting pair: %s \n" % str(random_choice_id))

def trying_superimpose(other, structure_id): #attemp_add
    sys.stderr.write("\nAttempting to add: %s \n" % str(other))
    sys.stderr.write("Superimposing structure: %s \n" % str(structure_id))

def complex_built():
    sys.stderr.write("\nThe complex was built!\n")
    sys.stderr.write("In case you would like to run modelais again with other parameters,"
                     + "remember to use --resume or -r to avoid unnecessary computation.\n\n")
