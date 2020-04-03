#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import re
import os.path

def read_args():
	parser = argparse.ArgumentParser(description="Recreates a macrocomplex using protein and DNA sequences,"
												+ " and the structures of protein-protein and protein-DNA "
												+ "interactions. For more --help.\n")

	parser.add_argument('-i', '--pdb',
						dest = "input",
						action = "store",
						default = None,
						help = "Directory with all the base pdb files to use.\n")

	parser.add_argument('-fa', '--fasta',
						dest = "fasta",
						action = "store",
						default = None,
						help = "Fasta file with the sequences of the complex.\n")

	parser.add_argument('-o', '--output',
						dest = "output",
						action = "store",
						default = None,
						help = "Name of the output.\n")

	parser.add_argument('-fo', '--format',
						dest = "format",
						action = "store",
						default = None,
						help = "Format of the output (pdb or mmcif).\n")

	parser.add_argument('-v', '--verbose',
						dest = "verbose",
						action = "store_true",
						default = False,
						help = "Verbose mode: print log in stderr.\n")

	parser.add_argument('-opt','--optimization',
						dest = "optimize",
						action = "store_true",
						default = False,
						help = "Optimizing the model of the macrocomplex.\n")

	parser.add_argument('-sit', '--sequence_identity_threshold',
						dest = "sequence_identity_threshold",
						action = "store",
						default = None,
						type = float,
						help = "Sequence identity threshold to use in the fasta filtering.\n")

	parser.add_argument('-r', '--resume',
    					dest="resume",
                        action="store_true",
                        default=False,
                        help="Resume the program after a crash or when using a different number of chains.\n")

	parser.add_argument('-c', '--chimera',
    					dest="open_chimera",
                        action="store_true",
                        default=False,
                        help="Open models in Chimera at the end of the execution of the program.\n")

	parser.add_argument('-gui', '--graphic_interface',
    					dest="gui",
                        action="store_true",
                        default=False,
                        help="Graphic user interface mode.\n")

	options = parser.parse_args()
	return options
