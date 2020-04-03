#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import argparse
from gooey import Gooey, GooeyParser

desc = """ This program recreates a macrocomplex using different pdb files containing
interacting protein pairs. In case you would like to have a graphic interface,
use the -gui argument in the commandline. """


@Gooey(program_name='MODELAIS')
def gui():
	"""Parse the commandline arguments and returns the namespace."""
	parser = GooeyParser(description=desc)
	parser = argparse.ArgumentParser(description="Recreates a macrocomplex using protein and DNA sequences,"
												+ " and the structures of protein-protein and protein-DNA "
												+ "interactions. The input can be the protein and DNA sequences"
												+ " or the structures of protein-protein and protein-DNA interactions."
												+ " Also, you can give a fasta file with sequences you would like to"
												+ " include in the macrocomplex. "
												+ "If you need information about the possible"
												+ " arguments, type --help.\n")

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
						dest = "outfile",
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
						help = "Optimize the macrocomplex.\n")

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


	options = parser.parse_args()
	return options
