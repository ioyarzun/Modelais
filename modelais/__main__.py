#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import arguments as arguments
import interface as interface
import optimize as optimize
import utils as utils
import copy
import random
import messages as msg
import os

# Step 1: parse the arguments
options = arguments.read_args()
if options.gui:
    options = interface.gui()

files = utils.get_PDB_input_files(options.input)

if options.verbose:
    utils.options = options
    msg.write_welcoming(files)

# Step 2: get possible structures for macrocomplex construction and skip others
if options.resume:
    (chains, pairs, similar_chains, structures) = utils.resume(options)
else:
    (chains, pairs, similar_chains, structures) = utils.get_information(files, options)


if options.fasta:
    fasta_filename=options.fasta

    if options.verbose:
        sys.stderr.write("\n Fasta sequenceces id parsed:\n")
        for identifier, fasta_sequence in utils.FASTA_iterator(fasta_filename):
            sys.stderr.write("  - "+identifier+"\n")
        sys.stderr.write("\n Fasta filtering in progress ...")

    sequence_identity_threshold=0.30 #Default value >>>>>>>>>>
    if options.sequence_identity_threshold:
        sequence_identity_threshold=float(options.sequence_identity_threshold)

    (chains, pairs, similar_chains, structures) = utils.fasta_filtering(fasta_filename, chains, pairs, similar_chains, structures, sequence_identity_threshold,options)

# Step 3: Macrocomplex reconstruction
complexes_found = []

if options.verbose:
    sys.stderr.write("\n# Beginning to construct the complex\n\n")

def construct_complex(current_complex_real, similar_chains,
                      structures, used_pairs_real=[], clashing_real=[],
                      old_complex_real=[]):
    # Making copies to avoiding modifying objects concurrently
    current_complex = copy.deepcopy(current_complex_real)
    used_pairs = copy.deepcopy(used_pairs_real)
    clashing = copy.deepcopy(clashing_real)
    old_complex = copy.deepcopy(old_complex_real)

    # FIRST ROUND: random pair of chains
    if current_complex is None:
        random_choice_id = random.choice(list(structures.keys()))
        random_choice = structures[random_choice_id]
        used_pairs.append(random_choice_id)
        if options.verbose:
            msg.beginning(random_choice_id)
        construct_complex(random_choice, similar_chains, structures, used_pairs, clashing, old_complex)
        return
    # OTHER ROUNDS
    else:
        for chain_in_cc in utils.get_chain_ids_from_structure(current_complex):
            ps = utils.get_possible_structures(chain_in_cc, similar_chains, structures, used_pairs, clashing)
            global complexes_found

            # Saving the current complex if it reached the max number of chains
            if len(ps) == 0: #>>>>>>>>> si no hay similar chains lo devuelve tal cual
                if options.verbose:
                    msg.complex_built()
                complexes_found.append(copy.deepcopy(current_complex))
                return

            for similar_chain_id in ps:
                structure_id = ps[similar_chain_id]
                structure_to_superimpose = structures[structure_id]
                other = [tuple_id for tuple_id in structure_id if tuple_id != similar_chain_id][0]
                if options.verbose:
                    msg.trying_superimpose(other, structure_id)
                matrix = utils.superimpose_chains(utils.get_chain(current_complex, chain_in_cc), utils.get_chain(structure_to_superimpose, similar_chain_id))
                matrix.apply(utils.get_chain(structure_to_superimpose, other))
                chain_to_add = utils.get_chain(structure_to_superimpose, other)
                if not utils.clashing_filter(current_complex, chain_to_add):
                    if options.verbose:
                        sys.stderr.write("Modelais is adding a chain to the complex \n")
                        sys.stderr.write("__________________________________________\n")
                    utils.add_chain(current_complex, chain_to_add)
                    used_pairs.append(structure_id)

                else:
                    clashing.append(structure_id)

                if current_complex is not None:
                    construct_complex(current_complex, similar_chains, structures, used_pairs,clashing, old_complex)
                    # if one suitable model is found the program stops
                    if len(complexes_found) == 1:
                        return

macrocomplex = construct_complex(None, similar_chains, structures)

# Step 6: write output file
format="pdb"
if options.format is not None:
    format = str(options.format)
    if format != "pdb" and format != "mmcif":
        print("Error! The output format can be either pdb or mmcif. \n")
    if format == "mmcif":
        format="cif"


dirName="results/"+str(options.output)
if not os.path.exists(dirName):
    os.mkdir(dirName)

output_files=[]
index_file = 0
for complex in complexes_found:
    for chain in complex[0]:
        chain_name=chain.id.split("_")[-1]
        chain.id=chain_name
    index_file += 1
    outname = dirName+"/"+str(options.output)+str(index_file)+"."+format
    utils.export_structure(complex, outname, format)

    # List of output files for use in optimization
    output_files.append(outname)

# Step 7 : plotting and optional optimization
for output_file in output_files:
    optimize.optimization(output_file,options)

# Step 8 (optional): open models in Chimera
if options.open_chimera:
    utils.open_in_chimera(options,dirName)
