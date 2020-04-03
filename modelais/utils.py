#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import Bio.PDB as pdb
from Bio.PDB.PDBIO import PDBIO
from Bio import pairwise2
import itertools
import copy
import __main__ as main
import pickle

# ________________________
# _______ Get input ______
# ________________________
def get_PDB_input_files(input):
	"""The input object is a directory, therefore all the files
	terminating in ".pdb" are collected and returned in a list object, or a list of files."""
	if input[-1] == '/':
		input = input[:-1]
	input_files = []
	if os.path.isdir(input):
		for filename in os.listdir(input):
			if filename.endswith(".pdb") and "chain" not in filename:
				input_files.append(os.path.join(input, filename))
	else:
		input_files = input

	return input_files


def FASTA_iterator( fasta_filename ):
    if fasta_filename.endswith((".fa", ".fasta",".fa.txt", ".fasta.txt")):
        fd=open(fasta_filename,"r")
        seq=""
        for line in fd:
            line=line.strip("\n")
            if line.startswith(">"):
                if seq:
                    yield (identifier,seq)
                    seq=""
                identifier=line.strip(">")
            else:
                seq=seq+line
        yield (identifier,seq)
    else:
        sys.stderr.write("\nThe input parsed is not a fasta file.\n For more information use the argument --help")


# ____________________________
# __ Resume related methods __
# ____________________________
def resume(options):
    """
    The data needed is returned in a list in order to execute the program.
    This is possible because of the binary files created by pickle. This will avoid unnecessary computation.
    """
    print("Input parsed using resume")
    print(options.input)
    prefix = 'resume/' + options.input.strip('/').split('/')[-1]
    chains = pickle.load(open(prefix + "_chains.p", "rb"))
    pairs = pickle.load(open(prefix + "_pairs.p", "rb"))
    similar_chains = pickle.load(open(prefix + "_similar_chains.p", "rb"))
    structures = pickle.load(open(prefix + "_structures.p", "rb"))

    if options.verbose:
        sys.stderr.write("The program has recovered the following data from previous computations:\n")
        sys.stderr.write("\tChains: %s\n" % len(chains))
        sys.stderr.write("\tPaired chains: %s\n" % len(pairs))
        sys.stderr.write("\tSimilar Chains: %s\n" % len(similar_chains))
        sys.stderr.write("\tStructures: %s\n\n" % len(structures))

    return (chains, pairs, similar_chains, structures)


def get_information(input_files, options):
    """
    For each input file this function takes the structure to get its chains.
    Each chain has an index, the chain Id, which is unique and is inserted into a
    dictionary using as key the chain Id.
    Furthermore, the chain Ids of each structure are appended into another list called
    pairs. After doing this with each input file, the function takes a dictionary of
    similar_chains using the chains dictionary.

    The function saves averything in binary files for the user to be able to use the resume argument.
    If the argument --verbose is used by the user, the function prints on stderr.

    Returns:
    - chains dictionary: indexed chain Ids as keys, chain objects as values.
    - pairs list: with sublists containing all the chain ids of the structures of the input files.
    - similar_chains dictionary: for each chain, the chain Id is set as key and the
      similar chains as values (objects) according to a threshold argument.
    - structures dictionary: tuple of chains Ids in the structure as key, structure object as value.
    """
    chains = {}
    pairs = []
    structures = {}
    structure_index = 1
    for file in input_files:
        paired_chains = []
        structure = get_structure_from_file(file, remove_het_atoms=True)
        for chain in get_chains_from_structure(structure, remove_het_atoms=True):
            chain_id = str(structure_index)+"_"+str(chain.id)
            chains[chain_id] = chain
            chain.id = chain_id
            paired_chains.append(chain_id)
        pairs.append(paired_chains)
        structures[tuple(paired_chains)] = structure
        structure_index += 1
    similar_chains = get_similar_chains(chains)

    prefix = 'resume/' + options.input.strip('/').split('/')[-1]
    chains_backup = open(prefix + "_chains.p", "wb")
    pairs_backup = open(prefix + "_pairs.p", "wb")
    similar_chains_backup = open(prefix + "_similar_chains.p", "wb")
    structures_backup = open(prefix + "_structures.p", "wb")
    pickle.dump(chains, chains_backup)
    pickle.dump(pairs, pairs_backup)
    pickle.dump(similar_chains, similar_chains_backup)
    pickle.dump(structures, structures_backup)

    if options.verbose:
        sys.stderr.write("The analysis of the input results in:\n")
        sys.stderr.write("\tChains: %s\n" % len(chains))
        sys.stderr.write("\tPaired chains: %s\n" % len(pairs))
        sys.stderr.write("\tSimilar Chains: %s\n" % len(similar_chains))
        sys.stderr.write("\tStructures: %s\n\n" % len(structures))
        sys.stderr.write("Resume files have been saved in %s/resume\n"
                         % os.getcwd())
        sys.stderr.write("________________________________________________\n\n")

    return (chains, pairs, similar_chains, structures)


# _________________________
# ____ Utils functions ____
# _________________________
def get_chain_sequence(chain):
    """
    This function, given a chain object, returns a string with the sequence
	of the polypeptide or the nucleotide sequence of the chain.
    """
    sequence=""
    ppb = pdb.PPBuilder()
    for pp in ppb.build_peptides(chain):
        sequence=sequence+pp.get_sequence()

    if not sequence:
        for residue in chain.get_residues():
            res= residue.get_resname()[2]
            if res in 'ATGCU':
                sequence=sequence+res

    return sequence


def get_structure_from_file(input_file, remove_het_atoms=False):
    """
    Given a PDB file, this function gets the structure from the file. If
    remove_het_atoms argument is set to true, the structure is returned
    as an object without the hetero atoms.
    """
    parser = pdb.PDBParser(QUIET=True)
    structure = parser.get_structure('X', input_file)

    if remove_het_atoms:
        for model in structure:
            for chain in model:
                remove_hetero_atoms(chain)
        return structure
    else:
        return structure


def get_chains_from_structure(structure, remove_het_atoms=False):
    """
	Creates a list of chains out of a given structure. Also, if remove_het_atoms
    argument is set to true, the hetero atoms are removed from the chains.
    Returns the list of chains objects.
    """
    chains = []
    for model in structure:
        for chain in model:
            if remove_het_atoms:
                remove_hetero_atoms(chain)
                chains.append(chain)
            else:
                chains.append(chain)
    return chains


def remove_hetero_atoms(chain):
    """
    Removes the heteroatoms from a given chain.
    """
    hetero_atoms = list(filter(lambda x: x.id[0] != " ",
                              chain.get_residues()))
    for hetero_atom in hetero_atoms:
        chain.detach_child(hetero_atom.id)


def get_chain(structure, chain_id):
    """
    Gets a specific chain based on the chain_id from a given structure.
    Returns the chain with that chain_id as an object.
    """
    for chain in get_chains_from_structure(structure):
        if chain.id == chain_id:
            return chain


def get_chain_ids_from_structure(structure):
    """
    Extracts all chain ids from a given sructure and joins them into a list.
    Returns the list of ids as strings.
    """
    chains_ids = []
    for model in structure:
        for chain in model:
            chains_ids.append(chain.id)
    return chains_ids


def add_chain(structure, chain_to_add):
    """
    Appends a given chain into a structure object modifying the chain Id if it's necessary.
    """
    chain = copy.deepcopy(chain_to_add)
    chains_ids = get_chain_ids_from_structure(structure)
    while chain.id in chains_ids:
        chain.id = chain.id+"_"+chain.id
    structure[0].add(chain)


# _____________________________________
# __ Macrocomplex building functions __
# _____________________________________
def get_similar_chains(chains, sequence_identity_threshold=0.95):
    """
    Given a dictionary with chain objects as values and chain Ids as keys, the
    function checks which chains are similar to other chains in the dictionary.
    Returns a dictionary (similar_chains) with chain Ids as keys and, as values,
    all the chains that are similar to that chain Id according to the threshold argument.
    """
    if main.options.verbose:
        sys.stderr.write("\nComputing similar chains ... \n")
    Nuc_chains={}
    Prot_chains={}
    Nucleotide_alphabet = 'GAUTC'
    for chain_id in chains:
        sequence=get_chain_sequence(chains[chain_id])
        counter=0
        for letter in set(sequence):
            if letter not in Nucleotide_alphabet:
                counter+=1
        if counter == 0:
            Nuc_chains[chain_id]=chains[chain_id]
        if counter > 0:
            Prot_chains[chain_id]=chains[chain_id]

    similar_chains = {}
    for chain_1, chain_2 in itertools.combinations(Nuc_chains, 2):
        sequence_1=get_chain_sequence(chains[chain_1])
        sequence_2=get_chain_sequence(chains[chain_2])
        if sequence_1 and sequence_2:
        	alignments = pairwise2.align.globalxx(sequence_1, sequence_2)
	        pairwise_sequence_identity = alignments[0][2]/len(alignments[0][0])
	        if pairwise_sequence_identity >= sequence_identity_threshold:
		        if chain_1 not in similar_chains:
		            similar_chains[chain_1] = []
		        if chain_2 not in similar_chains:
		            similar_chains[chain_2] = []
		        similar_chains[chain_1].append(chain_2)
		        similar_chains[chain_2].append(chain_1)

    for chain_1, chain_2 in itertools.combinations(Prot_chains, 2):
        alignments = pairwise2.align.globalxx(get_chain_sequence(chains[chain_1]), get_chain_sequence(chains[chain_2]))
        pairwise_sequence_identity = alignments[0][2]/len(alignments[0][0])
        if pairwise_sequence_identity >= sequence_identity_threshold:
	        if chain_1 not in similar_chains:
	            similar_chains[chain_1] = []
	        if chain_2 not in similar_chains:
	            similar_chains[chain_2] = []
	        similar_chains[chain_1].append(chain_2)
	        similar_chains[chain_2].append(chain_1)

    return similar_chains


def clashing_filter(chain_one, chain_two, max_clashes=30, contact_distance=2.0):
    """
    Given a maximum number of clashes and a contact distance, the function checks for clashes
    in the alpha carbons (CA) of two different chains (chain_one and chain_two).

    Returns true if the maximum number of clashes is reached, else returns false.
    """
    atoms_one = [atom for atom in chain_one.get_atoms() if atom.get_id() == 'CA' or atom.get_id() =="C5'"]
    atoms_two = [atom for atom in chain_two.get_atoms() if atom.get_id() == 'CA' or atom.get_id() =="C5'"]
    ns = pdb.NeighborSearch(atoms_one)
    clashes = 0

    for atom_two in atoms_two:
        for atom in ns.search(atom_two.get_coord(), contact_distance, 'A'):
            clashes += 1
            if clashes == max_clashes:
                if main.options.verbose:
                    sys.stderr.write("Clash Found!\n")
                return True

    return False


def fasta_filtering(fasta_filename, chains, pairs, similar_chains, structures, sequence_identity_threshold,options):
    """
    Given a Fasta file this function will keep the chains that, according to a
    sequence identity threshold, are similar with any of the chains from the Fasta
    file. The other chains will be removed.
	"""
    chain_ids_to_keep=[]
    for identifier, fasta_sequence in FASTA_iterator(fasta_filename):
        for chain_id in chains: #for chain_id in similar_chains:
            sequence_from_pdb=get_chain_sequence(chains[chain_id])
            if sequence_from_pdb:
	            alignments = pairwise2.align.globalxx(sequence_from_pdb, fasta_sequence)
	            pairwise_sequence_identity = alignments[0][2]/len(alignments[0][0])
	            if pairwise_sequence_identity >= sequence_identity_threshold:
	                if chain_id not in chain_ids_to_keep:
	                    chain_ids_to_keep.append(chain_id)
    #
    all_chains_ids=[]
    for chain_identifier in chains:
        all_chains_ids.append(chain_identifier)
    for chain_id in all_chains_ids:
        if chain_id not in chain_ids_to_keep:
            del chains[chain_id]
            for listp in pairs:
                if chain_id in listp:
                    listp.remove(chain_id)
            if chain_id not in similar_chains:
                similar_chains[chain_id]=[]
            del similar_chains[chain_id]
            for chain in similar_chains:
                if chain_id in similar_chains[chain]:
                    similar_chains[chain].remove(chain_id)
            for tuplest in structures:
                if chain_id in tuplest:
                    structure=structures[tuplest]
                    list_from_tuple=list(tuplest)
                    list_from_tuple.remove(chain_id)
                    del structures[tuplest]
                    for model in structure:
                        for chain in model:
                            if chain.id == chain_id:
                                del model[chain.id]
                                if options.verbose:
                                    print("\n",chain.id, " chain was deleted.")
                    structures[tuple(list_from_tuple)]=structure


    return (chains, pairs, similar_chains, structures)


def superimpose_chains(chain_structure_one, chain_structure_two):
    """
    Superimpose two chains or structures returning a superimposer object.
    """
    chain_one = copy.deepcopy(chain_structure_one)
    chain_two = copy.deepcopy(chain_structure_two)
    super_imposer = pdb.Superimposer()
    atoms_one = sorted(list(chain_one.get_atoms()))
    atoms_two = sorted(list(chain_two.get_atoms()))

    min_len = min(len(atoms_one), len(atoms_two))
    atoms_one = atoms_one[:min_len]
    atoms_two = atoms_two[:min_len]
    super_imposer.set_atoms(atoms_one, atoms_two)

    return super_imposer


def get_possible_structures(chain_in_current_complex, similar_chains,structures, used_pairs, clashing):
    """
    Taking into account the clashing, previous attemps (used_pairs) and chains
    similarity, this function selects the structures that may be included in
    the current complex.
    """
    possible_structures = {}
    if chain_in_current_complex in similar_chains:
        for similar_chain in similar_chains[chain_in_current_complex]:
            for tuple_key in structures:
                if similar_chain in tuple_key:
                    if tuple_key not in used_pairs:
                        if tuple_key not in clashing:
                            possible_structures[similar_chain] = (tuple_key)

    return possible_structures

# _________________________________
# __ Writing and showing outputs __
# _________________________________
def export_structure(structure, name, format):
    """
    Writes the strucuture into a file. The file can be either a pdb or a mmcif.
    """
    if format == "pdb":
        io = PDBIO()
    elif format == "cif":
        io = pdb.MMCIFIO()
    io.set_structure(structure)
    io.save(name)


def open_in_chimera(options,dirName):
    """
    If actived the --chimera argument by the user, open the models in chimera.
    It doesn't work on Windows.
    """
    for file in os.listdir(os.getcwd()+'/'+dirName):
        if file.endswith('.cif'):
            if options.verbose:
                sys.stderr.write('Opening model %s in Chimera' % file)
            if sys.platform == 'darwin':
                os.system('/Applications/Chimera.app/Contents/MacOS/chimera results/' + file)
            else:
                os.system('chimera results/' + file)
        if file.endswith('.pdb'):
            if options.verbose:
                sys.stderr.write('Opening model %s in Chimera' % file)
            if sys.platform == 'darwin':
                os.system('/Applications/Chimera.app/Contents/MacOS/chimera results/' + file)
            else:
                os.system('chimera results/' + file)
