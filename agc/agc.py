#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
#import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
#import nwalign3 as nw
from Bio import SeqIO

__author__ = "Hollier Laëtitia"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Hollier Laëtitia"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Hollier Laëtitia"
__email__ = "laetitia-hollier@outlook.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()



##########################################################
######## 1. Dé-duplication en séquence “complète" ########
##########################################################


def read_fasta(amplicon_file, minseqlen):
    """la fonction read_fasta prend:
    - amplicon_file: un fichier fasta.gz (str)
    - minseqlen: longueur minimale des séquences (int)
    retourne un générateur de sequences de longueur
    l >= minseqlen : yield sequence
    """
    with gzip.open(amplicon_file, "rt") as filin:
        for record in SeqIO.parse(filin, "fasta"):
            sequence = str(record.seq)
            if len(sequence) >= minseqlen:
                yield sequence


def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    """ https://docs.python.org/2/library/collections.html
    https://docs.python.org/fr/3/library/collections.html
    la fonction dereplication_fulllength prend:
    - amplicon_file: un fichier fasta.gz (str)
    - minseqlen: longueur minimale des séquences (int)
    - mincount: Comptage minimum des séquences (int)
    retourne des séquences (avec O>=mincount) par ordre
    décroissant d'occurence (O)
    """
    sequence = read_fasta(amplicon_file, minseqlen)
    for seq in Counter(sequence).most_common():
        if seq[1] >= mincount:
            yield seq


################################################################
##2. Recherche de séquences chimériques par approche “de novo”##
################################################################


def get_chunks(sequence, chunk_size):
    """La fonction get_chunks prend:
    - sequence: une sequence sous forme d'une chaine de caractères (str)
    - chunk_size: la longueur l de segment (int)
    renvoit une liste de sous-séquence de taille l non chevauchant
    NB: 4 segments doivent être obtenus par séquence !
    """
    liste_segments = []
    i = 0
    while i  < (len(sequence)-chunk_size+1):
        segm = sequence[i:i+chunk_size]
        liste_segments.append(segm)
        i += chunk_size

    return liste_segments


def cut_kmer(sequence, kmer_size):
    """La fonction cut_mer prend:
    - sequence : une sequence sous forme d'une chaine de caractères (str)
    - kmer_size : la taille du k-mer (int)
    renvoit les k-mers uniques présents dans la séquence
    """
    for i in range(len(sequence)-kmer_size+1):
        yield sequence[i:i+kmer_size]


def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    """La fonction prend get_unique_kmer:
    - kmer_dict: ayant pour clé un index de kmer et pour valeur
    une liste d’identifiant des séquences dont ils proviennent
    - sequence: 1sequence sous forme d'une chaine de caractères (str)
    - id_seq: identifiant de la séquence (int)
    - kmer_size: la taille du k-mer (int)
    https://moonbooks.org/Articles/Cl%C3%A9-dun-dictionnaire-avec-plusieurs-valeurs-associ%C3%A9es-sous-python/
    """
    kmers = cut_kmer(sequence, kmer_size)
    for kmer in kmers:
        if kmer not in kmer_dict:
            kmer_dict[kmer] = id_seq
        else:
            kmer_dict[kmer].append(id_seq)
    return kmer_dict



def search_mates(kmer_dict, sequence, kmer_size):
    """La fonction search_mates prend:
    - kmer_dict: un dictionnaire avec pour un clé un index et
    pour valeur une liste d'identifiant des séquences
    - sequence: une sequence sous forme d'une chaine de caractères (str)
    - kmer_size: la taille du kmer (int)
    """
    return [i[0] for i in Counter([ids for kmer
        in cut_kmer(sequence, kmer_size) if kmer in kmer_dict
        for ids in kmer_dict[kmer]]).most_common(8)]


def get_identity(alignment_list):
    """La fonction get_identity prend:
    - alignment_list: liste de 2 séquences (liste de sous
    forme d'une chaine de caractères (str))
    => calcule le pourcentage d'identité entre deux séquences
    selon id = nb nucleotides identiques / longueur de l'alignement
    """
    seq1 = alignment_list[0]
    seq2 = alignment_list[1]
    count = 0
    for base1 in seq1:
        for base2 in seq2:
            if base1 == base2:
                count += 1

    identity = count / len(seq1)
    return identity


##########################################################
################ 3. Regroupement glouton #################
##########################################################


def fill(text, width=80):
    """
    Split text with a line return to respect fasta format
    """
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))



def write_OTU(OTU_list, output_file):
    """ La fonction write_OTU prend:
    une liste d'OTU
    et un chemin vers un fichier de sortie
    affiche les OTU au format fasta
    """
    with open(output_file, "w") as filout:
        for i,seq in enumerate(OTU_list):
            filout.write(">OTU_{} occurence:{}\n".format(i+1,seq[1]))
            filout.write(fill(seq[0]))
            filout.write("\n")
        return filout


#=========================================================
# Main program
#=========================================================

def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()


    #lecture du fichier fasta
    amplicon_file = sys.argv[2]
    minseqlen = int(sys.argv[4])
    mincount = int(sys.argv[6])
    lecture = dereplication_fulllength(amplicon_file, minseqlen,mincount)
    print(lecture)


if __name__ == '__main__':
    main()
