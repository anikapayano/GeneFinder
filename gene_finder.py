# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Anika Payano

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

            nucleotide: a nucleotide (A, C, G, or T) represented as a string
            returns: the complementary nucleotide
        >>> get_complement('A')
        'T'
        >>> get_complement('C')
        'G'
    """
    if nucleotide == 'A':
        return 'T'
    if nucleotide == 'C':
        return 'G'
    if nucleotide == 'G':
        return 'C'
    if nucleotide == 'T':
        return 'A'

get_complement('C')




def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
            sequence

            dna: a DNA sequence represented as a string
            returns: the reverse complementary DNA sequence represented as a string
        >>> get_reverse_complement("ATGCCCGCTTT")
        'AAAGCGGGCAT'
        >>> get_reverse_complement("CCGCGTTCA")
        'TGAACGCGG'
    """
    complement= ''
    for nucleotide in dna:
        newtide = get_complement(nucleotide)
        complement = complement + newtide
    return complement[::-1]


get_reverse_complement('ATGCCCGCTTT')




def rest_of_ORF(dna):
    """" Takes a DNA sequence that is assumed to begin with a start
            codon and returns the sequence up to but not including the
            first in frame stop codon.  If there is no in frame stop codon,
            returns the whole string.

            dna: a DNA sequence
            returns: the open reading frame represented as a string
        >>> rest_of_ORF("ATGTGAA")
        'ATG'
        >>> rest_of_ORF("ATGAGATAGG")
        'ATGAGA'
    """
    index = len(dna)
    newdna = dna[3:]
    substrings = ['TAG', 'TAA', 'TGA']
    for i in range(3, len(dna)-2, 3):
        codon = dna[i: i+3]
        if codon in substrings:
            index = i
            return dna[:index]
    return dna[:index]

rest_of_ORF('ATGAAATAAAAGAAAGAAAAAAAA')



def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    orf = []
    i = 0
    start = 'ATG'
    while i < len(dna) - 2:
        codon = dna[i: i+3]
        if codon == start:
            frame = rest_of_ORF(dna[i:])
            i = i + len(frame)
            orf.append(frame)
        else:
            i = i + 3
    return orf

find_all_ORFs_oneframe('ATGCATGAATGTAGATAGATGTGCCC')



def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
            all 3 possible frames and returns them as a list.  By non-nested we
            mean that if an ORF occurs entirely within another ORF and they are
            both in the same frame, it should not be included in the returned list
            of ORFs.

            dna: a DNA sequence
            returns: a list of non-nested ORFs

        >>> find_all_ORFs("ATGCATGAATGTAG")
        ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    one =find_all_ORFs_oneframe(dna)
    two = find_all_ORFs_oneframe(dna[1:])
    three = find_all_ORFs_oneframe(dna[2:])
    return (one + two + three)

find_all_ORFs('ATGCATGAATGTAG')


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    one = find_all_ORFs(dna)
    two = find_all_ORFs(get_reverse_complement(dna))
    return(one + two)

find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    maxlength = max(len(orf) for orf in find_all_ORFs_both_strands(dna))
    for orf in find_all_ORFs_both_strands(dna):
        if len(orf) == maxlength:
            return orf

longest_ORF("ATGCGAATGTAGCATCAAA")


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF
    """
    longestorf = []
    for i in range(num_trials):
        longestorf.append(longest_ORF(dna))
        shuffle_string(dna)
    maxlengthorf = max(len(orf) for orf in longestorf)
    for orf in longestorf:
        if len(orf) == maxlengthorf:
            return orf

longest_ORF_noncoding("ATGCGAATGTAGCATCAAA", 5)

def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    AA = []
    i = 0
    for i in range(0, len(dna)-2, 3):
        codon = dna[i: i+3]
        AA.append(aa_table[codon])
        i = i+3
    return ''.join(AA)

coding_strand_to_AA("ATGCGA")


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    sequences = []
    threshold = longest_ORF_noncoding(dna, 1500)
    for i in find_all_ORFs_both_strands(dna):
        if len(i) >= len(threshold):
            sequences.append(coding_strand_to_AA(i))
    return sequences


gene_finder("ATGCGAATGTAGCATCAAA")

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    from load import load_seq
    dna = load_seq("./data/X73525.fa")

    print(gene_finder(dna))
