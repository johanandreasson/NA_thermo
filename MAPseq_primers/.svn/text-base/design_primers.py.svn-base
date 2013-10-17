#!/usr/bin/python

from nupack_wrappers import base_pair_exposure_string, get_base_pair_exposure
from random import randint

illumina1 = 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'
tail2 = 'GTTGTTGTTGTTGTTTCTTT'
barcode_length = 12
num_barcodes = 4
dna_chars = ['A','C','G','T']

def evaluate_base_pair_exposure( sequences ):
    v_tot = 0
    for n in range( len( sequences ) ):
        sequence = sequences[n]
        v =  get_base_pair_exposure( sequence )
        print sequence
        print base_pair_exposure_string( v ), sum(v)
        v_tot += sum( v )*sum(v)  # really penalize large values.
    return v_tot

def do_a_move( sequences ):
    nseq = len( sequences )

    pos = randint( 0, barcode_length-1 ) + len( illumina1 )

    partner1 = randint( 0, nseq-1)
    partner2 = partner1
    while partner2 == partner1: partner2 = randint( 0, nseq-1 )

    print "Switching %d<-->%d at position %d." % (partner1, partner2, pos )
    x = sequences[ partner1 ][ pos ]

    sequences[ partner1 ]  = sequences[ partner1 ][ :pos ] +  sequences[ partner2 ][ pos ] +  sequences[ partner1 ][ (pos+1): ]

    sequences[ partner2 ]  = sequences[ partner2 ][ :pos ] +  x  +  sequences[ partner2 ][ (pos+1): ]



def initialize():
    # initialize barcode combo:
    sequences = []
    for m in range( num_barcodes ):
        sequence  = illumina1
        for n in range( barcode_length ):
            sequence += dna_chars[ m ]
        sequence += tail2
        sequences.append( sequence )

    return sequences

sequences = initialize()

NITER = 1000
for n in range( NITER ):
    do_a_move( sequences )


NITER = 100
v_tot_current = 100000
for n in range( NITER ):

    print "Iteration: %d " % n

    sequences_test = []
    # force deep copy
    for n in range( len( sequences ) ): sequences_test.append( sequences[n] )

    do_a_move( sequences_test )
    v_tot = evaluate_base_pair_exposure( sequences_test )
    print v_tot, " compared to ", v_tot_current

    if (v_tot < v_tot_current):
        sequences = sequences_test
        v_tot_current = v_tot
        print "ACCEPTED"

    print

v_tot = evaluate_base_pair_exposure( sequences )
print v_tot
print "FINAL"
