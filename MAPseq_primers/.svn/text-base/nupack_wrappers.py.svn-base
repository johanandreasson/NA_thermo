#!/usr/bin/python

from string import *
from os import system
from random import randint
from os.path import exists

EXE_DIR = '/Users/rhiju/src/nupack3.0/bin/'
if not( exists(EXE_DIR) ): EXE_DIR = '/home/rhiju/src/nupack3.0/bin/'
assert( exists( EXE_DIR ) )

def get_base_pair_exposure( sequence ):

    tag = "%d" % randint( 0, 6553612312 )

    system( "echo %s > %s.in" % (sequence, tag) )

    system( EXE_DIR+"pairs -T 25 -material dna %s" % tag )

    seqlength = 0

    v = []

    for line in open( "%s.ppairs" % tag ):
        if len( line )< 1: continue
        if line[0] == '%': continue
        cols = line.split()

        if len( cols ) == 1:
            seqlength = int( cols[0] )
            v = []
            for i in range( seqlength ): v.append( 1.0 )

        if len( cols ) == 3 and int(cols[1]) ==  seqlength + 1:
            v[ int(cols[0])-1 ] = 1.0 - float( cols[2] )

    system( "rm -rf %s.*" % tag )

    return v

def base_pair_exposure_string( v ):
    s = ''
    for n in range( len( v ) ):
        d = int( v[n] * 10 )
        if (d > 9 ): d = 9
        s += '%d' % d
    if len( s) != len( v ): exit()
    return s
