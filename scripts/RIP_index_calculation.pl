#!/usr/bin/perl -w
use strict;

# $Id$

=head1 NAME

RIP index calculation

=head1 USAGE

RIP_index_calculation genome.fa > genome.RIP_stats.dat

=head1 DESCRIPTION

This script calculates RIP index for windows across a genome
 
RIP index types

   o "RIP Product Index" (TpA / ApT) (Margolin et al 1998, Selker et
     al 2003) 
     Sequences that have not been subjected to RIP exhibit
     values less than 0.8 while values of 1.1 indiate RIP'd sequence.

   o "RIP Substrate Index" (CpA + TpG / ApC + GpT) (Margolin et al
     1998, Selker et al 2003) 
     Sequences that have not been subject to RIP have values greater 
     than 1.1 while values less than 0.9 indicate RIP'd sequence.

   o "Composite Index" (CRI) = [substrate - product] (Lewis et al
      2009) "A positive CRI value implies that the DNA has been
      subjected to RIP" (Lewis et al).


=head1 OPTIONS

 -t or --type    RIP index type   ('product','substrate','composite' or CRI)
 -w or --window  Window Size       integer value, default is 500 bp
 -s or --step    Window step size  integer value < windowsize, default 100

=head2 References

Lewis ZA, Honda S, Khlafallah TK, Jeffress JK, Freitag M, Mohn F,
Schübeler D, Selker EU. Relics of repeat-induced point mutation direct
heterochromatin formation in Neurospora crassa. Genome Res 2009
19(3):427-37.

Margolin BS, Garrett-Engele PW, Stevens JN, Fritz DY, Garrett-Engele
C, Metzenberg RL, Selker EU. (1998) A methylated Neurospora 5S rRNA
pseudogene contains a transposable element inactivated by
repeat-induced point mutation. Genetics 149:1787–1797.

Selker EU, Tountas NA, Cross SH, Margolin BS, Murphy JG, Bird AP,
Freitag M. (2003) The methylated component of the Neurospora crassa
genome. Nature 422:893–897.

=head1 AUTHORS

Jason Stajich, jason.stajich [AT] ucr.edu


=cut

# BioPerl Modules
use Bio::SeqIO;
 
# System Modules
use Getopt::Long;
#use Statistics::Descriptive;

# variables to be set by cmdline options
my ($index_type);  # one of 'product', 

# get the cmdline options

