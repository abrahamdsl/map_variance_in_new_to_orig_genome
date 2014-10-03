#!/usr/bin/perl
=doc
  a.llave@irri.org 01OCT2014-1040
 
  This is a program to map variances that were retrieved on a mostly similar genome, to the
  origin/base genome. Takes in variances (SNPs only for now) in VCF format, gets its nearby X sequences
  to the left and right and these are BLAST-ed into the original genome.
  
  Best hit is selcted and the sequence is reduced again to just one base pair to infer the location in
  new genome.
=cut




