#!/usr/bin/perl
=doc
  a.llave@irri.org 01OCT2014-1040
 
  This is a program to map variances that were retrieved on a created, mostly similar genome, to the
  origin/base genome. Takes in variances (SNPs only for now) in VCF format, gets its nearby X sequences
  to the left and right and these are BLAST-ed into the original genome.
  
  Best hit is selcted and the sequence is reduced again to just one base pair to infer the location in
  new genome.
=cut

use warnings;
use strict;

# For measuring our performance
use Benchmark;
use DateTime;

# For database accesses
use Bio::DB::SeqFeature;
use Bio::DB::SeqFeature::Store;
use DBD::mysql;
use DBI;

# For parameters
use Getopt::Long;

# Now our own libraries...
use lib '.';
# Contains general functions for output handling, etc.
use generalUtilities;
# Contains the activity-centric common subroutines. Make sure, it is always synced with 
#   https://github.com/abrahamdsl/blast_and_infer_subfeatures/blob/master/InferAnnotationUtilities.pm
use InferAnnotationUtilities;

# Debugging levels
# 0 - no
# 1 - only those benchmarks
# 2 - reserved
# 3 - all
my $debug = 3;

# For BLAST results selection, manipulation
# Let's forget about the ports for now
my $blast6Table = 'blast6out';
my $databaseName = '';
my $databaseUsername = '';
my $databasePassword = '';
my $databaseHost = 'localhost';

# the variables

# benchmarking!
my $benchmark_grandStart = new Benchmark;


# This is the number added to the left and right of the variance (i.e. SNPs) so that we will have 
# a BLAST-able sequence. Default 100
my $endBuffer = 100;

# Indicator if user wants to run the program just for help.
my $help;

# New source to specify, for the new data to be generated. Dot as default.
my $newSourceName = '.';

# return destination of GetOptions(), for parameter processing
my $opt_success;

# The original genome and its related variables.
my $communism;
my $communismHandle;
my $communismGFF;
my $communismGFFHandle;
my $communismDBHandle;

# The new genome and its related variables
my $fascism;
my $fascismHandle;
my $fascismAnnotation;
my $fascismAnnotationHandle;
my $fascismDBHandle;
my %fascistFeatures;

# The ultimate file we want
my $finalOutputFileName;


# start procedures ...
splash();
$opt_success = GetOptions(
  'help' => \$help,
  'genome1=s' => \$communism,
  'genome2=s' => \$fascism,
  'gff1=s' => \$communismGFF,
  'gff2=s' => \$fascismAnnotation,
  'endbuffer=s' => \$endBuffer,
  'debug=s' => \$debug,
  'newsource=s' => \$newSourceName
);
die guide() if ( $help || ! $opt_success );
$benchmark_grandStart = new Benchmark;
handle_message( 'NOTICE', 'File Open', "Opening genome1 : $communism ...");
handle_message( 'FATAL', 'cannot_open_file', $communism ) unless open( $communismHandle, $communism );
handle_message( 'NOTICE', 'File Open', "Opening genome2 : $fascism ...");
handle_message( 'FATAL', 'cannot_open_file', $fascism ) unless  open( $fascismHandle, $fascism );
handle_message( 'NOTICE', 'File Open', "Opening gff1 : $communismGFF ...");
handle_message( 'FATAL', 'cannot_open_file', $communismGFF ) unless open( $communismGFFHandle, $communismGFF );
handle_message( 'NOTICE', 'File Open', "Opening gff2 : $fascismAnnotation ...");
handle_message( 'FATAL', 'cannot_open_file', $fascismAnnotation ) unless  open( $fascismAnnotationHandle, $fascismAnnotation );

# Load original genome and its annotation
$communismDBHandle = loadBioDBSeqFeature( $communism, $communismGFF );

# Now, determine if annotation for new genome is GFF3 or internal type
##  cancel that  -  										# Oops, cancel that, let's concentrate on using the internal type.
if ( isGFF3File( $fascismAnnotation ) ) {
   handle_message( 'DEBUG', 'debug', "$fascismAnnotation is GFF3" );
}else{
   handle_message( 'DEBUG', 'debug', "$fascismAnnotation is not GFF3. Might be internal" );
   %fascistFeatures =  %{ loadInternalVK_SIMP_1( $fascismAnnotation ) };
   hanlde_message( 'FATAL', 'debug', "Sorry, $fascismAnnotation needs to be GFF3. Please run the accompanying conversion tool." );
}

##test code START
my @communismGenes = $communismDBHandle->get_features_by_type('gene');
my ($temp_arr_ref2, $queryHashRef2) = convert_BioDBSeqFeature_To_Hash( \@communismGenes, 1 );
my %queryFeaturesFast2 = %$queryHashRef2;
my @queryFeatureKeys2 = @$temp_arr_ref2;
foreach my $hashKey (@queryFeatureKeys2) {
  my $obj = $queryFeaturesFast2{ $hashKey };
  print $obj->name . "\n";
}
##test code END

#
# end of main, start of subroutines
#
sub getSingleNonSharpLineFromFile {
=doc
  Gets the chomp()-ped, first non-comment line of a text file, where the comment indicator is a #,
   as used in Perl and Linux Bash scripts.

  Arguments:
    0 - file name

  Returns
    string. The line.
=cut
  my $fileHandle;
  my $strtemp;

  open( $fileHandle, $_[0] )
    or die( handle_message('FATAL', 'File open error', "Can't open $_[0] for processing.", ) );
  while( <$fileHandle> ) {
    chomp $_;
    if ( $_ =~ /^(\s+)?#/ ) { next; };
    $strtemp = $_;
    last;
  }
  close( $fileHandle );

  return $strtemp;
} #sub


sub guide {
  print "

  Synopsis:

  $0 [options]

  Options:
  --endbuffer Optional. Number of extra bases to include from left and end of gene sequence in original. Default $endBuffer.
  --genome1   Required. The original genome, in FASTA format.
  --genome2   Required. The newly created genome (the one you're generating a new annotation for), in FASTA format.
  --gff1      Required. The GFF3 for --genome1. The subfeatures are required too. As in, complete.
  --gff2      Required. For --genome2.
                 Can be GFF3 or internal format ( insider tag: VK_SIMP_1 )
  --newsource Optional. Specify if we have to have a new source to be put in the GFF3 column 2 (1-based indexing).
                If not specified, uses the original source as specified in --gff1
  ";
} # sub

sub isGFF3File {
=doc
  Checks if the file in question is a GFF3-compliant one. Actual implementation, just
    checks the first non-comment line if it is such.

  Arguments:
    0 - string. File name of the file.
 
  Returns.
    Boolean.
=cut
  return ( getSingleNonSharpLineFromFile( $_[0] ) =~ /^\S+\s+\S+\s+\w+\s+\d+\s+\d+\s+\S+\s+\S+\s+[\d|\.]\s+.*ID=.*$/ );
} # sub

sub loadBioDBSeqFeature {
=doc
  @devnote Move to InferAnnotationUtilities.pm when finalized.

  Loads FASTA and GFF to a in-memory database, in the form of Bio::DB::SeqFeature

  Arguments:
    0 - string. File name of genome in FASTA format.
    1 - string. The file name of the annotation in GFF3 format.

  Returns:
    Bio::DB::SeqFeature::Store object
=cut
  my $genomicDB_handle;

  eval {
    $genomicDB_handle = Bio::DB::SeqFeature::Store->new( 
      -adaptor => 'memory',
      -fasta   => $_[0],
      -gff     => $_[1]
    );    
  };

  if ( $@ ) {
    print $@;
    handle_message( 'FATAL', 'File open error', "Error opening $_[0] or $_[1]. Program exiting." );
  }
  return $genomicDB_handle;
} # sub


sub loadInternalVK_SIMP_1 {
=doc
  Reads from text file in the format, which we internally call 'VK_SIMP_1', ex:

    #target	POS	ORIG	CHANGE	AF   
    Chr01_pilonV9   2347043 C       T       0.527777777777778

  ... puts it into a hash and returns hash address.

  Arguments:
    0 - string. The file name of the file.
   
  Returns:
    Reference to a hash. Keys are ordinal number of the feature scanned.
=cut
  my $fileHandle;
  my %hushHash;
  my $x = 0;

  open ( $fileHandle, $_[0] ) or 
    handle_message( 'FATAL', 'File open error', "Error opening $_[0]. Program exiting." );
  while( <$fileHandle> ) {
    my @temparr;

    # There should be at least one pattern of <word><whitespace><word>
    if( $_ !~ /\w+\s+\w+/ ) { next; }

    chomp $_; 
    @temparr = split( /\s+/, $_ );
    $hushHash{ $x } = {
      'target' => $temparr[0],
      'pos' => $temparr[1],
      'orig' => $temparr[2],
      'change' => $temparr[3],
      'af' => $temparr[4]
    };
    $x++;
  }
  close( $fileHandle );

  return \%hushHash;
} # sub

sub splash {
  print "\nVariance Mapper to Original Genome Script\n";
  print "a.llave\@irri.org 01OCT2014-1654\n";
  print "---------------------------------------------\n\n";
}
