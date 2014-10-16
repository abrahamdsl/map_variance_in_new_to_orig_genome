#!/usr/bin/perl
=doc
  a.llave@irri.org 01OCT2014-1040
 
  This is a program to map variances that were retrieved on a created, mostly similar genome, to the
  origin/base genome. Takes in variances (SNPs only for now) in VCF format, gets its nearby X sequences
  to the left and right and these are BLAST-ed into the original genome.
  
  Best hit is selcted and the sequence is reduced again to just one base pair (SNPs only for now ) to
  infer the location in new genome.
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

# For debugging
use Data::Dumper;

# For parameters
use Getopt::Long;

# Now our own libraries...
use lib '.';
# Contains general functions for output handling, etc.
use generalUtilities;
# Contains the activity-centric common subroutines. Make sure, it is always synced with 
# https://github.com/abrahamdsl/blast_and_infer_subfeatures/blob/master/InferAnnotationUtilities.pm
use InferAnnotationUtilities;
# For BLAST activities
use Blast6BestResultsUtilities;

# Debugging levels
# 0 - no
# 1 - only those benchmarks
# 2 - reserved
# 3 - all
my $debug = 3;

# For BLAST results selection, manipulation
# Let's forget about the ports for now

# Temporary table to store more than one BLAST 6 results
my $blast6Table = 'blast6out';
my $databaseName = '';
my $databaseUsername = '';
my $databasePassword = '';
my $databaseHost = 'localhost';

# the variables

# benchmarking!
my $benchmark_grandStart = new Benchmark;

# Describes shortly what the features are about. Should be min length of 1.
my $datasetTag = 'null';

# This is the number added to the left and right of the variance (i.e. SNPs) so that we will have 
# a BLAST-able sequence. Default 100
my $endBuffer = 100;

# Indicator if user wants to run the program just for help.
my $help;

# If set to 1, warns if the new target locus does not match via regex, the original locus of the
# submitted annotation/$fascistAnnotation.
# Test case: Chr01 in Chr01_genometitle_programtitlev99
my $detectMismatch = 1;

# New source to specify, for the new data to be generated. Dot as default.
my $newSourceName = '.';

# return destination of GetOptions(), for parameter processing
my $opt_success;

# some tag to be assigned to files outputted by this program
# generator sourced from http://stackoverflow.com/questions/801347/how-should-i-generate-a-random-alphanumeric-initial-password-for-new-users
my $programTag = join '', map +(0..9,'A'..'Z')[rand(10+26*2)], 1..8;

# Used for argument, you can tell the program to skip processing first n lines
# Useful if for some reason, program stops working or crashes due to OS issues,
# so you don't have to start from start again. ( LOLWHUT? )
my $startAt = -1;

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
my $finalOutput_fileHandle;

# Counters, in loops
my $ax;
my $bx;

# Take note variables...
my $probableLocusMismatch = 0;

# start procedures ...
splash();
$opt_success = GetOptions(
  'help' => \$help,
  'blasttb=s' => \$blast6Table,
  'string_tag=s' => \$datasetTag,
  'genome1=s' => \$communism,
  'genome2=s' => \$fascism,
  'gff1=s' => \$communismGFF,
  'gff2=s' => \$fascismAnnotation,
  'endbuffer=s' => \$endBuffer,
  'debug=s' => \$debug,
  'dsn=s' =>\$databaseName,
  'dbhost=s' => \$databaseHost,
  'dbuser=s' => \$databaseUsername,
  'dbpass=s' => \$databasePassword,
  'detectmismatch=i' => \$detectMismatch,
  'newsource=s' => \$newSourceName,
  'start=i' => \$startAt
);
die guide() if ( $help || ! $opt_success );
$benchmark_grandStart = new Benchmark;
handle_message( 'NOTICE', 'Info', "Run ID : $programTag" );
handle_message( 'NOTICE', 'File Open', "Opening genome1 : $communism ...");
handle_message( 'FATAL', 'cannot_open_file', $communism ) unless open( $communismHandle, $communism );
handle_message( 'NOTICE', 'File Open', "Opening genome2 : $fascism ...");
handle_message( 'FATAL', 'cannot_open_file', $fascism ) unless  open( $fascismHandle, $fascism );
handle_message( 'NOTICE', 'File Open', "Opening gff1 : $communismGFF ...");
handle_message( 'FATAL', 'cannot_open_file', $communismGFF ) unless open( $communismGFFHandle, $communismGFF );
handle_message( 'NOTICE', 'File Open', "Opening gff2 : $fascismAnnotation ...");
if( $startAt > -1 ) {
  handle_message('NOTICE','Skipping','Processing will start at line ' . $startAt );
}
handle_message( 'FATAL', 'cannot_open_file', $fascismAnnotation ) unless  open( $fascismAnnotationHandle, $fascismAnnotation );

# Load original genome and its annotation
$communismDBHandle = loadBioDBSeqFeature( $communism, $communismGFF );

# Now, determine if annotation for new genome is GFF3 or internal type
##  cancel that  -  										# Oops, cancel that, let's concentrate on using the internal type.
if ( isGFF3File( $fascismAnnotation ) ) {
   handle_message( 'DEBUG', 'debug', "$fascismAnnotation is GFF3" );
   handle_message( 'FATAL', 'debug', "GFF3 not yet supported. Program exiting." );
}else{
   handle_message( 'DEBUG', 'debug', "$fascismAnnotation is not GFF3. Might be internal" );
   %fascistFeatures =  loadInternalVK_SIMP_1( $fascismAnnotation ) ;
}

# Load new genome, without its annotation, it's not that important anyway
$fascismDBHandle = loadBioDBSeqFeature( $fascism, '' );

# Determine out final output filename
( my $sourceAnnotSansExt = $fascismAnnotation ) =~ s/\.\w+$//; # remove last extension
my @tempx = split( /\//, $sourceAnnotSansExt );
$finalOutputFileName = $tempx[ scalar( @tempx ) - 1 ] . ".variance-mapped_run-$programTag.gff3";

# Now iterate throughout each variance
$ax = 0;
my $averageTimeSoFar = 0;
my $remainingTimeSoFar=0;
my @aliciaKeys = sort keys %fascistFeatures;
$bx = scalar( @aliciaKeys );
my $blastDBIndexCMD = 'makeblastdb -dbtype nucl -in ';
my $blastCommand = "blastn -task blastn -num_threads 8 -outfmt 6 ";
my $searchProperLoopSecs = 0;

# Open our final output handle
handle_message( 'NOTICE', 'File Open', "Opening final output file $finalOutputFileName ..." );
open ( $finalOutput_fileHandle, "> $finalOutputFileName" ) or
  handle_message(
    'FATAL',
    'File open error',
    "Cannot open $finalOutput_fileHandle for outputting of final GFF3 output"
  );
# Some headers
print $finalOutput_fileHandle "##gff-version 3\n";
print $finalOutput_fileHandle "# Generated " . localtime() . " via https://github.com/abrahamdsl/map_variance_in_new_to_orig_genome \n";
print $finalOutput_fileHandle "# Sourced from $fascismAnnotation\n";
print $finalOutput_fileHandle "#\n";

foreach my $hashKey ( sort keys %fascistFeatures ) {
  my %bestResult;

  # The line to be written on the GFF file.
  my $finalGFFLine = "";
  
  my %thisVariance = %{ $fascistFeatures{ $hashKey } };
  my $thisVarianceBLASTfile; 
  my $thisVarianceFeatureID;
  my $thisVarianceFASTA;
  my $varianceType; 
  my $start_time = new Benchmark;
#debug  print Dumper %thisVariance;
  # The x-1 to the left and x+ to the right base pairs/sequence of this variance
  my $blastThisSeq;
  
  # Coordinates of the sequence we want to retrieve and BLAST on the old genome.
  my $seqStart;
  my $seqEnd;

  # Coordinates on the older genome
  my $inferredStart;
  my $inferredEnd;

  print "[notice] " . localtime() . " Progress " . ($ax + 1) . " \/ $bx " . (($ax+1)/$bx) * 100 ;
  if ($debug > 0) { print "\n"; }else{ print "\r" };

  # If $startAt argument was specified, skip $startAt lines
  if( $startAt > -1 ) {
    if ( $ax != $startAt-1 ) {
      if ( $debug > 1 ) {
        handle_message('DEBUG','Skipping','Line ' . ($ax + 1) . ' skipped' );
      }
      $ax++;   #one of the increments here, there's another one down the road, if this
               # is not reached
      next;
    }else{
      # we have reached the specified start line, reset this variable so this if statement
      # won't be executed
      $startAt = -1; 
    }
  }

  $varianceType = ( exists( $thisVariance{ 'pos2' } ) ) ? 'INDEL' : 'SNP';

  $seqStart =  $thisVariance{ 'pos' } - ( $endBuffer - 1 );
  if ( $seqStart < 0 ) { $seqStart = 0 };
  
  # We don't know yet how to check if this exceeds length of target (chromosome)
  $seqEnd = ( $varianceType eq 'INDEL'  ) ?  $thisVariance{ 'pos2' }
    : $thisVariance{ 'pos' } + $endBuffer;

  # Time to retrieve sequence 
  $blastThisSeq = $fascismDBHandle->fetch_sequence(
    -seq_id => $thisVariance{ 'target' },
    -start => $seqStart,
    -end =>  $seqEnd
  );
#debug  print $blastThisSeq . "\t" .  $thisVariance{ 'target' } . "\t" . $seqStart . "\t" . $seqEnd . "\n";
  $ax++;
  
  # generate feature ID
  $thisVarianceFeatureID =  constructFeatureID (
    1, '_', substr( $varianceType, 0, 0 ), $datasetTag, $ax
  );

  # Generate variance fasta file name.
  $thisVarianceFASTA = generateVarianceNameFile(
    $thisVarianceFeatureID,
    'seq',
    '.fa'
  );
 
  # Now write to file
  if (  writeSingleFASTA(
     $thisVarianceFASTA,
     $thisVarianceFeatureID . " " . $thisVariance{ 'target' } . ":" . $seqStart . "-" . $seqEnd,
     $blastThisSeq     
    ) != 0
  ) {
    exit 1;
  }

  # now index
  if ( execCommandEx(
      $blastDBIndexCMD . ' ' . $thisVarianceFASTA . ' > /dev/null',
      "Indexing FASTA for $thisVarianceFASTA",
      "Something went wrong with BLAST indexing"
    ) != 0
  ) {
    handle_message(
     'FATAL', 
     'Indexing error',
     "Cannot index $thisVarianceFASTA, therefore can't proceed" );
  }

  # generate BLAST 6 results file
  $thisVarianceBLASTfile = generateVarianceNameFile(
    $thisVarianceFeatureID,
    'blastresult',
    '.tmp'
  );

  # then blast
  if ( execCommandEx(
      $blastCommand . "-query $thisVarianceFASTA -db $communism > $thisVarianceBLASTfile",
       "BLAST process launched ", " Something went wrong with BLAST!"
    ) != 0
  ) {
     handle_message( 'FATAL', 'BLASTn Error', 'Something went wrong with BLAST' );
     exit(2);
  }

  # Parse results
  %bestResult = parseBLAST6Results(
    $thisVarianceBLASTfile,
    $thisVarianceFeatureID,
    0,
    $databaseName,
    $blast6Table,
    $databaseUsername,
    $databasePassword,
    $databaseHost,
    $programTag
  );  
  
  if ( ! keys %bestResult ) {
    handle_message(
      'ERROR',
      'Blast 404',
      'Return 404 from parseBLAST6Results(), no BLAST match for ' . $thisVarianceFeatureID . '?'
    );
    next;
  }else{
    if( $bestResult{ 'target_start' } < 0 or $bestResult{ 'target_end' } < 0 ) {
      handle_message(
        'ERROR',
        'BLAST ERROR',
        "No BLAST results at all for $thisVarianceFeatureID :'( "
      );
    }else{
      # We have a match yey!
      # Now restore to original
      $inferredStart = $bestResult{ 'target_start' } + ( $endBuffer - 1 );
      $inferredEnd = $bestResult{ 'target_end' } - $endBuffer;

      # Now, generate GFF3 for that
      my $origSeq  =  $communismDBHandle->fetch_sequence(
        -seq_id =>  $bestResult{ 'target' },
        -start => $bestResult{ 'target_start' },
        -end =>   $bestResult{ 'target_end' }
      );
      my $snpx = $blastThisSeq = $communismDBHandle->fetch_sequence(
        -seq_id =>  $bestResult{ 'target' },
        -start => $inferredStart,
        -end => $inferredEnd
      );
      my $snp_len =  length( $snpx );
      if( $snp_len > 1 or $snp_len == 0 ) {
        # We only handle SNPs for now
        handle_message(
          'ERROR',
          'BLAST Error',
          "Top BLAST match for $thisVarianceFeatureID not a SNP!"
        );
      }else{
        if ( $detectMismatch == 1 ) {
	  my $quotemetaedx = quotemeta( $bestResult{ 'target' } );
	  if ( $thisVariance{ 'target' } !~ /$quotemetaedx/ ) {
	    handle_message(
	      'WARNING',
	      'Probable locus mismatch',
	      'Orig ' .  $thisVariance{ 'target' } . ' does not match new target ' . $bestResult{ 'target' } 
	    );
            $probableLocusMismatch++;
	  }
        }
#debug      print '[result] seq=' . ( ( lc($blastThisSeq) eq lc($origSeq) ) ? "eq" : "ne" ) . ' snp=' . ( ( uc($snpx) eq uc( $thisVariance{ 'change' })  )  ?  "no" : "yes" )    .  ' ' . $thisVariance{ 'orig' } . '>' .  $thisVariance{ 'change' } . ' ? ' . $snpx . "\n";

        # Compose the final GFF3 lines
        $finalGFFLine = $bestResult{ 'target' } . "\t$newSourceName\tSNP\t" .  $inferredStart . "\t" . $inferredEnd;
        $finalGFFLine .= "\t" .  $bestResult{ 'evalue' }  ."\t+\t.\t";
        $finalGFFLine .= "ID=$thisVarianceFeatureID;Name=$thisVarianceFeatureID;Note=$snpx > " . $thisVariance{ 'change' }  ;
        $finalGFFLine .= ", AF ";
        $finalGFFLine .= sprintf("%.2f\n",  $thisVariance{ 'af' } );
#        $finalGFFLine .= ",ORIG " . $thisVariance{ 'target' } . "\n";
        if ( $debug >= 3 ) { print "[gff3] $finalGFFLine"; }
        print $finalOutput_fileHandle $finalGFFLine;
        deleteBlastnRelatedFA( $thisVarianceFASTA );
      }
    } 
    # Benchmarking matters
    my $end_time = new Benchmark;
    my $difference = timediff( $end_time, $start_time );
    my $timeNotif =  "[debug] All processes for gene " . $thisVarianceFeatureID . " took " . timestr( $difference ) . "\n"; 
    $timeNotif =~ /\s+(\d+)\s+wallclock\s+secs/;
    $searchProperLoopSecs += $1;
    $averageTimeSoFar = ( $searchProperLoopSecs / $ax );
    $remainingTimeSoFar = ( $averageTimeSoFar * ( $bx - $ax ) );
    if ( $debug > 0 ) { 
      my $dt = DateTime->now();
      my $expectedArrival = $dt->add( seconds => $remainingTimeSoFar );
      print $timeNotif;
      handle_message(
        'DEBUG',
        'Benchmark',
        "$ax features took $searchProperLoopSecs for an average of $averageTimeSoFar. Estimated arrival at " . $expectedArrival->datetime() . "\n"
      );
    }
  } # ???
} # main for loop
if( $probableLocusMismatch > 0 ) {
  handle_message('NOTICE','End notice', "There might be $probableLocusMismatch probable locus mismatches.");
}
handle_message( 'NOTICE', 'Finished', 'Program finished.' );

#
# end of main, start of subroutines
#

sub constructFeatureID {
=doc
  @devnote Sourced from https://github.com/abrahamdsl/vcf_to_gff3/blob/master/vcf_to_gff3.pl as of 03OCT2014-1046
    Please have this moved to a new library, perhaps a "InferVarianceUtilities.pm" 
    or BioDBSeqFeatureUtilities.pm ?
  Constructs the appropriate feature ID for the feature concerned.
  
  Arguments:
    0 - int. Use simple names or not?
    1 - string. Separator, if applicable.
    2 - string. "I" or "S" or "" ( indel, SNP, none )
    3 - string. Data set tag
    4 - int. The ordinal number of this feature
    5 - string. Source of the feature ( VCF 2nd column, 1-index)
  Returns:
    String. The constructed feature ID.
=cut
  if ( $_[0] ) {
    return $_[3] . $_[1] . $_[2] .  $_[4];
  }else{
    return $_[5] . $_[1] . $_[2] . $_[3] . $_[4];
  }
} # sub

sub deleteBlastnRelatedFA {
=doc
  Deletes the 3 other files created by blast when indexing a fasta file, for the feature in question.

  Arguments:
   0 - string. Required. The filename of the feature's FASTA file.
=cut

  if ( execCommandEx(
      "rm $_[0].n*",
       "Deleting associated files of $_[0] ", " Deletion fail."
    ) != 0
  ) {
     handle_message( 'FATAL', 'Cleanup error', 'Something went wrong with deletion.' );
     return -1;
  }
} #deleteBlastn...

sub execCommandEx {
=doc
  @devnote Sourced from https://raw.githubusercontent.com/abrahamdsl/blast_and_infer_subfeatures/master/blast_and_infer_genes_and_subfeatures.pl
  A subroutine to execute outside programs.

  Arguments:
    0 - string. Required. The command (and arguments to execute).
    1 - string. Required. Message to output upon start of execution.
    2 - string. Optional. Message to output upon execution error.
=cut
 my $pidx;
 my $executeCommand = $_[0];
 my $executeStartMsg = $_[1];
 my $executeErrMsg  =  ! defined $_[2]  ?  $_[2] : 'Execution error.\n';
 
 eval {
   $pidx = fork || exec $executeCommand;
   handle_message( 'NOTICE', $executeStartMsg . ' | PID ' . $pidx, '' );
 };
 if ($?){
   handle_message( 'ERROR', $executeErrMsg  . ' | PID ' . $pidx , '' );
   return $?;
 }
 waitpid( $pidx, 0 );
 
 return 0;
} #sub

sub generateVarianceNameFile {
=doc
   Generates the approriate file name for the variance in question.
 
   Arguments:
     0 - string. A plausible feature ID, returned by constructFeatureID(..)
     1 - string. Short description of the file.
     2 - string. Extension. It can be '.fa', '.gff3' for FASTA and GFF3 An-
       notation formats respectively

   Returns:
     String. The generated file name.
=cut
  return "_tmp_$programTag" . "_$_[1]" . '_' . $_[0] . $_[2];
} # sub

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
  --blasttb        Optional. The table within --dsn to store blast 6 results. Default `blast6out`.
  --endbuffer      Optional. Number of extra bases to include from left and end of gene sequence in original. Default $endBuffer.
  --string_tag     Required. A string of minlength 1 to describe the GFF3-able data that will be generated.
  --genome1        Required. The original genome, in FASTA format.
  --genome2        Required. The newly created genome (the one you're generating a new annotation for), in FASTA format.
  --gff1           Required. The GFF3 for --genome1. The subfeatures are required too. As in, complete.
  --gff2           Required. For --genome2.
                     Can be GFF3 or internal format ( insider tag: VK_SIMP_1 )
  --dsn            Required. The MySQL database name for BLAST 6 results filtering.
  --dbuser         Required. The user for --dsn
  --dbpass         Required. The password for --dbuser
  --dbhost         Optional. The host for --dsn, --dbuser and --dbpass. Default 'localhost'
  --detectmismatch Optional. Check via regex if the detected new location's target/locus matches that of the original. 
                     Ex. \"Chr01\" in \"Chr01_xxx_yyyyv10\". Default 1 ( TRUE ). Set to 0 for FALSE.
  --newsource      Optional. Specify if we have to have a new source to be put in the GFF3 column 2 (1-based indexing).
                   If not specified, uses the original source as specified in --gff1
  --start          Optional. The line number of -gff2 to start processing. Useful if program suddenly terminated on  a previous
                     run and you don't wanna start from top again.

  Ex: $0 --endbuffer=100 -string_tag=rand --dsn=taskXXX --start=2402 ...
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
    $genomicDB_handle = ( $_[1] ne '' ) ?  Bio::DB::SeqFeature::Store->new( 
        -adaptor => 'memory',
        -fasta   => $_[0],
        -gff     => $_[1]
      ) 
      :
      Bio::DB::SeqFeature::Store->new(
        -adaptor => 'memory',
        -fasta   => $_[0]
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
    $hushHash{ $x } =  {
      'target' => $temparr[0],
      'pos' => $temparr[1],
      'orig' => $temparr[2],
      'change' => $temparr[3],
      'af' => $temparr[4]
    } ;
    $x++;
  }
  close( $fileHandle );
  #print Dumper %hushHash;
  return %hushHash;
} # sub

sub writeSingleFASTA {
=doc
  Writes to a text file, FASTA sequence of only a single feature.

  Arguments:
    0 - string. File name to write.
    1 - string. Feature name.
    2 - string. The DNA/protein sequence.

  Returns:
    0  - on success
    -1 - any failure   
=cut
  my $fileHandle;
  my $failureDispatch = {
    'default' => sub {
      handle_message( 'ERROR', 'FASTA writing', "Cannot open $_[0]" );
      return -1;
    }
  };

  open ( $fileHandle,"> $_[0]" ) or return $failureDispatch->{ 'default' }->( $_[0] );
  print $fileHandle '>' . $_[1] . "\n";
  print $fileHandle $_[2] . "\n";
  close $fileHandle;

  return 0;
} #sub

sub splash {
  print "\nVariance Mapper to Original Genome Script\n";
  print "a.llave\@irri.org 01OCT2014-1654\n";
  print "---------------------------------------------\n\n";
}

