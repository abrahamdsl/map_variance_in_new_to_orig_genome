
=doc
 a.llave@irri.org 10APR2014-0924

 This is a module that is required by Perl Scripts present in this directory:

=cut

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(convert_BioDBSeqFeature_To_Hash);

sub appendLineToFile {
=doc
  Appends a single line to a text file.

  Arguments:
    0 - string. The text to be appended.
    1 - string. The file name.
  
  Returns:
    Nothing
=cut
  my $fileHandle;

  if( ! open( $fileHandle, ">> $_[1]" ) ) { 
    print '[fatal] Cannot open ' . $_[0] . " for writing.\n";
    return -1;
  }
  print $fileHandle $_[0] . "\n";
  close( $fileHandle ) ? return 0 : return -1;
} # sub

sub convert_BioDBSeqFeature_To_Hash
{
  #
  # Turns an array of  Bio::DB::SeqFeature objects to hash.
  #
  # Arguments:
  #  0 - array reference
  #  1 - int. 0 for by primary_id, 1 for load_id
  # Returns (ARRAY ref, HASH ref). 
  #  Array ref - elements are the sorted load_ids.
  #  For hash ref, Key: load_id of the  Bio::DB::SeqFeature object.
  #
  my @features = @{ $_[0] };
  my @sortedKeys = ();
  my %tempHash = ();
  my $sortBy = $_[1];
  my $totalCount = scalar( @features );
  my $counter = 1;
  my $updateEvery = 0;

=consider-in-the-future
  if ( scalar(@features) < 1 ){
    print "[error] " . localtime() . " sub convert_BioDBSeqFeature_To_Hash: Features array has no content!\n\n";
    exit (1);
  }
=cut

  foreach my $feature( @features )
  {
    if ( $sortBy == 0 ) {
      $tempHash{ $feature->primary_id } = $feature;
    }elsif ( $sortBy == 1) {
      $tempHash{ $feature->load_id } = $feature;
    }
  }
  if ( $sortBy == 0 ) {
    foreach my $hashKey ( sort { $a <=> $b} keys %tempHash ){
      push (@sortedKeys, $hashKey );
    }
  }elsif( $sortBy == 1){
    foreach my $hashKey ( sort { $a cmp $b} keys %tempHash ){
      push (@sortedKeys, $hashKey );
    }
  }
  return ( \@sortedKeys, \%tempHash );
}

sub convertGFF3Column9_To_Hash
{
  # @devnote Copied from "adjust_feature_coordinates.pl"
  #
  # Converts the attribute part of a Bio::DB::SeqFeature object into a hash
  #
  # Arguments:
  #  0 - The Bio::DB::SeqFeature object
  #  1 - Int. 0 or 1. 0 - Do not include 'load_id', else do include it.
  #
  # Returns:
  #  A duple of the references of the array containing the keys of the
  #  hash in the order as they appear in the GFF3 file, then the reference
  #  to the hash.
  #
  my $featureObj = shift;
  my $keepLoadID = shift;
  my @attributes;
  my @attributesKeyInGFF3Order = ();
  my $x = 0;
  my $lastVal = '';
  my %thisHash = ();

  if( !defined( $featureObj ) ) {
    print "[error] " . localtime() . " Feature object not defined.\n";
    return {};
  }
  @attributes = $featureObj->attributes;
  if( scalar(@attributes) < 1 ){
    return ( (),() );
  }
  if($debug > 1){
    if( $featureObj->type =~ /gene/ ) {
      print Dumper $featureObj;
    }
  }
  foreach my $element( @attributes )
  {
    if($debug) { print "\t\t\t\t\t$element \n"; }
    if( $x % 2 == 0 ){
      $lastVal = $element;
    }else{
      if( $lastVal eq 'load_id' and !$keepLoadID ){
        last;
      }
      $thisHash{ $lastVal } = extractSingleArrayElement( $element );
      push( @attributesKeyInGFF3Order, $lastVal );
      if($debug) { print "\t\t\t\t\t$lastVal : " . $thisHash{ $lastVal } . "\n"; }
    }
    $x++;
  }

  return (\@attributesKeyInGFF3Order, \%thisHash);
} #sub

sub extractSingleArrayElement
{
  # @devnote  Copied from "adjust_feature_coordinates.pl"
  #
  # Returns the single element within an array.
  #
  # Such arrays are common in a Bio::DB::SeqFeature object.
  #
  # Arguments:
  # 0 - array
  #
  my @thisArray = @{ $_[0] };
  my $singleValue;

  # this could spell a problem. `defined()` is deprecated however...
  if ( !( exists( $thisArray[0] ) ) ){
    print "[error] " . localtime() . " extractSingleArrayElement: Array is not defined.\n";
    return '';
  }
  foreach my $x (@thisArray)
  {
    $singleValue = $x;
  }
  return $singleValue;
}

sub getLeftInColonPairing {
=doc
 @devnote Copied from "adjust_feature_coordinates.pl"::simpifyFeatureType() and evolved to here

 Gets the left part in a <left-data>:<right-data> pairing 
 ( 
   ex: <feature-type>:<source> of a Bio::DB::SeqFeature object's 'type' field )

   There are times that the Bio::DB::SeqFeature loader makes the object type in the format <feature-type>:
   <source> instead of just <feature-type>. Example: `gene:phytozome_8.0` instead of just `gene`.So this
   was built to accomodate that.
 )

 Arguments:
   0 - string. The  <left-data>:<right-data> pairing
=cut
 my @tempx = split( /:/, $_[0] );
 return $tempx[0];
}

sub getLiteralStrand {
=doc
  Returns the textual representation of the strand of a feature in Bio::DB::SeqFeature object format.
 
  Arguments:
    0 - Bio::DB::SeqFeature object. The feature in question.
=cut
  if ( $_[0]->strand == 1 ) {
    return '+';
  }elsif ( $_[0]->strand == -1 ) {
    return '-';
  }else{
    print "[warning] Invalid strand " .  $_[0]->strand . "\n";
    return '.';
  }
} #sub

sub getSingleLineFromFile {
=doc
  Gets the chomp()-ped first line of a text file.

  Arguments:
    0 - file name

  Returns
    string. The line.
=cut
  my $fileHandle;
  my $strtemp;
 
  open( my $fileHandle, $_[0] ) or die( handle_message('ERROR', "Can't open $_[0] for processing.", '') );
  $strtemp = <$fileHandle>;
  chomp( $strtemp );
  close( $fileHandle );

  return $strtemp;
}


