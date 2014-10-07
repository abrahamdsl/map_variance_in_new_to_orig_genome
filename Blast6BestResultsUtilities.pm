=doc
  a.llave@irri.org 03OCT2014-1605

  A module that contains the BLAST outfmt 6 results processing subroutines for our
  mission.

  Some subroutines need a MySQL database and so the related schema should be accompanying this file.
=cut

require Exporter;
@ISA = qw( Exporter );
#@EXPORT = qw();

sub buildBLAST6LineHash {
=doc
  @devnote Sourced from https://github.com/abrahamdsl/blast_and_infer_subfeatures/blob/master/blast_and_infer_genes_and_subfeatures.pl
  Accepts a line from BLAST result outfmt 6 and builds a hash accoridngly.

  Arguments:
    0 - string. The line/one result or match from BLAST results file.

  Returns:
    Reference to a hash. The hash then is accessed by BLAST6/Table labels
=cut
  my @temparr;
  my %thisHash;
  chomp $_[0];
  @temparr = split( /\s+/, $_[0] );
  $thisHash{ 'query_label' } = $temparr[ 0 ];
  $thisHash{ 'target' } = $temparr[ 1 ];
  $thisHash{ 'percent_identity' } = $temparr[ 2 ];
  $thisHash{ 'alignment_length' } = $temparr[ 3 ];
  $thisHash{ 'mismatch' } = $temparr[ 4 ];
  $thisHash{ 'gap' } = $temparr[ 5 ];
  $thisHash{ 'query_start' } = $temparr[ 6 ];
  $thisHash{ 'query_end' } = $temparr[ 7 ];
  $thisHash{ 'target_start' } = $temparr[ 8 ];
  $thisHash{ 'target_end' } = $temparr[ 9 ];
  $thisHash{ 'evalue' } = $temparr[ 10 ];
  $thisHash{ 'bitscore' } = $temparr[ 11 ];

  return \%thisHash;
} # sub buildBLAST6LineHash

sub insertBLAST6ToMySQL {
=doc
  @devnote Sourced from https://raw.githubusercontent.com/abrahamdsl/blast_and_infer_subfeatures/master/blast_and_infer_genes_and_subfeatures.pl
    and changed as needed.
  Populates our main table in our own-devised MySQL database for BLAST 6 results manipulation.
 
  Arguments:   
    0 - string. The dtabase name.
    1 - string. The table name.
    2 - reference to a HASH. The hash contains the BLAST matches, and whose keys are the line numbers
      from BLAST outfmt 6 results file.

  Returns:
    Nothing
=cut
  my $connection = ConnectToMySql( $_[0] );  
  my $query;
  my %resultsHash = %{ $_[2] };
  my $statement;
  print join ( ' ', @_, "\n" );
 
  $query = "INSERT INTO `$_[1]` (`id`,`query_label`,`target`,`percent_identity`,`alignment_length`,`mismatch`,`gap`,";
  $query .= "`query_start`,`query_end`,`target_start`,`target_end`,`evalue`,`bitscore`)";
  $query .= "  VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?);";  
  $statement = $connection->prepare( "TRUNCATE TABLE `$_[1]`;" );
  $statement->execute();
  foreach my $lineNum ( sort keys %resultsHash ) {
    my %localHash = %{  $resultsHash{ $lineNum } };
    $statement = $connection->prepare( $query );
    $statement->execute(
      $lineNum,
      $localHash{ 'query_label' },
      $localHash{ 'target' },
      $localHash{ 'percent_identity' },
      $localHash{ 'alignment_length' },
      $localHash{ 'mismatch' },
      $localHash{ 'gap' },
      $localHash{ 'query_start' },
      $localHash{ 'query_end' },
      $localHash{ 'target_start' },
      $localHash{ 'target_end' },
      $localHash{ 'evalue' },
      $localHash{ 'bitscore' }
    )  || print "$DBI::errstr \n" ;
  }
} # sub insertBLAST6ToMySQL

sub parseBLAST6Results {
=doc
  Parses BLAST results where output format is 6 (tab-delimited, ala CSV), gets the best matching
  one as per our crteria and returns it to the caller.
  
  Arguments:
    0 - string. BLAST results output file
    1 - string. Feature (i.e. Gene) name
    2 - int. The zero-index, i-th time this function was called for the particular feature. 
          Useful for re-BLASTing short sequences where there were no results in the regular BLAST 
          (with default parameters).
    3 - string. The database name for BLAST6 results filtering.
    4 - string. The table name for the blast6 results.

  Returns:
    Hash accessed by BLAST6/Table labels

    Special note, the 'target_start' and 'target_end' are both assigned the following values as per
      the description:

    -4 - There is no BLAST match at all, for first time this was called
    -5 - There is no BLAST match at all, the second time this was called
=cut
  my @temparr;
  my %thatHash;
  my $lineCount = getFileNumLines( $_[0] );
  my $featureName = $_[1];
  my $timesCalled = $_[2]; 

  print '[debug]' . join( ' ' , @_, "\n" );
# Is result just one line?
  if ( $lineCount == 1 ) {
    return %{
      buildBLAST6LineHash( getSingleLineFromFile( $_[0] ) )       
    };
  }elsif ( $lineCount == 0 ){   
      handle_message('ERROR', "Zero BLAST match for $featureName !", '');
      $thatHash{ -1 } = {        
	'target_start' => ( $timesCalled == 0 ) ? -4 : -5,
	'target_end' => ( $timesCalled == 0 ) ? -4 : -5
      };
      return %{ $thatHash{ -1 } };
  }else{
    if ( $debug > 1 ) { handle_message('NOTICE', 'More than one BLAST match detected for ' . $featureName, ''); }
    my %theResults;
    my $x = 0;
    my @arraytemp;

    # scan matches and put into hash
    open( my $fileHandle, $_[0] ) or die( handle_message('ERROR', "Can't open $_[0] for processing.", '') );
    while ( <$fileHandle> ) {
       $theResults{ $x++ } = buildBLAST6LineHash( $_ );
    }
    close( $fileHandle );
    # insert into database for easier selection
    insertBLAST6ToMySQL( $_[3] , $_[4], \%theResults );
    # We reuse $x !
    # Now eventually filter matches according to our criteria, progressively
    $x =  searchBestBLAST6Match( 0, $featureName, \@arraytemp, '', $_[3], $_[4] ); 
    if ( $x == -1 ) {
      handle_message('ERROR', 'No suitable BLAST match found for ' . $featureName, '' );
    } 
    if( ! exists( $theResults{ $x } ) ) {
       handle_message('ERROR', 'Incorrect BLAST match ' . $x . ' found for ' . $featureName, '' );
    }
    if ( $debug > 1 ) { print "[notice] Successfully found match $x for $featureName\n"; }
    return %{ $theResults{ $x } };
  }
} # parseBLAST6Results

sub searchBestBLAST6Match{
=doc
 @devnote Sourced from https://raw.githubusercontent.com/abrahamdsl/blast_and_infer_subfeatures/master/blast_and_infer_genes_and_subfeatures.pl

 This is a recursive function, progressively filtering out rows. If at the upper level, there are still more than 
 one row, the topmost one will be returned.
 
=cut
  my $recursionLevel = $_[0] ;
  my $geneName = $_[1];  
  my @previousValues = @{ $_[2] };
  my $previousQuerySpecific = $_[3];  
  my $previousTable = $_[5];   
  my $sqlQueryOrderBy; 
  my $sqlQuerySpecific;
  my $dispatchTable;
  my $connection = ConnectToMySql($_[4]);
  my $statement;
  my $thisRecursionFeaturedValue;
  my $thisRecursionFeaturedColumn;
  my $thisTempTable;
  my $returnThis;

  if ( $debug > 1 ) { print "\t Recursion Level $recursionLevel \n"; }
  $dispatchTable = {
    0 => sub {
       $thisTempTable = "lessstrict_0";
       $thisRecursionFeaturedColumn = 'evalue';
#       $previousTable = $blast6Table;
       $sqlQueryOrderBy = "SELECT DISTINCT($thisRecursionFeaturedColumn) FROM $previousTable ORDER BY `$thisRecursionFeaturedColumn`;";
       $sqlQuerySpecific = "SELECT * FROM `specifictable` WHERE `$thisRecursionFeaturedColumn` = ? ";
      },
    1 => sub {
       $thisTempTable = "lessstrict_1";
       $thisRecursionFeaturedColumn = 'alignment_length';
       $sqlQueryOrderBy = "SELECT DISTINCT($thisRecursionFeaturedColumn) FROM $previousTable t ORDER BY t.$thisRecursionFeaturedColumn DESC; ";
       $sqlQuerySpecific = $previousQuerySpecific . " AND $thisRecursionFeaturedColumn = ? "; 
    },
    2 => sub {
       $thisTempTable = "lessstrict_2";
       $thisRecursionFeaturedColumn = 'percent_identity';
       $sqlQueryOrderBy = "SELECT DISTINCT($thisRecursionFeaturedColumn) FROM $previousTable u ORDER BY u.$thisRecursionFeaturedColumn DESC; ";
       $sqlQuerySpecific =  $previousQuerySpecific . " AND $thisRecursionFeaturedColumn = ? ";
    },
    3 => sub {
       $thisTempTable = "lessstrict_3";
       $thisRecursionFeaturedColumn = 'bitscore';
       $sqlQueryOrderBy = "SELECT DISTINCT($thisRecursionFeaturedColumn) FROM $previousTable v ORDER BY v.$thisRecursionFeaturedColumn DESC; ";
       $sqlQuerySpecific =  $previousQuerySpecific . " AND $thisRecursionFeaturedColumn = ? ";
    },
    4 => sub {
       $thisTempTable = "lessstrict_4";
       $thisRecursionFeaturedColumn = 'mismatch';
       $sqlQueryOrderBy =  "SELECT DISTINCT($thisRecursionFeaturedColumn) FROM $previousTable w ORDER BY w.$thisRecursionFeaturedColumn; ";
       $sqlQuerySpecific = $previousQuerySpecific . " AND $thisRecursionFeaturedColumn = ? ";
    },
    5 => sub {
       $thisTempTable = "lessstrict_5";
       $thisRecursionFeaturedColumn = 'gap';
       $sqlQueryOrderBy =  "SELECT DISTINCT($thisRecursionFeaturedColumn) FROM $previousTable x ORDER BY x.$thisRecursionFeaturedColumn; ";
       $sqlQuerySpecific = $previousQuerySpecific . " AND $thisRecursionFeaturedColumn  = ? ";
    },
    6 => sub {
      $thisTempTable = "lessstrict_6";
      print "[fatal] Exhausted all filters for $geneName at recursionLevel $recursionLevel ,still duplicated .. | $previousQuerySpecific \n";
      return "";
    },
   'default' => sub {
      print "[error] Unrecognized $recursionLevel \n";
      return "";
   }
  };

  # Execute dispatch table
  $dispatchTable->{ $recursionLevel }->();
  $statement = $connection->prepare(  $sqlQueryOrderBy  );
  $statement->execute() || print "$DBI::errstr \n";
  
  # Get the topmost row  
  my $rows =  $statement->rows;
  if( $rows < 1 ) {
    print "[fatal] At gene $geneName recursionLevel $recursionLevel No rows at ORDER by part \n";
    $returnThis = -1;
  } 
  # save row in a hash
  my $ref = $statement->fetchrow_hashref() ;
  # get the best value for this recursion
  $thisRecursionFeaturedValue = $ref->{ $thisRecursionFeaturedColumn };

  # some string changes so we can substitute the SQL query, and select only the appropriate data
  # to be progressively filtered
  my $localQuerySpecific =  $sqlQuerySpecific;
  my $quotemetaED = quotemeta( $previousTable );
  $localQuerySpecific =~ s/specifictable/$quotemetaED/;
  my $statement2 = $connection->prepare( "INSERT INTO `$thisTempTable` ( $localQuerySpecific  ) ;"  );
  push( @previousValues, $thisRecursionFeaturedValue );
  $statement2->execute( @previousValues )  || print "$DBI::errstr \n";;
  
  my $statement3 =  $connection->prepare("SELECT * FROM $thisTempTable ;");	
  $statement3->execute();	

  # zero condition not possible (this function won't be called if that's the case), it's either 1 or > 1
  my $s3rows =  $statement3->rows;
  if( $s3rows == 1 )
  {
    my $finalRef = $statement3->fetchrow_hashref();
    $returnThis = $finalRef->{ 'id' };
  }else{
    if ( $recursionLevel < 5 ) {
      if ( $debug > 1) { print "[notice] For $geneName $thisRecursionFeaturedColumn  $thisRecursionFeaturedValue is still not enough. Leveling up...\n"; }
      $returnThis = searchBestBLAST6Match( $recursionLevel + 1,  $geneName, \@previousValues, $sqlQuerySpecific, $_[4], $thisTempTable);
    }else{
       # STill there are more results in the topmost criteria, return the top most
      my $finalRef = $statement3->fetchrow_hashref();
      if ( $debug > 1) { handle_message('WARNING', 'There are still more than one results in recursionLevel 5/max criterion. Returning ' . $finalRef->{ 'id' } . " for $geneName", ''); }
      $returnThis = $finalRef->{ 'id' };
    }
  }
  # TRUNCATE ANY TABLE TEMPORARY
  my $statement5 = $connection->prepare("TRUNCATE TABLE $thisTempTable ;");
  $statement5->execute()  || print "$DBI::errstr \n";;
  return $returnThis;
}#sub

