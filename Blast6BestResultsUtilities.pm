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
    0 - string. The database name.
    1 - string. The table name.
    2 - string. The user for arg 0.
    3 - string. The password for arg 2.
    4 - string. The hostname for args 0-3.
    5 - string. Reference to a HASH. The hash contains the BLAST matches, and whose keys are the line numbers
      from BLAST outfmt 6 results file.
    6 - string. Program/Run Tag for this process.

  Returns:
    Nothing
=cut
  my $connection = ConnectToMySql( $_[0], $_[2], $_[3], $_[4] );  
  my $query;
  my %resultsHash = %{ $_[5] };
  my $statement;
 
  $query = "INSERT INTO `$_[1]` (`run_id`,`id`,`query_label`,`target`,`percent_identity`,`alignment_length`,`mismatch`,`gap`,";
  $query .= "`query_start`,`query_end`,`target_start`,`target_end`,`evalue`,`bitscore`)";
  $query .= "  VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?);";  
  $statement = $connection->prepare( "TRUNCATE TABLE `$_[1]`;" );
  $statement->execute();
  foreach my $lineNum ( sort keys %resultsHash ) {
    my %localHash = %{  $resultsHash{ $lineNum } };
    if( exists( $localHash{ 'chr_mismatch_indicator' }  ) ) {
      # The BLAST result for this line in the results file is crossed out as it does not match
      # the chromosome we are targeting.
#      handle_message( 'DEBUG', 'chr_mismatch_indicator present', "$0:insertBLAST6ToMySQL  chromosome mismatch check present, skipping line" ); 
      next;
    }else{
      $statement = $connection->prepare( $query );
      $statement->execute(
	$_[6],
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
    5 - string. The username for arg 3.
    6 - string. The password for arg 3.
    7 - string. The hostname for arg 3.
    8 - string. Program Run ID
    9 - string. Chromosome to loosely match via regex, with the BLAST results, to avoid mismatch of chromosomes.
          If not needed, and by default ''.
  Returns:
    Hash accessed by BLAST6/Table labels
      If "anti-chromosome mismatch" is enabled by supplying orig chromosome of feature in question in parameter 9 (0-index)
        then, the key of the hash contains only a hash whose only key is 'chr_mismatch_indicator'

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
  my $origChromosome = $_[9];

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
       # Here we check for mismatch chromosome in BLAST results, if caller specified.
       # So far I think, this needs editing, the chromosome name length in the original genome being mapped to should
       # be lesser than the new one.
       # Ex: Original: Chr02
       #     New: Chr02_pilonv9
       if( $origChromosome ne '' ) {
         chomp $_;
         # Find the second word from beginning of line. In BLATS output format 6, that's where the chromosome is.
         $_ =~ /^\S+\s+(\S+)/;
         my $new_genome_chr = $1;
         my $quotemetaed_blastresult_target = quotemeta( $new_genome_chr ); # quotemeta the matched from earlier line ( first word from start of line -  a chromosome perhaps )
#         handle_message( 'DEBUG', 'Mismatch chromosome check', $featureName . ' quotemetaed ' . $new_genome_chr . ' orig ' . $origChromosome );
         if ( $origChromosome !~ $quotemetaed_blastresult_target ) {
           if ( $debug > 1 ) {
             handle_message( 'NOTICE', 'Mismatch chromosome excluded for ' . $featureName, 'BLAST Result #' . ( $x + 1 ) . ' excluded due to chromosome mismatch' );
           }
           $theResults{ $x }{ 'chr_mismatch_indicator' } = -10  ;  # mismatch chromosome code 
           $x++;
           next;
           
         }
       }
       $theResults{ $x++ } = buildBLAST6LineHash( $_ );
    }
    close( $fileHandle );
    # insert into database for easier selection
    insertBLAST6ToMySQL( 
      $_[3],
      $_[4],
      $_[5],
      $_[6],
      $_[7],
      \%theResults,
      $_[8]
    ); 
    # We reuse $x !
    # Now eventually filter matches according to our criteria, progressively
    $x = searchBestBLAST6Match(
      0,
      $featureName,
      \@arraytemp,
      '',
      $_[3],
      $_[4],
      $_[5],
      $_[6],
      $_[7],
      $_[8]
    ); 
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

 Arguments:
   0 - int. Recursion level.
   1 - string. Feature name.
   2 - Reference to an Array. Array contains the values for SQL selection/fitlering.
   3 - string. The specific previous query.
   4 - string. database name.
   5 - string. Database table containing the previously filtered results.
   6 - string. User name for arg 4.
   7 - string. Password for arg 6.
   8 - string. Host for args 4,6-8.
   9 - string. Program run tag. 


=cut
  my $recursionLevel = $_[0] ;
  my $geneName = $_[1];  
  my @previousValues = @{ $_[2] };
  my $previousQuerySpecific = $_[3];  
  my $previousTable = $_[5];   
  my $programTag = $_[9];
  my $sqlQueryOrderBy; 
  my $sqlQuerySpecific;
  my $dispatchTable;
  my $connection = ConnectToMySql( $_[4], $_[6], $_[7], $_[8] );
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
       $sqlQueryOrderBy = "SELECT DISTINCT($thisRecursionFeaturedColumn) FROM $previousTable WHERE `run_id` = '$programTag' ORDER BY `$thisRecursionFeaturedColumn`;";
       $sqlQuerySpecific = "SELECT * FROM `specifictable` WHERE `run_id` = '$programTag' AND `$thisRecursionFeaturedColumn` = ? ";
      },
    1 => sub {
       $thisTempTable = "lessstrict_1";
       $thisRecursionFeaturedColumn = 'alignment_length';
       $sqlQueryOrderBy = "SELECT DISTINCT($thisRecursionFeaturedColumn) FROM $previousTable t WHERE t.`run_id` = '$programTag' ORDER BY t.$thisRecursionFeaturedColumn DESC; ";
       $sqlQuerySpecific = $previousQuerySpecific . " AND $thisRecursionFeaturedColumn = ? "; 
    },
    2 => sub {
       $thisTempTable = "lessstrict_2";
       $thisRecursionFeaturedColumn = 'percent_identity';
       $sqlQueryOrderBy = "SELECT DISTINCT($thisRecursionFeaturedColumn) FROM $previousTable u  WHERE u.`run_id` = '$programTag' ORDER BY u.$thisRecursionFeaturedColumn DESC; ";
       $sqlQuerySpecific =  $previousQuerySpecific . " AND $thisRecursionFeaturedColumn = ? ";
    },
    3 => sub {
       $thisTempTable = "lessstrict_3";
       $thisRecursionFeaturedColumn = 'bitscore';
       $sqlQueryOrderBy = "SELECT DISTINCT($thisRecursionFeaturedColumn) FROM $previousTable v  WHERE v.`run_id` = '$programTag' ORDER BY v.$thisRecursionFeaturedColumn DESC; ";
       $sqlQuerySpecific =  $previousQuerySpecific . " AND $thisRecursionFeaturedColumn = ? ";
    },
    4 => sub {
       $thisTempTable = "lessstrict_4";
       $thisRecursionFeaturedColumn = 'mismatch';
       $sqlQueryOrderBy =  "SELECT DISTINCT($thisRecursionFeaturedColumn) FROM $previousTable w  WHERE w.`run_id` = '$programTag' ORDER BY w.$thisRecursionFeaturedColumn; ";
       $sqlQuerySpecific = $previousQuerySpecific . " AND $thisRecursionFeaturedColumn = ? ";
    },
    5 => sub {
       $thisTempTable = "lessstrict_5";
       $thisRecursionFeaturedColumn = 'gap';
       $sqlQueryOrderBy =  "SELECT DISTINCT($thisRecursionFeaturedColumn) FROM $previousTable x  WHERE x.`run_id` = '$programTag' ORDER BY x.$thisRecursionFeaturedColumn; ";
       $sqlQuerySpecific = $previousQuerySpecific . " AND $thisRecursionFeaturedColumn  = ? ";
    },
    6 => sub {
      # Actually this will not be called, but keeping it here for legacy 
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
    print "[fatal] At feature $geneName recursionLevel $recursionLevel No rows at ORDER by part \n";
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
      $returnThis = searchBestBLAST6Match(
        $recursionLevel + 1,
        $geneName,
        \@previousValues,
        $sqlQuerySpecific,
        $_[4],
        $thisTempTable,
        $_[6],
        $_[7],
        $_[8],
        $programTag
      );
    }else{
      # STill there are more results in the topmost criteria, return the top most
      my $finalRef = $statement3->fetchrow_hashref();
      if ( $debug > 1) { 
        handle_message('WARNING', 'There are still more than one results in recursionLevel 5/max criterion. Returning ' . $finalRef->{ 'id' } . " for $geneName", '');
      }
      $returnThis = $finalRef->{ 'id' };
    }
  }
  # TRUNCATE ANY TABLE TEMPORARY
  my $statement5 = $connection->prepare("DELETE FROM $thisTempTable WHERE `run_id` = '$programTag';");
  $statement5->execute()  || print "$DBI::errstr \n";;
  return $returnThis;
}#sub

