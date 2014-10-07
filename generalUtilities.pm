=doc
  a.llave@irri.org 02OCT2014-0930

  A Perl module for general needs, like printing/output handling.
=cut

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(handle_message);

sub ConnectToMySql {
=doc
  Used for connecting to MySQL database.
  Disclaimer: This function was sourced from somewhere else, I forgot where. All credits go to
  original poster.
  Arguments
  0 - database name
  Returns
  DBI object for MySQL connection.
=cut
  my ($db) = @_;

  # assign the values in the accessDB file to the variables
  my $host = "localhost";
  my $userid = "root";
  my $passwd = "e\$nodenpri\$m";
  # assign the values to your connection variable
  my $connectionInfo="dbi:mysql:$db;$host";
  # make connection to database
  my $l_connection = DBI->connect($connectionInfo,$userid,$passwd);
  # the value of this connection is returned by the sub-routine
  return $l_connection;
} # Sub

sub getFileNumLines {
=doc
  @devnote Sourced from https://raw.githubusercontent.com/abrahamdsl/blast_and_infer_subfeatures/master/blast_and_infer_genes_and_subfeatures.pl

  Arguments:
    0 - string. The filename of the file.
 
  Returns:
    int. The number of lines in a file. 
=cut
  my $lines = 0;
  my $fileHandle;
  open ( $fileHandle, $_[0] ) or die('ERROR',"Can't open $_[0] for line number checking!");
  $lines++ while( <$fileHandle> );
  close $fileHandle;
  return $lines;
} # getFileNumLines

sub handle_message {
=doc
  Our unified function for printing messages.

  Arguments:
    0 - string. Message level. Could be:
      NOTICE, DEBUG, BENCHMARK, ERROR, FATAL
      The last one causes program to halt.
    1 - string. Short description, maybe a word or two regarding the message.
    2 - string. The full text of the message.

  Returns:
    Nothing
=cut
  my ($level, $code, $sentMsg) = @_;
  $sentMsg ||= "No message.";
  chomp $sentMsg;
  my $message .= ( localtime() . ' ' . $sentMsg . "\n" );
  if ($level eq 'FATAL') {
    die join( ' : ', ( "[$level]", $code, $message) );
  }else {
    print STDERR join( ' : ', ( "[$level]", $code, $message) );
  }
} # handle_message

