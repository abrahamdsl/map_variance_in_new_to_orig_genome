=doc
  a.llave@irri.org 02OCT2014-0930

  A Perl module for general needs, like printing/output handling.
=cut

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(handle_message);

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

