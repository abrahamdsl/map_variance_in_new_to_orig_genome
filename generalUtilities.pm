=doc
  a.llave@irri.org 02OCT2014-0930

  A Perl module for general needs, like printing/output handling.
=cut

require Exporter;
@ISA = qw(Exporter);
#@EXPORT;

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

