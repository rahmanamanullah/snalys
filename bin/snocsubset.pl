#!/usr/bin/perl
#
# Extract a subset of events, satisfying given conditions, from an output
# file of SNOC.
#
# Rahman Amanullah (rahman@physto.se), 2001-12-17
#
# HISTORY
#   2005-02-16 Added the possibility to give a list of possible
#              acceptable strings to the 'EQ'.
#

#
# The predefined limits that will be used of none are
# specified
#
$minLimit = -1.E30;
$maxLimit = 1.E30;
$ll = 'minLim';
$ul = 'maxLim';
$lls = 'LL';
$uls = 'UL';
$eqs = 'EQ';

&help unless @ARGV;
$file = shift @ARGV;

open ( FP, $file ) ||
  die "Unable to open $file: $!\n";

#
# Perhaps not the fanciest way to read a text file, but
# it works.
#
my $buffer;
while ( $tmp = <FP> ) { $buffer = $buffer.$tmp; }

my @all = split( /event/ , $buffer ); # each event are separated by 'event'
&help unless @all;
my $header = shift( @all );           # the first item is the header
for ( $n = 0; $n <= $#all; $n++ ) { $all[$n] = 'event'.$all[$n]; }

#
# Build a list of all the keywords in the the data file.
#
my @firstevent = split( /\n/ , $all[0] );
my @DataKeys = undef;
foreach $line ( @firstevent ) {
  my ( $key , $value ) = split( /[ \t\n]+/ , $line );
  push @DataKeys , $key;
}
pop @DataKeys;                        # the 'end' is not a valid key
$DataKeys = join ( " " , @DataKeys);

#
# Take care of the rest of the argument list to
# see which keywords that should be considered.
#
&help unless ( ($#ARGV + 1) % 3 == 0 );
my $nrarg = $#ARGV;
for ( $n = 0; $n <= $nrarg; $n += 3 ) {
  my $key = shift @ARGV;
  &help unless $DataKeys =~ /$key/;
  &addKey($key);
  my $first = shift @ARGV;
  my $second = shift @ARGV;

  #
  # Build up the hash with conditioning keywords.
  #
  if ( $first eq $lls && &IsNr($second) ) { &minLim( $key , $second  ); }
  elsif ( $first eq $uls && &IsNr($second)) { &maxLim( $key, $second ); }
  elsif ( $first eq $eqs ) { &Limits( $key, $second, $second ) }
  elsif ( &IsNr($first) && &IsNr($second) ) {
    if ( $second < $first ) { ( $first , $second ) = ( $second , $first ); }
    &Limits( $key, $first, $second );
  }
  else { &help;}
}

#
# Sort out the events that should be kept.
#
my @kept;
foreach $event ( @all ) {
  my $keep = 1;
  my @entries = split( /\n/ , $event );
  
  #
  # Store the values for each entry in a hash table.
  #
  my %sn = undef;
  foreach my $entry ( @entries ) {
    my ( $key , $value ) = split ( /[ \t\n]+/ , $entry );
    $sn{$key} = $value;
  }

  #
  # Run through the keyword conditions and check them.
  #
  foreach my $key ( keys %{&GetHash} ) {
    my $llim = &minLim($key);
    my $ulim = &maxLim($key);

    # If this is text we might have given a list of alternative
    # acceptable values.
    #
    if ( (! &IsNr($llim)) ) {
      $keep = 0;
      my @strings = split /[ \t\n]+/ , $llim;
      foreach my $string ( @strings ) {
	$keep = 1 if $sn{$key} eq $string;
      }
    }
    elsif ( ! ( $sn{$key} >= $llim &&
		$sn{$key} <= $ulim ) ) { $keep = 0; }
  }

  #
  # Keep the event?
  #
  if ( $keep ) { push @kept, $event; }
}

#
# Change the number of events in the header.
#
my $nrall = $#all+1;
my $nrkept = $#kept+1;
$header =~ s/NUMSNE      $nrall/NUMSNE      $nrkept/;

#
# Print the kept events
#
print STDOUT $header, @kept if $nrkept > 0;


#
# The help (and exit) routine.
#
sub help {
print STDOUT <<EOHelp;

  snocsubset.pl -- extracts a subset of a SNOC output file that
                   fulfills some given conditions.

  Usage
		     
     snocsubset.pl <snocfile.out> [<key> <op/val> <val> [..]]

  Description
    
     The script extracts a subset of events from a SNOC output file
     that satisfies some given conditions and redirects these to
     standard output.
    
     A condition is defined by three words, the keyword to which the
     condition will be applied, an operator or a value, and a value.
     The operator could be 'UL', 'LL' or 'EQ'. For example, 'UL' means
     that only supernovae whos <key> are less than <val> should be kept.

     The keyword 'EQ' is only recommended to be used when <val> is a
     string though it should work also in the numerical case. The other
     operators only work for numerical values.

     Instead of an operator, it is possible to give two values defining
     an interval in which the value of <key> must lie for the event to
     be kept.

  Example

     snocsubset.pl file.out linsdm UL 0.01    # only keep low lensed SNe
     snocsubset.pl file.out sntype EQ Ia      # only interested in Ia:s
     snocsubset.pl file.out sntype EQ "Ia Ic" # only keep Ia and Ic events
     snocsubset.pl file.out zs 0.0 2.0        # defines a redshift range
     snocsubset.pl file.out zs 0.0 2.0 sntype EQ Ia
       
EOHelp
  
  exit 1;
}

#
# Boolean check if a variable is a number
#
sub IsNr {
  my $bad = 0;
  local $SIG{__WARN__} = sub { $bad++ };
  local $^W = 1;
  my $guess = shift;
  $guess += 0;
  return not $bad;
}

#########################################################
#                                                       #
# Structure to handle the keyword condition for keeping #
# an event.                                             #
#                                                       #
#########################################################

#
# Add a key to the hash table.
#
sub addKey {
  my $key = shift;
  my $keys = &GetHash;
  my %limits = ( $ll, $minLimit ,
		 $ul, $maxLimit );
  $keys->{$key} = \%limits;
}

#
# Add or get a lower limit for a key.
#
sub minLim {
  return unless @_;
  my $key = shift;
  my $keys = &GetHash;
  if ( @_ ) { $keys->{$key}->{$ll} = shift; }
  return $keys->{$key}->{$ll};
}

#
# Add or get an upper limit for a key.
#
sub maxLim {
  return unless @_;
  my $key = shift;
  my $keys = &GetHash;
  if ( @_ ) { $keys->{$key}->{$ul} = shift; }
  return $keys->{$key}->{$ul};
}

#
# Add or get the boundary interval for a key.
#
sub Limits {
  return unless @_;
  my $key = shift;
  my $keys = &GetHash;
  if ( scalar(@_) > 1 ) {
    $keys->{$key}->{$ll} = shift;
    $keys->{$key}->{$ul} = shift;
  }
  return ( $keys->{$key}->{$ll},  $keys->{$key}->{$ul});
}

#
# Return the hash table. 
#
sub GetHash { return \%KeyWords; }

# EOF
