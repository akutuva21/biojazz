use strict;
use warnings;
use lib "/tmp/anc/base";
use lib "modules";

package Named;
sub new { bless {}, shift }

package main;
use Scoring;

print "Successfully loaded Scoring.pm\n";
