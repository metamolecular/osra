#!/usr/bin/perl
#
# Simple test script that traverses the given folder and rund OSRA for every TIFF file found.
#

use strict;

use File::Find;
use File::Path;
use File::Basename;
use File::Spec::Functions qw(abs2rel);

if ($#ARGV != 0) {
	print "Usage: $0 <dir>\n";
	exit 0;
}

my $test_dir = $ARGV[0];
my $out_dir = 'run5';

local $| = 1;
local $/;

find({ no_chdir => 1, wanted => sub {
	return if -d $_ || (-f $_ && !/\.tif/);
	
	my $location = $out_dir . '/' . abs2rel($File::Find::dir, $test_dir);
	
	mkpath($location) or die "Can't create the directory $location: $?" unless -d $location;
	
	$location .= '/' . basename($_) . '.out';

	print "process $_ --> $location\n";
	
	open IN, "../src/osra -c -p -b -f sdf '$_' |" or die;
	my $mol = <IN>;
	close IN;
	
	open OUT, ">", $location or die; 
	print OUT $mol;
	close OUT;
} }, $test_dir);
