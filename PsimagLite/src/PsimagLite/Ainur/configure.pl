#!/usr/bin/perl
# SPDX-FileCopyrightText: Copyright (c) 2009-2017, UT-Battelle, LLC
# SPDX-License-Identifier: GPL-3.0-or-later
#
# [PsimagLite, Version 1.]

use warnings;
use strict;

use lib "../../../PsimagLite/scripts";
use Make;

my @drivers = ();

createMakefile();

sub createMakefile
{
	Make::backupMakefile();
	if (!(-r "Config.make")) {
		my $cmd = "cp Config.make.sample Config.make";
		system($cmd);
		print STDERR "$0: Executed $cmd\n";
	}

	my $fh;
	open($fh, ">", "Makefile") or die "Cannot open Makefile for writing: $!\n";

	local *FH = $fh;
	my @units = qw(test testReal);
	my $combinedUnits = combine("",\@units,".o ");
	my $combinedUnits2 = combine("./",\@units,".cpp ");

	print FH<<EOF;
include Config.make
all: @units
EOF

	foreach my $unit (@units) {
		my $doth = "../src/".ucfirst($unit).".h";
		my $tmp = (-r "$doth") ? "$doth" : "";
		print FH<<EOF;
$unit.o: ./$unit.cpp $tmp Makefile Makefile.dep
\t\$(CXX) \$(CPPFLAGS) -c -I../src -I..  ./$unit.cpp

$unit: $unit.o $tmp Makefile Makefile.dep
\t\$(CXX) -o $unit $unit.o \$(LDFLAGS)
EOF
	}

print FH<<EOF;
Makefile.dep: $combinedUnits2
\t\$(CXX) \$(CPPFLAGS) -I../src -I.. -MM  $combinedUnits2  > Makefile.dep

clean: Makefile.dep
\trm -f core* *.o *.dep *.a @units

include Makefile.dep

EOF

	close($fh);
	print STDERR "File Makefile has been written\n";
}

sub combine
{
	my ($pre,$a,$post) = @_;
	my $n = scalar(@$a);
	my $buffer = "";
	for (my $i = 0; $i < $n; ++$i) {
		$buffer .= $pre.$a->[$i].$post;
	}

	return $buffer;
}

