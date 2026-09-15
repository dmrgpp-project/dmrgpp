#!/usr/bin/perl
# SPDX-FileCopyrightText: Copyright (c) 2009-2015, UT-Battelle, LLC
# SPDX-License-Identifier: GPL-3.0-or-later
#
# [Lanczos++, Version 1.0]

use warnings;
use strict;

use Getopt::Long qw(:config no_ignore_case);
use lib "../../PsimagLite/scripts";
use NewMake;
use lib ".";
use PsiTag;

my ($flavor, $lto) = (NewMake::noFlavor() , 0);
my $usage = "USAGE: $0 [-f flavor] [-lto] [-c config]\n";
my $config;

GetOptions('f=s' => \$flavor,
           'lto' => \$lto,
           'c=s' => \$config) or die "$usage\n";

my $gccdash = "";
if ($lto == 1) {
	$gccdash = "gcc-";
	$lto = "-flto";
} else {
	$lto = "";
}

my @configFiles = ("../TestSuite/inputs/ConfigBase.psiTag");
push @configFiles, $config if (defined($config));

my %thermalD = (name => "thermal", dotos => "thermal.o");
my %lorentzianD = (name => "lorentzian", dotos => "lorentzian.o");
my %dynamics1 = (name => "dynamics1", dotos => "dynamics1.o ProgramGlobals.o LanczosDriver0.o");
my %quasiparticleWeightZ = (name => "quasiparticleWeightZ", dotos => "quasiparticleWeightZ.o ProgramGlobals.o LanczosDriver0.o");
my %edD = (name => "ed", dotos => "ed.o ProgramGlobals.o LanczosDriver0.o");
my %ld0 = (name => "LanczosDriver0", aux => 1);
my %ld1 = (name => "LanczosDriver1", aux => 1);
my %ld2 = (name => "LanczosDriver2", aux => 1);
my %ld3 = (name => "LanczosDriver3", aux => 1);
my %ld4 = (name => "LanczosDriver4", aux => 1);
my %ld5 = (name => "LanczosDriver5", aux => 1);
my %progGlob = (name => "ProgramGlobals", aux => 1);
my $dotos = "lanczos.o LanczosDriver0.o LanczosDriver1.o LanczosDriver2.o";
$dotos .= " LanczosDriver3.o LanczosDriver4.o LanczosDriver5.o ProgramGlobals.o";
my %lanczosD = (name => "lanczos", dotos => $dotos);

my @drivers = (\%thermalD,
	\%lorentzianD,
	\%dynamics1,
	\%quasiparticleWeightZ,
	\%edD,
	\%ld0,
	\%ld1,
	\%ld2,
	\%ld3,
	\%ld4,
	\%ld5,
	\%progGlob,
	\%lanczosD);

my %args;
$args{"CPPFLAGS"} = $lto;
$args{"LDFLAGS"} = $lto;
$args{"flavor"} = $flavor;
$args{"code"} = "Lanczos++";
$args{"configFiles"} = \@configFiles;

createMakefile(\@drivers, \%args);

sub createMakefile
{
	my ($drivers, $args) = @_;
	NewMake::backupMakefile();

	my $fh;
	open($fh, ">", "Makefile") or die "Cannot open Makefile for writing: $!\n";

	NewMake::main($fh, $args, \@drivers);

	close($fh);
	print STDERR "$0: File Makefile has been written\n";
}

