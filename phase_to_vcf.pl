
###########################
# Code by Erick C. Castelli 
# www.castelli-lab.net
# January, 2015
# Version 2.5
# This Perl script converts a PHASE .out_pairs file fo a VCF file 
###########################


use strict;
use Getopt::Std;

our ($opt_i, $opt_v, $opt_o);
getopts('i:v:o:');

if ($opt_i eq "") {help();}
if ($opt_v eq "") {help();}
if ($opt_o eq "") {help();}

if (! -e $opt_i) 
{
	print "PHASE input file not found or not indicated!\n";help();
	exit;
}
if (! -e $opt_v) 
{
	print "VCF input file not found or not indicated!\n";help();
	exit;
}


open (IN, $opt_i);
my %hash = ();
my $sample;
my $prob;
my $new_prob;
my @samples = ();
my $initial_index;

while (<IN>)
{

chomp $_;
$_ =~ s/ , /,/g;

if ($_ =~ /IND:/)
{
	$sample = $_;
	$sample =~ s/IND: //;
	$prob = 0;
	next;	
}

my @values = split (",",$_);
$new_prob = $values[2];


if ($new_prob > $prob) {
	$prob = $new_prob;
	$values[0] =~ s/ //g;
	$values[1] =~ s/ //g;
	$hash{$sample}{1} = $values[0];
	$hash{$sample}{2} = $values[1];
	$hash{$sample}{"prob"} = $values[2];
	next;
}	

}
close (IN);





open (IN, $opt_v);
open (OUT, ">$opt_o");


my $index = 0;

while (<IN>)
{

	chomp $_;
	chomp $_;
	
	if ($_ =~ /\#\#/) {
		if ($_ =~ /fileformat/) {print OUT $_ . "\n";}
		if ($_ =~ /ID=GT/) {print OUT $_ . "\n";}
		if ($_ =~ /contig=/) {print OUT $_ . "\n";}
		if ($_ =~ /reference=file/) {print OUT $_ . "\n";}
		next;
		}


	if (substr($_,0,1) eq "#") 
	{
		print OUT $_ . "\n";
		
		my @temp = split ("\t",$_);
		for my $line (@temp)
		{
			if ($line =~ m/(CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT)/) {
				$initial_index++;
				next;
				}
			push (@samples,$line);
		}
		next;
	}
	

	

	my @temp = split ("\t",$_);
	my @slice = @temp[0..4];
	print OUT join("\t",@slice);
	print OUT "\t";		

	print OUT ".\t.\t.\tGT";

	for my $sample_name (@samples)
	{
		my @temp1 = split("",$hash{$sample_name}{1});
		my @temp2 = split("",$hash{$sample_name}{2});

				
			if ($temp1[$index] eq "") {$temp1[$index] = ".";}
			if ($temp2[$index] eq "") {$temp2[$index] = ".";}
		

			print OUT "\t" . "$temp1[$index]|$temp2[$index]";

	}
	print OUT "\n";
	$index++;
}

sub help
{
	print "-i [the PHASE out_pairs file]\n-v [the VCF to be used as a draft]\n-o [the output VCF file]\n\n";
	exit;
}

