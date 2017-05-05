
###########################
# Code by Erick C. Castelli 
# www.castelli-lab.net
# January, 2015
# Version 2.5
# This script uses the phased data from GATK readbackedphasing VCF file
#      to create a fragmented .known file for the PHASE algorithm
#      and compares PHASE results from multiple runs 
#
# Dependencies: PHASE (http://stephenslab.uchicago.edu/software.html#phase)
#               VCFX (www.castelli-lab.net/apps/vcfx)
#
###########################

use strict;
use threads;
use Thread::Semaphore;
use Getopt::Std;
use FindBin '$Bin';

#***** System Variables ***************************************
my $program = "phase-using-known-from-readbackedphasing";
my $version = "2.5";
my $day;
my $month;
my $year;
my $hour;
my $min;
my $seg;
($day, $month,$year,$hour,$min,$seg) = (localtime)[3,4,5,2,1,0];
$year = $year+1900;
my $os = $^O;
#****************************************************************

my $opt_out = "$Bin/output/";
my $max_threads = 0;
our ($max_segments);

our ($opt_o, $opt_t, $opt_v, $opt_d, $opt_b, $opt_p, $opt_x);

$opt_x = `which vcfx`;
chomp $opt_x;

$opt_b = `which phase`;
chomp $opt_b;

getopts('o:t:v:d:b:p:x:');


my $model = "";
if ($opt_d ne "") {$model = $opt_d;}

my $phase = $opt_b;
if (! -e $phase) {$opt_b = ""; print "\nThe PHASE algorithm was not detected or not found"; help();}

my $vcfx = $opt_x;
if (! -e $vcfx) {$opt_x = ""; print "\nVCFx was not detected or not found"; help();}


my $phase_para = "";
if ($opt_p ne "") {$phase_para = $opt_p;}

$max_threads = $opt_t;
if ($opt_t eq "") {$max_threads = 4;};

if (! -e $opt_v) {print "\nVCF file not found"; help();}




# creating output structure 
	my $folder;
	if ($opt_o eq ""){
		$folder = "$program-" . $day . "_" . ($month+1) . "_" . $year . "-" . $hour . "_" . $min . "_" . $seg;
	}
	else {
		$folder = $opt_o;
	}
		mkdir($folder);
#*********************************************************

$folder = $folder . "/";



open (LOG, ">$folder" . "command.txt");
	print LOG $0 . ", version $version\n";
	print LOG "input file: " . $opt_v . "\n";
	print LOG "output folder: " . $folder. "\n";
	print LOG "Parameters used " . "\n";
	print LOG "     iteraction, thinning, burning: $phase_para" . "\n";
	print LOG "     others: $model" . "\n";
	
close (LOG);


print "\n";
print "$program (version $version, using PHASE version 2.1)\n";

print "Creating output folder at $folder\n";

print "Exporting VCF to .INP using VCFX...\n";
my $cmd = "$vcfx phase input=$opt_v output=$folder" . "/original --quiet --noupdate";
system($cmd);
my $opt_i = $folder . "/original.inp";


print "Creating fragmented KNOWN from VCF...\n";
create_original_known($opt_v, $folder);
my $opt_k = $folder . "/original.known";

print "Creating KNOWN file for each PHASE run...\n";
create_known_file($opt_k, $folder);

my $cmd = "cp $opt_v $folder/original.vcf";
system ($cmd);


my $sem = Thread::Semaphore->new($max_threads); # max 15 threads
my $current_segment = 0;

my @threads = map {
	# request a thread slot, waiting if none are available:
	$current_segment++; 
	$sem->down;
    
	if ($current_segment eq ($max_segments + 1))
	{
		threads->create(\&run_phase_unknown, $phase, $opt_i, $folder, $phase_para, $model)
	}
	else {
		threads->create(\&run_phase, $phase, $opt_i, $current_segment, $folder, $max_segments, $phase_para, $model)
	}
	
	
} 1..$max_segments+1;
$_->join for @threads;





my %data = ();
my $sample_index = 0;




print "Comparing data...\n";
for (my $a = 1; $a <= $max_segments; $a++)
{
	open (IN, $folder . "phase_" . sprintf("%04d",$a) . ".out_pairs");
	$sample_index = 0;
	my $sample_name = "";
	
	while (<IN>)
	{
		chomp $_;
		if ($_ =~ /IND:/)
		{
			$sample_name = $_;
			$sample_name =~ s/IND://;
			$sample_index++;
			next;
		}
		
		my @temp = split (" , ",$_);
		$data{$sample_index}{$sample_name}{$temp[0] . "," . $temp[1]}{$a} = $temp[2];
		$data{$sample_index}{$sample_name}{$temp[1] . "," . $temp[0]}{$a} = $temp[2];
		
	}
	close (IN);
}


open (IN, $folder . "phase_unknown.out_pairs");
$sample_index = 0;
my $sample_name = "";
	
while (<IN>)
{
	chomp $_;
	if ($_ =~ /IND:/)
	{
		$sample_name = $_;
		$sample_name =~ s/IND://;
		$sample_index++;
		next;
	}
	
	my @temp = split (" , ",$_);
	$data{$sample_index}{$sample_name}{$temp[0] . "," . $temp[1]}{"unknown"} = $temp[2];
	$data{$sample_index}{$sample_name}{$temp[1] . "," . $temp[0]}{"unknown"} = $temp[2];
	
}
close (IN);



open (OUT, ">$folder" . "haplotypes.csv");
print OUT "sep=,\n";

print OUT "sample,haplotype_1,haplotype_2";
for (my $b = 1; $b <= $max_segments; $b++)
{
	print OUT ",run_segment_$b";
}
print OUT ",run_unknown_phase,compatibility\n";



for ($a = 1; $a <= $sample_index; $a++)
{
	my @sample_name = keys ($data{$a});
	print OUT "\n";
	
	my @haplos = keys ($data{$a}{$sample_name[0]});
	my %haplo_duplicate = ();
	foreach (@haplos)
	{
		
		if ($haplo_duplicate{$_} eq 1) {next;}
		
		my @haplo_temp = split(",", $_);
		$haplo_duplicate{$haplo_temp[0] . "," . $haplo_temp[1]} = 1;
		$haplo_duplicate{$haplo_temp[1] . "," . $haplo_temp[0]} = 1;
		
		print OUT $sample_name[0] . "," . $_ . ",";
		for (my $b = 1; $b <= $max_segments; $b++)
		{
			if ($data{$a}{$sample_name[0]}{$_}{$b} eq "")
			{
				print OUT "not detected" . ",";
			}
			else {
				print OUT $data{$a}{$sample_name[0]}{$_}{$b} . ",";
			}
			
		}
		
		if ($data{$a}{$sample_name[0]}{$_}{"unknown"} eq "")
		{
			print OUT "not detected" . ",";
		}
		else {
			print OUT $data{$a}{$sample_name[0]}{$_}{"unknown"} . ",";
		}

		
		print OUT "\n";
	}
	
}
close (OUT);





open (OUT, ">$folder" . "haplotypes.checked.csv");
my $input_haps = "$folder" . "haplotypes.csv";

my %higher_p;
my %higher_h;

open (IN, $input_haps);
while (<IN>)
{
	chomp $_;
	if ($_ =~ /sep/) {next;}
	if ($_ =~ /compatib/) {next;}
	if ($_ eq "") {next;}
	
	my @temp = split(",",$_);

	for (my $a = 3; $a < scalar(@temp); $a++)
	{
		
		if ($temp[$a] eq "not detected") {$temp[$a] = 0;}
	
		if ($temp[$a] > $higher_p{$temp[0]}{$a})
		{
			$higher_p{$temp[0]}{$a} = $temp[$a];
			$higher_h{$temp[0]}{$a} = $temp[1] . "," . $temp[2];
		}
	}
	
}
close (IN);




open (IN, $input_haps);
while (<IN>)
{
	chomp $_;
	if ($_ =~ /sep=/) {print OUT $_ . "\n";next;}
	if ($_ =~ /compatib/) {print OUT $_ . "\n";next;}
	if ($_ eq "") {print OUT $_ . "\n"; next;}
	my @temp = split(",",$_);
	
	my $helper = 0;
	for (my $a = 3; $a < scalar(@temp); $a++)
	{
		my $h = $temp[1] . "," . $temp[2];
		if ($higher_h{$temp[0]}{$a} ne $h) {$helper = 1;}
	}
	
	if ($helper eq 0) {print OUT $_ . "ok\n";}
	else {
		print OUT $_ . "no\n";
	}	
	
	
}
close (IN);

close (OUT);









print "Checking PHASE log...\n";
for (my $a = 1; $a <= $max_segments; $a++)
{
	
	open (IN, $folder . "log_run_$a.txt");
	while (<IN>)
	{
		chomp $_;
		if ($_ =~ /Warning/)
		{
			if ($_ !~ /genotype missing/)
			{
				print "     run $a: " . $_ . "\n";
			}
		}
	}
	close (IN);
}
open (IN, $folder . "log_run_unknown.txt");
while (<IN>)
{
	chomp $_;
	if ($_ =~ /Warning/)
	{
		if ($_ !~ /genotype missing/)
		{
			print "     run phase unknown: " . $_ . "\n";
		}
	}
}
close (IN);

print "Computation done! Please check file $folder/haplotypes.checked.csv\n\n";

exit;





sub run_phase {
	 	my ($phase,$i,$run,$out,$max,$phase_para,$model) = @_;
		
		$model = $model . " ";
	
		print "Starting PHASE run for segment number $run...\n";

		my $known_file = $out . "known.$run.known";
		my $output = $out . "phase_" . sprintf("%04d",$run) . ".out";
		
		my $command = "";
		$command = "$phase -k$known_file $model-F0.01 -O0.01 $i $output $phase_para >$out" . "log_run_$run.txt 2>&1";
		#print $command . "\n";
		
		system ($command);
		$sem->up;
	}





sub run_phase_unknown {
	 	my ($phase,$i,$out,$phase_para,$model) = @_;

		$model = $model . " ";
		print "Starting PHASE run ignoring the known file...\n";
		my $output = $out . "phase_unknown.out";
					
		my $command = "";
		$command = "$phase $model-F0.01 -O0.01 $i $output $phase_para 1>$out" . "log_run_unknown.txt 2>&1";
		system ($command);
		$sem->up;
	
}




sub help {
	print "\n\n";
	print "$program (version $version, using PHASE version 2.1)\n";
	print "    -o [output_folder] (optional)\n";
	print "    -t [number of threads] (optional, default: 4)\n";
	print "    -v [input VCF file phased by GATK ReadBackedPhasing] (mandatory)\n";
	print "    -b [path do the PHASE binary] ( ";
		if ($opt_b eq "") {print "not detected, mandatory)\n";} else {print "located at " . $opt_b . " )\n";}
	

	print "    -x [path do the VCFX binary] ( ";
		if ($opt_x eq "") {print "not detected, mandatory)\n";} else {print "located at " . $opt_x . " )\n";}

	print "    -p [iteraction, thinning and burning parameters] (e.g. -p '1000 1 100')\n";
	print "    -d [other parameters] (such as -d '-N300 -l12 -d1')\n";

	print "    example: perl $0 -t 6 -v file.vcf -b /usr/local/bin/phase -x /usr/local/bin/vcfx\n";

	exit;	




}








sub create_known_file {
 	my ($k, $out) = @_;

	$max_segments = 0;
	my %segments = ();
	my $line_number = 1;
	
	open (KNOWN, $k);
	while (<KNOWN>)
	{
		chomp $_;
		my @data = split("\:",$_);
		$_ =~ s/$data[0]//g;
		$_ =~ s/\|/,/g;
		$_ =~ s/\://g;
		my @temp = split(",",$_);
		if (scalar(@temp) > $max_segments) {$max_segments = scalar(@temp);}
		$segments{$line_number} = scalar(@temp);
		$line_number++;
	}	
	close (KNOWN);


	for (my $a = 1; $a <= $max_segments; $a++)
	{
		open (OUT, ">" . $out . "/known.$a.known");
		close (OUT);
	}


	my $line_number = 1;
	open (KNOWN, $k);
	while (<KNOWN>)
	{
		chomp $_;
		my $line = $_;
		
		my @data = split("\:",$line);
		$line =~ s/$data[0]//g;
		$line =~ s/\://g;
				
		for ($a = 1; $a <= $max_segments; $a++)
		{
			open (OUT, ">>$out/known.$a.known");

			if ($segments{$line_number} eq 1) 
			{
				print OUT $line . "\n";
			}
			
			else {
				$line =~ s/\|/,/g;
				my @temp = split(",",$line);
				my $current_number_of_segments = scalar(@temp);
			
				if ($a <= $current_number_of_segments)
				{
					for ($b=1; $b <= $max_segments; $b++)
					{
							if ($a eq $b) {print OUT $temp[$b-1];}  # segmento printado
							if ($a ne $b) {print OUT ("*" x length($temp[$b-1]));}
					}
				}
				
				
				if ($a > $current_number_of_segments)
				{
					my $random = int(rand($current_number_of_segments+1));
					if ($random eq 0) {$random = int(rand($current_number_of_segments+1));}
					if ($random eq 0) {$random = int(rand($current_number_of_segments+1));}
					if ($random eq 0) {$random = int(rand($current_number_of_segments+1));}
					if ($random eq 0) {$random = 1;}
					
					for ($b=1; $b <= $max_segments; $b++)
					{
							if ($random eq $b) {print OUT $temp[$b-1];}  # segmento printado
							if ($random ne $b) {print OUT ("*" x length($temp[$b-1]));}
					}
				}
				
				
				
				print OUT "\n";
			}			
		
			close (OUT);

		}
			
	$line_number++;
	}	


}






sub create_original_known
{
 	my ($k, $out) = @_;

	my @samples = ();
	my $sample_index = 0;
	my %data = ();
	my %known = ();

	open (IN, $k);
	while (<IN>)
	{
		chomp $_;
		if (substr($_,0,2) eq "##") {next;}
		
		if ($_ =~ /CHROM/) 
		{
			@samples = split ("\t", $_);
			foreach (@samples)
			{
				$sample_index++;
				if ($_ eq "FORMAT") {last;}
			}
		next;
		}
		my @line = split("\t", $_);

		for (my $a = $sample_index; $a < scalar(@line); $a++)
		{
			$data{$samples[$a]}{$line[1]} = $line[$a];
		}
		

	}
	close (IN);



	for (my $a = $sample_index; $a < scalar(@samples); $a++)
	{
		my @position = (sort keys $data{$samples[$a]});
		my $segment = 0;
		foreach (@position)
		{
			
			my $position = $_;
			my @snp_data = split(":",$data{$samples[$a]}{$position});
			my @alleles = split("/",$snp_data[0]);
			
			if ($alleles[0] eq $alleles[1]) {$known{$samples[$a]}{$position} = "*";next;}
			if (($alleles[0] eq ".") || ($alleles[1] eq ".")) {$known{$samples[$a]}{$position} = "*";next;}
			if ($snp_data[4] !~ /-/) {$known{$samples[$a]}{$position} = "*";next;}

			if ($snp_data[4] =~ /-/) 
			{
				my @vetor = split (",",$snp_data[4]);
				my @pos = split("-",$vetor[0]);
				my $pos = $pos[0];
				my $fase = $pos[1];
				
				if ($position eq $pos) {$known{$samples[$a]}{$position} = "*";next;}
				if ($position ne $pos)
				{
					if ($segment eq 0) {$segment = $pos;}
					
					if ($segment eq $pos)
					{
						$known{$samples[$a]}{$pos} = "0";
						$known{$samples[$a]}{$position} = ($fase - 1);
					}
					else 
					{
						$known{$samples[$a]}{$pos} = "|0";
						$known{$samples[$a]}{$position} = ($fase - 1);

					}
					next;
				}
				
				
			}

		}
		
	}

	open (OUT, ">$out/original.known");
	for (my $a = $sample_index; $a < scalar(@samples); $a++)
	{
		print OUT $samples[$a] . ":";
		my @position = (sort keys $data{$samples[$a]});
		foreach (@position)
		{
			
			my $position = $_;
			print OUT $known{$samples[$a]}{$position};			

		}
		print OUT "\n";
	}
	close (OUT);


}


