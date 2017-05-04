# phase-using-known-from-readbackedphasing.pl
Perl script
Written by Erick C. Castelli
Version 2.5
Molecular Genetics and Bioinformatics Laboratory, School of Medicine, Unesp, Botucatu-SP, Brazil (www.castelli-lab.net). 

This script uses the phased data from a GATK ReadBackedPhasing VCF file to create a fragmented .known file to be used with the PHASE algorithm.  Then, it compares the PHASE results from multiple runs considering this fragmented .known file.

# References

This methodology was described and used at the following manuscripts. If you use this script, please cite them both.

Castelli EC et al. HLA-G variability and haplotypes detected by massively parallel sequencing procedures in the geographicaly distinct population samples of Brazil and Cyprus. Mol Immunol. 2017 Mar;83:115-126. doi: 10.1016/j.molimm.2017.01.020. 

Lima TH et al. HLA-F coding and regulatory segments variability determined by massively parallel sequencing procedures in a Brazilian population sample. Hum Immunol. 2016 Oct;77(10):841-53. doi: 10.1016/j.humimm.2016.07.231. 

# How to use this script

First, you need working copies of 
- the PHASE algorithm (download it at http://stephenslab.uchicago.edu/software.html#phase) 
- VCFx software (download it at www.castelli-lab.net/apps/vcfx).

To use the script, at your terminal, 'perl phase-using-known-from-readbackedphasing.pl'. This will display a list of options to set the script, as it follows:

- -o [output_folder] (optional)
- -t [number of threads] (optional, default: 4)
- -v [input VCF file phased by GATK ReadBackedPhasing] (mandatory)
- -b [path do the PHASE binary] ( not detected, mandatory)
- -x [path do the VCFX binary] ( located at /usr/local/bin/vcfx )
- -p [iteraction, thinning and burning parameters] (e.g. -p '1000 1 100')
- -d [other parameters] (such as -x '-N300 -l12 -d1')
- example: perl phase-fragmented-known-version_2.5.pl -t 6 -v file.vcf 

Option -o

This is optional. You may indicate the output folder here. If you don't, a new folder with the date and time will be created in the same folder where the script is located.

Option -t

This is optional. You may indicate the number of threads to be used. Each thread perform a parallel PHASE run.

Option -v

This is the VCF file. This file should be phased using the GATK ReadBackedPhasing algorithm.

Option -b

Path to the PHASE binary. If not installed in the system, please indicate the full path for the binary (including the binary name) here.

Option -x

Path to the VCFx binary. If not installed in the system, please indicate the full path for the binary (including the binary name) here.

Option -p

You may set the interaction, thinning and burning parameters here. E.g., -p '1000 1 100'. The default value is -p '100 1 100'

Option -d

You may set additional PHASE parameters here, such as the segment size, mutation model, etc. E.g., -d '-l12. 


# The script output

At the output folder, you will find:
- original.inp (this is the VCF file converted to the .INP file, by using VCFx phase)
- original.known (this is the fragmented .known file, keeping the phase information from ReadBackePhasing). The presence of "|" means that we don't known the phase between these two fragments.
- original.vcf (this is a copy of the VCF file)
- known.@.known files (these files represent the .known that will be used by the PHASE algorithm in each @ run). There will be as @ as the maximum number of fragments observed in the original.known file.
- log_run_@.txt (the log file for each PHASE run)
- phase_@.out (the PHASE outputs for each run)
- haplotypes.checked.csv (This is the final file. In here you will find the haplotype pairs for each sample in each run, and the haplotype pairs when not considering the .known file, as well as the compatibility among runs. All pair of haplotypes marked with "ok" under compatibility indicate that, for that sample, the same pair of haplotypes were detected presenting the highest P-value for all runs. Thus, this pair is also compatible with ReadBackedPhasing information.)


