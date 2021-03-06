# EXOME SEQUENCING ALIGNMENT: CONTROL SAMPLES

## Background:

I need some exome sequenceing data from control samples to determine ethnicity of our exome data and for burden tests etc. The 1000 Genomes project has freely available exome fastq data from different ethnic populations, which I can download and analyse. Analysing data from fastq files, rather than per-prepared variant lists, will allow a better comparison with our data, as the analysis methods are identical.


### Step 1. Find appropriate samples

I need exome fastq files that are comparable to our data, i.e.:

* Sequenced on the Illumina HiSeq
* Paired end
* >50 million reads
* Suitable exome targets (e.g. BGI used an old exome enrichment kit called SureSelect V1, which does not capture RYR2)
* Single large files, as opposed to smaller ‘pond’ library files

Go to ncbi sra > search > [SRA object](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=search_obj)

Search for ‘1000 Genomes Exome’ and retrieve the Accession number (e.g. SRP004060)

Go to ncbi > sra run selector > [search for SRP004060](http://trace.ncbi.nlm.nih.gov/Traces/study/?go=home)

On left hand column (facets), select 

* Centre >bcm, bi
* Platform > illumina
* MBases > 5000+

Select + to select all runs, then download RunInfo Table. Now go to [EBI ena](http://www.ebi.ac.uk/ena) and search for the Run_s (e.g. SRR765989) view > text > copy the ftp link. I did not download BGI data as they used an old exome enrichment design that is small and, for example, does not capture the RYR2 gene. I did not download  WUGSC data as they used the illumine GAIIx, and I’m not sure how well that will compare to the HiSeq. Anything less than 5000 MBases is verging on too shallow coverage compared to my data. If I want to increase my EUR samples, I could test the WUGSC data, and download fastq data with less than 5000 MBases, if there are multiple files available per sample. 

### Step 2. Download samples

Use aspera to download samples. Put the ftp links in a ftp_link_file.txt file and use Rob Middleton’s `aspera-download.pl ftp_link_file.txt` script. 

### Step 3. Align and calculate the coverage metrics

Analyse samples in batches, according to ethnic super population (AFR, EAS, EUR, SAS) and sequencing centre (BCM use vcrome2.1; BI use sureselectv2). Use the Control_exome_alignment.pl script, which will

* Align with bwa
* Handle the space in the fastq read header lines
* Clip the sam file @PG line for GATK compatability
* Stream sorting with novosort
* Calculate per base coverage
* Plot coverage and determine sex of samples

See the [ExomeCoverage script](https://github.com/rbagnall/ExomeCoverage)
