# GVF2VCF_Convertor_TestBlu
Super fast conversion of GVF file format to VCF file format

In order to accurately identify variants, particularly somatic variants, it is necessary to have knowledge of the germline variants present in the organism being studied. 
This information can be obtained from databases such as NCBI for human genomes, but for non-human organisms, it is only available in a limited form in some databases. For example, the <a href="http://www.ensembl.org/info/data/ftp/index.html/">Ensembl database</a> provides germline variants in the GVF (genome variant format) file format, which differs from the VCF format used by tools such as GATK.

For this reason, Ensembl developers have made an API available to convert GVF files to VCF, but this API has several shortcomings, including:
1. slow convergence times (it took about 6 hours in one tested sample!)
2. potential incompleteness due to multiple users accessing the database
3. limited support for reference genome header patterns (It only supports the header pattern of the Ensembl reference genome) 
4. lack of some variants reported in the GVF file in the output VCF file.

To address these issues, a Python script was developed to perform an offline conversion. This script is:
1. significantly faster (It performs the conversion in 20 minutes (about six gigs of RAM is required))
2. can report significant variants in a separate file (i.e. reported variants in articles or seen in many samples)
3. supports both the UCSC and Ensembl reference genome header patterns. 


In addition, comparison of the output VCF file from this script and the Ensembl API shows that all variants from the Ensembl output are included, as well as variants not found in the Ensembl output but present in the GVF file.
