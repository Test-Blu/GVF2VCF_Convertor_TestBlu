# GVF2VCF_Convertor_TestBlu
 Super fast conversion of GVF file format to VCF file format
 
Finding robust variants (predominantly somatic variants) takes some additional resources.

To detect somatic variants, it is preferable to have germline variants known in the studied organism so that germline variants can be removed from the VCF output file of somatic variants.
</br>
These resources can be accessed at NCBI or other relevant databases for human genomes, but for non-human organisms, they are only available in a limited form in some databases.
</br>
Germline variants are only reported in GVF (genome variant format) files for some non-human organisms in the Ensembl database

As the structure of this file differs from the VCF file, it should be converted to the VCF format for use with GATK.

For this reason, Ensembl developers provided an API to convert GVF files to VCF in 2014, but this API has some severe problems:

1. Convergence is very slow (it took about 6 hours in one tested sample!)
2. Because many people connect to the database through different Ensembl APIs, the conversion may be done incompletely.
3. By comparing the GVF, and VCF files of Ensembl's API output, I realized that some germline variants were not reported in the VCF output!
4. It only supports the header pattern of the Ensembl reference genome.

For these reasons, we wrote a script in Python that performs an offline conversion.
Its features include the following:
1. It performs the conversion in 20 minutes (about six gigs of RAM is required)
2. In addition to converting GVF to VCF, significant variants, i.e. reported variants in articles or seen in many samples, can be reported in a separate file.
3. It supports the header pattern of the UCSC and Ensembl reference genome.

In addition, by comparing the output VCF file of this script and the output VCF file of Ensembl API, we found that the output VCF file of our script reports all the variants in the output VCF file of Ensembl, as well as variants not found in the Ensembl output VCF file (but found in the GVF file).
