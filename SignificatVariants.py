import gzip


def write_vcf_file(vcf_file, header, body):
    # write header
    f = gzip.open(vcf_file.split('vcf')[0] + 'significant.vcf.gz', "wb")
    for line in header:
        f.write(line)
    for line in body:
        f.write(line)
    f.close()
    print('VCF file successfully has been generated!')


def find_significant_variants(vcf_file):
    """
    goal: finding significant variants in VCF file
    :param vcf_file VCF file directory
    :return: VCF file includes: significant variants
    """

    # which variants are significant?
    # 1. The variant has multiple independent dbSNP submissions,
    #  i.e. submissions with a different submitter handles or different discovery samples. (E_Multiple_observations)
    # 2. The variant is cited in a PubMed article (E_Cited)

    header = []
    body = []
    with gzip.open(vcf_file, 'rb') as f:
        for line in f:
            if not str(line)[2:-5].startswith("#"):
                if str(line).count('E_Multiple_observations') > 0 or str(line).count('E_Cited'):
                    body.append(line)
            else:
                header.append(line)
    write_vcf_file(vcf_file, header, body)


if __name__ == '__main__':
    vcf_file = 'input/canis_lupus_familiaris_CanFam3.1.vcf.gz'
    find_significant_variants(vcf_file)
