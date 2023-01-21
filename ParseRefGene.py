import gzip


def parse_ensembl_reference_genome(reference_genome, selected_chr):
    """
    goal: parsing reference genome provided by Ensembl
    @param reference_genome str reference genome directory and file name
    @param selected_chr str selected chromosome
    @return: genome  dictionary key: selected chromosome value: list of selected chromosome sequences
    (each sequence line is stored in a list)
    """
    genome = {selected_chr: []}
    read = False
    with gzip.open(reference_genome, 'rb') as f:
        for line in f:
            if str(line)[2:-3].startswith('>' + selected_chr.lstrip().rstrip()):
                read = True
            elif str(line)[2:-3].startswith('>'):
                if read:
                    break
                else:
                    pass
            else:
                if read:
                    genome[selected_chr].append(str(line)[2:-3])
    return genome


def parse_ucsc_reference_genome(reference_genome, selected_chr):

    """
    goal: parsing reference genome provided by UCSC
    @param reference_genome str reference genome directory and file name
    @param selected_chr str selected chromosome
    @return: genome  dictionary key: selected chromosome value: list of selected chromosome sequences
    (each sequence line is stored in a list)
    """
    genome = {selected_chr: []}
    read = False
    with gzip.open(reference_genome, 'rb') as f:
        for line in f:
            if str(line)[2:-3].startswith('>chr' + selected_chr.lstrip().rstrip()):
                read = True
            elif str(line)[2:-3].startswith('>chr'):
                if read:
                    break
                else:
                    pass
            else:
                if read:
                    genome[selected_chr].append(str(line)[2:-3])
    return genome


def parse_reference_genome(reference_genome, selected_chr, ref_genome_db):

    """
    goal: calling genome parser function according to reference genome db
    @param reference_genome str reference genome directory and file name
    @param selected_chr str selected chromosome
    @param ref_genome_db str 'ensembl' or 'ucsc'
    @return: genome  dictionary key: selected chromosome value: list of selected chromosome sequences
    (each sequence line is stored in a list)
    """

    if ref_genome_db == 'ensembl':
        return parse_ensembl_reference_genome(reference_genome, selected_chr)
    elif ref_genome_db == 'ucsc':
        return parse_ucsc_reference_genome(reference_genome, selected_chr)

