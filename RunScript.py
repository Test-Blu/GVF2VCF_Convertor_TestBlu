import time
from GVF2VCF import parse_gvf_file


if __name__ == '__main__':

    # input files
    gvf_file = '../CanFam3/canis_lupus_familiaris_CanFam3.1.gvf.gz'
    reference_genome_file = '../CanFam3/canFam3.fa.gz'
    chr_name_list = [str(k) for k in list(range(1, 39))] + ['X']
    start_time = time.time()
    # run script
    # ref_genome_db can be: 'ensembl' or 'ucsc'
    ref_genome_db = 'ucsc'
    parse_gvf_file(gvf_file, reference_genome_file, ref_genome_db, chr_name_list)
    print("--- %s seconds ---" % (time.time() - start_time))


