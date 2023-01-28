import gzip
import pandas as pd
from ParseRefGene import parse_reference_genome
from datetime import date


def seq(x):
    chromosome = x['#CHROM']
    pos = x['POS'] - 1
    quotient = (pos - 1) // len(ref_gene[chromosome][0])
    remainder = (pos - 1) % len(ref_gene[chromosome][0])
    x['POS'] = pos

    ref_nuc = ref_gene[chromosome][quotient][remainder]

    if x['TSA'] == 'deletion':
        x['REF'] = ref_nuc.upper() + x['REF']
        x['ALT'] = ref_nuc.upper()
        return x

    elif x['TSA'] == 'insertion':
        x['REF'] = ref_nuc.upper()
        x['ALT'] = ','.join([ref_nuc.upper() + k for k in x['ALT'].split(',')])
        return x

    elif x['TSA'] == 'sequence_alteration':
        x['REF'] = ref_nuc.upper() + x['REF']
        x['ALT'] = ','.join([ref_nuc.upper() + k.replace('-', '') for k in x['ALT'].split(',')])
        return x

    elif x['TSA'] == 'tandem_repeat':

        x['POS'] += 1
        alt = x['ALT'].split(',')
        alt_length = [len(k) for k in alt]
        min_length = min(min(alt_length), len(x['REF']))
        x['REF'] = x['REF'][min_length - 1:]
        x['ALT'] = ','.join([k[min_length - 1:] for k in alt])
        return x

    else:
        return x


def write_vcf_file(gvf_file, vcf_df, dbSNP_v, ref_genome_name, source):
    today = date.today()
    today.strftime("%Y%m%d")
    header = '##fileformat=VCFv4.1\n' \
             '##fileDate=' + str(today).replace("-", "") + '\n' \
                                                           '##' + source + '\n' \
                                                                           '##reference=' + ref_genome_name + '\n' \
                                                                                                              '##INFO=<ID=' + dbSNP_v + ',Number=0,Type=Flag,Description="Variants (including SNPs and indels) ' \
                                                                                                                                        'imported from dbSNP [Remapped to ' + ref_genome_name + ']">\n' \
                                                                                                                                                                                                '##INFO=<ID=TSA,Number=1,Type=String,Description="Type of sequence alteration. Child of term sequence_alteration as defined by the sequence ontology project.">\n' \
                                                                                                                                                                                                '##INFO=<ID=E_Cited,Number=0,Type=Flag,Description="Cited.https://www.ensembl.org/info/genome/variation/prediction/variant_quality.html#evidence_status">\n' \
                                                                                                                                                                                                '##INFO=<ID=E_Multiple_observations,Number=0,Type=Flag,Description="Multiple_observations.https://www.ensembl.org/info/genome/variation/prediction/variant_quality.html#evidence_status">\n' \
                                                                                                                                                                                                '##INFO=<ID=E_Freq,Number=0,Type=Flag,Description="Frequency.https://www.ensembl.org/info/genome/variation/prediction/variant_quality.html#evidence_status">\n' \
                                                                                                                                                                                                '##INFO=<ID=E_TOPMed,Number=0,Type=Flag,Description="TOPMed.https://www.ensembl.org/info/genome/variation/prediction/variant_quality.html#evidence_status">\n' \
                                                                                                                                                                                                '##INFO=<ID=E_Hapmap,Number=0,Type=Flag,Description="HapMap.https://www.ensembl.org/info/genome/variation/prediction/variant_quality.html#evidence_status">\n' \
                                                                                                                                                                                                '##INFO=<ID=E_Phenotype_or_Disease,Number=0,Type=Flag,Description="Phenotype_or_Disease.https://www.ensembl.org/info/genome/variation/prediction/variant_quality.html#evidence_status">\n' \
                                                                                                                                                                                                '##INFO=<ID=E_ESP,Number=0,Type=Flag,Description="ESP.https://www.ensembl.org/info/genome/variation/prediction/variant_quality.html#evidence_status">\n' \
                                                                                                                                                                                                '##INFO=<ID=E_gnomAD,Number=0,Type=Flag,Description="gnomAD.https://www.ensembl.org/info/genome/variation/prediction/variant_quality.html#evidence_status">\n' \
                                                                                                                                                                                                '##INFO=<ID=E_1000G,Number=0,Type=Flag,Description="1000Genomes.https://www.ensembl.org/info/genome/variation/prediction/variant_quality.html#evidence_status">\n' \
                                                                                                                                                                                                '##INFO=<ID=E_ExAC,Number=0,Type=Flag,Description="ExAC.https://www.ensembl.org/info/genome/variation/prediction/variant_quality.html#evidence_status">\n'

    # write header
    f = open(gvf_file.split('gvf')[0] + 'vcf', "w")
    f.write(header)
    f.close()
    # write VCF file
    # sort df
    # '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'
    vcf_df = vcf_df.sort_values(by=['#CHROM', 'POS'], ascending=True)
    vcf_df.to_csv(gvf_file.split('gvf')[0] + 'vcf', index=False, mode='a', sep='\t')
    print('VCF file successfully has been generated!')


def convert_vcf(gvf_file, gvf_df, reference_genome_file, dbSNP_v, ref_genome_name, ref_genome_db, source):
    # add new columns
    gvf_df['QUAL'] = '.'
    gvf_df['FILTER'] = '.'
    gvf_df['ref_N'] = ''
    gvf_df['INFO'] = 'TSA=' + gvf_df['TSA'] + ';' + gvf_df['DB'] + ';' + gvf_df['INFO']
    INDEL_df = gvf_df.loc[gvf_df['TSA'] != 'SNV']
    global ref_gene
    if INDEL_df.shape[0] > 0:
        # read reference genome
        for chr in INDEL_df['#CHROM'].unique():
            print('chromosome: ' + chr)
            ref_gene = parse_reference_genome(reference_genome_file, chr, ref_genome_db)
            chr_df = gvf_df.loc[(gvf_df['TSA'] != 'SNV') & (gvf_df['#CHROM'] == chr)]

            if chr_df.shape[0] > 0:
                # modify REF & ALT
                gvf_df.loc[(gvf_df['TSA'] != 'SNV') & (gvf_df['#CHROM'] == chr)] = \
                    gvf_df.loc[(gvf_df['TSA'] != 'SNV') & (gvf_df['#CHROM'] == chr)].apply(seq, axis=1)

    # select columns
    if genome_browser == 'ucsc':
        gvf_df['#CHROM'] = 'chr' + gvf_df['#CHROM']

    vcf_df = gvf_df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']]

    # write to csv file
    write_vcf_file(gvf_file, vcf_df, dbSNP_v, ref_genome_name, source)


def parse_gvf_file(gvf_file, reference_genome_file, ref_genome_db, chr_name_list):
    """
    goal: parsing GVF file format
    :param gvf_file GVF file directory
    @param reference_genome_file reference genome file directory
    @param ref_genome_db str 'ensembl' or 'ucsc'
    :return: VCF file
    """
    body = []
    ref_genome_name = ''
    source = ''
    global genome_browser
    genome_browser = ref_genome_db
    with gzip.open(gvf_file, 'rb') as f:
        for line in f:
            if not str(line)[2:-3].startswith("#"):
                splitted_line = str(line)[2:-3].split('\\t')
                info_list = {k[0]: k[1] for k in list(map(lambda x: x.split('='), splitted_line[-1].split(';')))}
                [db, db_id] = info_list['Dbxref'].split(':')
                info = ';'.join(['E_' + info_list[k] for k in info_list.keys()
                                 if k not in ['ID', 'Variant_seq', 'Dbxref', 'Reference_seq']])

                if splitted_line[0] in chr_name_list and splitted_line[2] in ['deletion', 'insertion',
                                                                              'sequence_alteration', 'SNV']:
                    body.append(splitted_line[:-1] + [info_list['ID'], info_list['Variant_seq'], db, db_id,
                                                      info_list['Reference_seq'], info])
            else:
                if str(line)[2:-3].startswith("##genome-build"):
                    ref_genome_name = str(line)[2:-3].split("##genome-build")[1].rstrip().lstrip()
                elif str(line)[2:-3].startswith("##data-source"):
                    source = 's' + str(line)[2:-3].split("data-source S")[1].rstrip().lstrip()
    dbSNP_v = db
    gvf_df = pd.DataFrame(body, columns=['#CHROM', 'db', 'TSA', 'POS', 'end_pos', 'unknown_flag', 'SNV_flag',
                                         'INDEL_flag', 'id', 'ALT', 'DB', 'ID', 'REF', 'INFO'])
    gvf_df['POS'] = gvf_df['POS'].astype(int)
    # convert column type
    convert_vcf(gvf_file, gvf_df, reference_genome_file, dbSNP_v, ref_genome_name, ref_genome_db, source)
