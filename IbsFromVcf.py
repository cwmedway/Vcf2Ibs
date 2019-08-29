"""
IbdFromVcf.py
Calculates pairwise IBS from VCF file
Aurthor: Christopher Medway
Created: 29th August 2019
"""

import argparse
from pysam import VariantFile
import numpy as np
import pandas as pd
import sys


def get_args():
    """
    uses argparse package to extract command line arguments
    """

    parser = argparse.ArgumentParser(
        description='calculates pairwise identity-by-state matrix for all samples in a given VCF file')

    parser.add_argument('vcf_file', help='path to input vcf')
    parser.add_argument('output_file', help='output matrix name')

    args = parser.parse_args()
    return args


def read_vcf(vcffile):
    """
    reads vcf file line by line
    keeps a running total of IBS as a 2D array (samples x samples)
    writes matrix once all vcf lines have been processed
    """

    vcf_in = VariantFile(vcffile)
    samples = vcf_in.header.samples

    # at least two samples or exit
    if (len(samples) < 2 ):
        sys.exit("only 1 sample available in VCF")

    # initialise 2D array based on sample numbers to store ibs running total
    ibs_rntotal = [[0 for i in range(len(samples))] for j in range(len(samples))]

    # i.e. number of processed variants
    count_row = 0

    # loop over lines in vcf
    for var in vcf_in.fetch():
        # extract genotypes as list of tuples
        genotypes = [s['GT'] for s in var.samples.values()]
        # get sample x sample ibd array for given SNP
        ibs_instance = get_var_ibs(genotypes)
        # add to running total
        ibs_rntotal = np.add(ibs_instance,ibs_rntotal)
        count_row += 1

    # write out IBS matrix once vcf is done
    pd.DataFrame(
        (ibs_rntotal/count_row).round(3),
        index=samples,
        columns=samples).to_csv(args.output_file)


def get_var_ibs(genotypes):
    """
    outer loop over genotypes/samples
    """

    # stores ibs for sample x sample (2D array)
    ibs = []

    for genotype in genotypes:
        alleleA = int(genotype[0])
        alleleB = int(genotype[1])

        # returns 1D array for alleleA
        a = check_allele_match(alleleA, genotypes)
        # returns 1D array for alleleB
        b = check_allele_match(alleleB, genotypes)
        # add arrays and /4 to get IBS and push to 2D array
        ibs.append(np.add(a,b)/4)

    return(ibs)


def check_allele_match(allele,genotypes):
    """
    inner loop over genotypes for each sample and checks for number of matches
    with given allele
    """

    shared = []
    for genotype in genotypes:
        if allele == int(genotype[0]) and allele == int(genotype[1]):
            shared.append(2)
        elif allele == int(genotype[0]) or allele == int(genotype[1]):
            shared.append(1)
        else:
            shared.append(0)
    return(shared)


if __name__ == '__main__':
    args = get_args()
    read_vcf(args.vcf_file)