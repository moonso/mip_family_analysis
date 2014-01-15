#!/usr/bin/env python
# encoding: utf-8
"""
variant_parser.py


Parse a file with variant info, this can be a .vcf file, an annotated annovar file, 
a annotated .txt cmms file, a annotated .txt cmms_ranked .

Create a variant objects and a dictionary with individuals that have a dictionary with genotypes for each variant.

Created by Måns Magnusson on 2013-01-17.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import argparse
import shelve
import tempfile
from Mip_Family_Analysis.Variants import genetic_variant, genotype
from collections import OrderedDict

def variant_parser(variant_batch, header_line, individuals):
    """Take a batch of variants and create variant ocbjects and returns them in a list"""
    # variant_type = file_type
    # individuals = []
    #Dictionary with <Gene ID> : [variant_id_1, variant_id_2, ...] for controlling the compound heterozygotes
    
    variants = []
    
    for variant_line in variant_batch:
            
        variant_info = OrderedDict()
        individual_genotypes = {} #DICT with {ind_id:{<genotype>}}
        counter = 0
        variant_line = variant_line.split('\t')
    
        for entry in range(len(variant_line)):
            
            if 'IDN' in header_line[entry]:
                # Looks like IDN:11-1-2A
                individual = header_line[entry].split(':')[-1]
                if individual not in individuals:
                    individuals.append(individual)
            variant_info[header_line[entry]] = variant_line[entry]
                
        chrom = variant_info['Chromosome']
        start = variant_info['Variant_start']
        stop = variant_info['Variant_stop']
        alternative = variant_info['Alternative_allele']
        reference = variant_info['Reference_allele']
        identity = variant_info['Dbsnp_rs_nr']
        my_variant = genetic_variant.Variant(chrom, start, stop, reference, alternative, identity, variant_info)
    
        # Add the genotypes to variant:
        
        for individual in individuals:
            genotype_arguments = {} # args for genotype class
            key = 'IDN:' + individual
            # gt_info looks like 11-1-2A:GT=0/1:PL=32,3,2:...
            for gt_info in variant_info[key].split(':')[1:]:
                value_pair = gt_info.split('=')
                genotype_arguments[value_pair[0]] = value_pair[-1]
            my_genotype = genotype.Genotype(GT=genotype_arguments.get('GT','./.'), AD=genotype_arguments.get('AD','.,.'), DP=genotype_arguments.get('DP','0'), GQ=genotype_arguments.get('GQ','0'))
            my_variant.genotypes[individual] = my_genotype
        variants.append(my_variant)
   
    return variants


def main():
    parser = argparse.ArgumentParser(description="Parse different kind of pedigree files.")
    parser.add_argument('variant_file', type=str, nargs=1 , help='A file with variant information.')
    args = parser.parse_args()
    infile = args.variant_file[0]
    file_type = 'cmms'
    my_parser = VariantParser(infile, file_type)
    number_of_variants = 0
    for chrom in my_parser.chrom_shelves:
        current_db = shelve.open(my_parser.chrom_shelves[chrom])
        print chrom, len(current_db)
        for variant in current_db:
            print current_db[variant]
        number_of_variants += len(current_db)
        current_db.close()
        os.remove(my_parser.chrom_shelves[chrom])
    print 'Number of shelved:', number_of_variants


if __name__ == '__main__':
    main()
