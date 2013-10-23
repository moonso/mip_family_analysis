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
from Mip_Family_Analysis.Variants import genetic_variant, genotype
from collections import OrderedDict

class VariantParser(object):
    """Parses a file with family info and creates a family object with individuals."""
    def __init__(self, infile, file_type):
        super(VariantParser, self).__init__()
        self.variant_type = file_type
        self.variants = OrderedDict() # DICT with {variant_id:<variant_object>}
        self.header_lines = []
        self.individuals = {} # DICT with {ind_id_1:{variant_id_1:genotype}}
        self.metadata = []
        
        with open(infile, 'r') as f:
            line_count = 0
            for line in f:
                line = line.rstrip()
                line_count += 1
                if line[0] != '#' and len(line) > 1:
                    if self.variant_type == 'vcf':
                        my_variant = self.vcf_variant(line)
                    elif self.variant_type == 'cmms':
                        my_variant = self.cmms_variant(line)
                    # elif self.variant_type == 'annovar':
                    #     my_variant = get_variant.annovar_variant(variant_line)
                elif line[:2] != '##':
                    # If necesary we can write something to parse the headers.
                    self.header_lines = line[1:].split()
                    self.metadata.append(line)
                elif line[:2] == '##':
                    self.metadata.append(line)
                    
    
    def cmms_variant(self, variant_line):
        """Returns a Variant objekt from the cmms variant format. Creates a list with Genotype objects that becomes a member of the variant."""        
        variant_line.rstrip()
        variant_info = OrderedDict()
        individual_genotypes = {} #DICT with {ind_id:{<genotype>}}
        counter = 0
        variant_line = variant_line.split('\t')
        
        for entry in range(len(variant_line)):
            
            if 'IDN' in self.header_lines[entry]:
                individual = self.header_lines[entry].split(':')[-1]
                if individual not in self.individuals:
                    self.individuals[individual] = {}
                variant_info[individual] = variant_line[entry]
            else:
                variant_info[self.header_lines[entry]] = variant_line[entry]
                
        counter += 1
        chrom = variant_info['Chromosome']
        start = variant_info['Variant_start']
        stop = variant_info['Variant_stop']
        alternative = variant_info['Alternative_allele']
        reference = variant_info['Reference_allele']
        identity = variant_info['Dbsnp_rs_nr']
        my_variant = genetic_variant.Variant(chrom, start, stop, reference, alternative, identity, variant_info)
        self.variants[my_variant.variant_id] = my_variant
        for individual in individual_genotypes:
            self.individuals[individual][my_variant.variant_id] = individual_genotypes[individual]
        
        # Check the genotypes:
        
        for individual in self.individuals:
            genotype_arguments = {} # args for genotype class
            for gt_info in variant_info[individual].split(':')[1:]:
                value_pair = gt_info.split('=')
                genotype_arguments[value_pair[0]] = value_pair[-1]
            my_genotype = genotype.Genotype(GT=genotype_arguments.get('GT','./.'), AD=genotype_arguments.get('AD','.,.'), DP=genotype_arguments.get('DP','0'), GQ=genotype_arguments.get('GQ','0'))
            self.individuals[individual][my_variant.variant_id] = my_genotype
        


def main():
    parser = argparse.ArgumentParser(description="Parse different kind of pedigree files.")
    parser.add_argument('variant_file', type=str, nargs=1 , help='A file with variant information.')
    parser.add_argument('-vcf', '--vcf', action="store_true", help='Variant file is in vcf format.')
    parser.add_argument('-cmms', '--cmms', action="store_true", help='Cmms annotated format.')
    args = parser.parse_args()
    infile = args.variant_file[0]
    file_type = 'vcf'
    if args.vcf:
        file_type = 'vcf'
    if args.cmms:
        file_type = 'cmms'
    my_parser = VariantParser(infile, file_type)
    for individual in my_parser.individuals:
        for genotype in my_parser.individuals[individual]:
            print individual, genotype, my_parser.individuals[individual][genotype]
    


if __name__ == '__main__':
    main()
