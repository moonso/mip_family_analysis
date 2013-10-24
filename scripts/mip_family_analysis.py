#!/usr/bin/env python
# encoding: utf-8
"""
mip_family_analysis.py

Script for annotating which genetic models that are followed for variants in the mip pipeline.

Created by Måns Magnusson on 2013-01-31.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import argparse
import shelve

from Mip_Family_Analysis.Family import family_parser
from Mip_Family_Analysis.Variants import variant_parser
from Mip_Family_Analysis.Models import genetic_models, score_variants


def main():
    parser = argparse.ArgumentParser(description="Parse different kind of ped files.")
    parser.add_argument('family_file', type=str, nargs=1, help='A pedigree file. Default is ped format.')
    parser.add_argument('variant_file', type=str, nargs=1, help='A variant file.Default is vcf format')
    
    parser.add_argument('-v', '--verbose', action="store_true", help='Increase output verbosity.')
    
    parser.add_argument('-o', '--output', type=str, nargs=1, help='Specify the path to a file where results should be stored.')

    args = parser.parse_args()
    
    new_headers = []    
        
    # Start by parsing at the pedigree file:
    family_type = 'ped'
    family_file = args.family_file[0]
        
    my_family_parser = family_parser.FamilyParser(family_file, family_type)
    
    preferred_models = my_family_parser.preferred_models
    
    # Stupid thing but for now when we only look at one family
    my_family = my_family_parser.families.popitem()[1]
        
    # Check the variants:
    
    var_file = args.variant_file[0]
    file_name, file_extension = os.path.splitext(var_file)
        
    var_type = 'cmms'        
        
    my_variant_parser = variant_parser.VariantParser(var_file, var_type)

    # Add info about variant file:
    new_headers = my_variant_parser.header_lines 
    
    # Add new headers:
    
    new_headers.append('Inheritance_model')
    new_headers.append('Rank_score')
    
    print '\t'.join(new_headers)
        
    # Check the genetic models
    for chrom in my_variant_parser.chrom_shelves:
        variants = shelve.open(my_variant_parser.chrom_shelves[chrom])
        if args.verbose:
            for variant in variants:
                print variants[variant]
        genetic_models.genetic_models(my_family, variants)
        for variant_id in variants:
            variant = variants[variant_id]
            score_variants.score_variant(variant, preferred_models)
            variants[variant_id] = variant
        for variant in sorted(variants.keys()):
            print '\t'.join(variants[variant].get_cmms_variant())
        variants.close()
        os.remove(my_variant_parser.chrom_shelves[chrom])

        
    
    #     for individual in my_family.individuals:
    #         print individual, my_family.individuals[individual].genotypes[variant]
    # 
    # print 'Score variants:'
    # my_family.score_variants()
    # print 'Variants scored!'
    # # for variant in my_family.variants:
    # #     my_family.variants[variant].print_model_info()


if __name__ == '__main__':
    main()

