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
from datetime import datetime

from Mip_Family_Analysis.Family import family_parser
from Mip_Family_Analysis.Variants import variant_parser
from Mip_Family_Analysis.Models import genetic_models, score_variants


def main():
    parser = argparse.ArgumentParser(description="Parse different kind of ped files.")
    parser.add_argument('family_file', type=str, nargs=1, help='A pedigree file. Default is ped format.')
    parser.add_argument('variant_file', type=str, nargs=1, help='A variant file.Default is vcf format')
    
    parser.add_argument('-v', '--verbose', action="store_true", help='Increase output verbosity.')
    
    parser.add_argument('-ga', '--gene_annotation', type=str, choices=['Ensembl', 'HGNC'], nargs=1, help='What gene annotation should be used, HGNC or Ensembl.')
    
    parser.add_argument('-o', '--output', type=str, nargs=1, help='Specify the path to a file where results should be stored.')

    parser.add_argument('-pos', '--position', action="store_true", help='If output should be sorted by position. Default is sorted on rank score')
    
    args = parser.parse_args()
    
    gene_annotation = 'Ensembl'
    
    if args.gene_annotation:
       gene_annotation = args.gene_annotation[0]
    
    new_headers = []    
        
    # Start by parsing at the pedigree file:
    family_type = 'ped'
    family_file = args.family_file[0]
        
    my_family_parser = family_parser.FamilyParser(family_file, family_type)
    
    preferred_models = my_family_parser.preferred_models
    
    # Stupid thing but for now when we only look at one family
    my_family = my_family_parser.families.popitem()[1]
        
    # Check the variants:
    
    start_time_variant_parsing = datetime.now()
    
    var_file = args.variant_file[0]
    file_name, file_extension = os.path.splitext(var_file)
        
    var_type = 'cmms'        
        
    my_variant_parser = variant_parser.VariantParser(var_file, var_type)


    if args.verbose:
        print 'Variants done!. Time to parse variants: ', (datetime.now() - start_time_variant_parsing)
        print ''

    # Add info about variant file:
    new_headers = my_variant_parser.header_lines 
    
    # Add new headers:
    
    new_headers.append('Inheritance_model')
    new_headers.append('Rank_score')
    new_headers.append('Compounds')
    
    start_time_genetic_models = datetime.now()

        
    # Check the genetic models
    variants = []
    for chrom in my_variant_parser.chrom_shelves:
        variant_db = shelve.open(my_variant_parser.chrom_shelves[chrom])
        print 'Length of ',chrom , len(variant_db)
        for var_id in variant_db:
            variants.append(variant_db[var_id])
        genetic_models.genetic_models(my_family, variants, gene_annotation)
        
        if args.verbose:
            print 'Models checked!. Time to check models: ', (datetime.now() - start_time_genetic_models)
            print ''

        for variant in variants:
            score_variants.score_variant(variant, preferred_models)
            variant_db[variant.variant_id] = variant
        for variant in variants:
            if len(variant.ar_comp_genes) > 0:
                for gene in variant.ar_comp_genes:
                    print 'Comp!', gene, variant.ar_comp_genes[gene]
                    for compound_variant in variant.ar_comp_genes[gene]:
                        variant.ar_comp_variants.append(variant_db[compound_variant.variant_id])
                variant.ar_comp_genes = {}
        for variant in variants:
            variant_db[variant.variant_id] = variant
        
        start_time_close_db = datetime.now()

        variant_db.close()
    
        if args.verbose:
            print 'db '+ my_variant_parser.chrom_shelves[chrom] +' closed!. Time to close db: ', (datetime.now() - start_time_close_db)
            print ''
    

    print '\t'.join(new_headers)
    variants = []
    for chrom in my_variant_parser.chrom_shelves:
        variant_db = shelve.open(my_variant_parser.chrom_shelves[chrom])
        for var_id in variant_db:
            variants.append(variant_db[var_id])
        
        if args.position:
            for variant in sorted(variants, key=lambda genetic_variant:genetic_variant.variant_id):
                print '\t'.join(variant.get_cmms_variant())
        
        else:
            for variant in sorted(variants, key=lambda genetic_variant:genetic_variant.rank_score, reverse = True):
                print '\t'.join(variant.get_cmms_variant())
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

