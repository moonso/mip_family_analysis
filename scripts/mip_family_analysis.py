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
import operator
import multiprocessing
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

    parser.add_argument('-tres', '--treshold', type=int, nargs=1, help='Specify the lowest rank score to be outputted.')
    
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
    new_headers.append('Compounds')
    new_headers.append('Rank_score')
    
    start_time_genetic_models = datetime.now()
    
    if args.verbose:
        print 'Checking genetic models...'
        print ''
    
    for data in my_variant_parser.metadata:
        print data
    
    print '#'+'\t'.join(new_headers)
    
    if not args.position:
        all_variants = {}
    # Check the genetic models
    for chrom in my_variant_parser.chrom_shelves:
        
        variants = []
        variant_db = shelve.open(my_variant_parser.chrom_shelves[chrom])
        variant_dict = {}
        for var_id in variant_db:
            variants.append(variant_db[var_id])
        
        variants = genetic_models.check_genetic_models(my_family, variants, gene_annotation, verbose = args.verbose)
        
        if args.verbose:
            for variant in variants:
                print variant.models
            print 'Models checked!. Time to check models: ', (datetime.now() - start_time_genetic_models)
            print ''
            
        # Score the variants
        for variant in variants:
            score_variants.score_variant(variant, preferred_models)
            variant_dict[variant.variant_id] = variant
    
        # Score the compound pairs:
        for variant in variants:
            if len(variant.ar_comp_variants) > 0:
                for compound_variant_id in variant.ar_comp_variants:
                    comp_score = variant.rank_score + variant_dict[compound_variant_id].rank_score
                    variant.ar_comp_variants[compound_variant_id] = comp_score
            if not args.position:
                all_variants[variant.variant_id] = variant.get_cmms_variant()
        
        # Print by position if desired
        if args.position:
            for variant in sorted(variants, key=lambda genetic_variant:genetic_variant.start):
                print '\t'.join(variant.get_cmms_variant())

        variant_db.close()
        os.remove(my_variant_parser.chrom_shelves[chrom])

    # Else print by rank score:
    if not args.position:
        for variant in sorted(all_variants.iteritems(), key=lambda (k,v): int(operator.itemgetter(-1)(v)), reverse=True):
            if args.treshold:
                rank_score = int(variant[-1][-1])
                if rank_score >= args.treshold[0]:
                    print '\t'.join(variant[1])
            else:
                print '\t'.join(variant[1])



if __name__ == '__main__':
    main()

