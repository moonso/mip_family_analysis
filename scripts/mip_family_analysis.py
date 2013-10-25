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
    
    var_file = args.variant_file[0]
    file_name, file_extension = os.path.splitext(var_file)
        
    var_type = 'cmms'        
        
    my_variant_parser = variant_parser.VariantParser(var_file, var_type)

    # Add info about variant file:
    new_headers = my_variant_parser.header_lines 
    
    # Add new headers:
    
    new_headers.append('Inheritance_model')
    new_headers.append('Rank_score')
    new_headers.append('Compounds')
    
        
    # Check the genetic models
    variants = []
    print 'Hej'
    for chrom in my_variant_parser.chrom_shelves:
        variant_db = shelve.open(my_variant_parser.chrom_shelves[chrom])
        for var_id in variant_db:
            variants.append(variant_db[var_id])
        if args.verbose:
            for variant in variants:
                print 'Before:', variant.variant_id, variant.models
        genetic_models.genetic_models(my_family, variants, gene_annotation)
        if args.verbose:
            for variant in variants:
                print 'After:',variant.variant_id, variant.models
                print ''
        print 'du'
        for variant in variants:
            if args.verbose:
                print 'Before:', variant.variant_id, variant.rank_score
            score_variants.score_variant(variant, preferred_models)
            if args.verbose:
                print 'After:',variant.variant_id, variant.rank_score
                print ''
            variant_db[variant.variant_id] = variant
        print 'glade'
        for variant in variants:
            if len(variant.ar_comp_genes) > 0:
                for gene in variant.ar_comp_genes:
                    print 'Variant_id', variant.variant_id,'Gene', gene, 'Variants', variant.ar_comp_genes[gene]
                    for compound_variant in variant.ar_comp_genes[gene]:
                        variant.ar_comp_variants.append(variant_db[compound_variant.variant_id])
        print 'kop'
        for variant in variants:
            variant_db[variant.variant_id] = variant
        print 'Dig'
        variant_db.close()
    
    

    print '\t'.join(new_headers)
    print 'en'        
    for chrom in my_variant_parser.chrom_shelves:
        variant_db = shelve.open(my_variant_parser.chrom_shelves[chrom])
        
        print 'Spade'
        if args.position:
            for variant in sorted(variant_db.keys()):
                print '\t'.join(variant_db[variant].get_cmms_variant())
        
        else:
            rank_scores = {}
            for variant_id in variant_db:
                rank_score = variant_db[variant_id].rank_score
                if rank_score in rank_scores:
                    rank_scores[rank_score].append(variant_db[variant_id])
                else:
                    rank_scores[rank_score] = [variant_db[variant_id]]
            
            for rank in sorted(rank_scores.keys()):
                for variant in rank_scores[rank]:
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

