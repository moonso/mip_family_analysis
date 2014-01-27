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
from multiprocessing import JoinableQueue, Queue, Lock, cpu_count
from datetime import datetime
from pprint import pprint as pp

from mip_family_analysis.family import family_parser
from mip_family_analysis.variants import variant_parser
from mip_family_analysis.models import genetic_models, score_variants
from mip_family_analysis.utils import variant_consumer, variant_sorter, header_parser

def get_family(args):
    """Return the family"""
    family_type = 'cmms'
    family_file = args.family_file[0]
    
    my_family_parser = family_parser.FamilyParser(family_file, family_type)
    # Stupid thing but for now when we only look at one family
    return my_family_parser.families.popitem()[1]

def get_header(variant_file):
    """Return a fixed header parser"""
    head = header_parser.HeaderParser(variant_file)
    head.header.append('Compounds')
    head.header.append('Rank_score')
    return head
    

def main():
    parser = argparse.ArgumentParser(description="Parse different kind of ped files.")
    parser.add_argument('family_file', type=str, nargs=1, help='A pedigree file. Default is cmms format.')
    parser.add_argument('variant_file', type=str, nargs=1, help='A variant file.Default is vcf format')
    
    parser.add_argument('-out', '--outfile', type=str, nargs=1, help='Specify the path to output, if no file specified the output will be printed to screen.')
    
    parser.add_argument('-v', '--verbose', action="store_true", help='Increase output verbosity.')
    
    parser.add_argument('-ga', '--gene_annotation', type=str, choices=['Ensembl', 'HGNC'], nargs=1, default=['HGNC'], help='What gene annotation should be used, HGNC or Ensembl.')
    
    parser.add_argument('-o', '--output', type=str, nargs=1, help='Specify the path to a file where results should be stored.')

    parser.add_argument('-pos', '--position', action="store_true", help='If output should be sorted by position. Default is sorted on rank score')

    parser.add_argument('-tres', '--treshold', type=int, nargs=1, help='Specify the lowest rank score to be outputted.')
    
    args = parser.parse_args()
    
    
    # If gene annotation is manually given:
    gene_annotation = args.gene_annotation[0]
    
    
    
    # Start by parsing at the pedigree file:
    my_family = get_family(args)
    preferred_models = my_family.models_of_inheritance
        
    # Check the variants:

    if args.verbose:
        print 'Parsing variants ...'
        print ''

    
    start_time_variant_parsing = datetime.now()
    
    var_file = args.variant_file[0]
    file_name, file_extension = os.path.splitext(var_file)
    
    # Take care of the headers from the variant file:
    head = get_header(var_file)
        
    # The task queue is where all jobs(in this case batches that represents variants in a region) is put
    # the consumers will then pick their jobs from this queue.
    tasks = JoinableQueue()
    # The consumers will put their results in the results queue
    results = Queue()
    # We will need a lock so that the consumers can print their results to screen
    lock = Lock()
    
    # Create a temporary file for the variants:
    temp_file = 'temp.tmp'
    with open(temp_file, 'w') as file_handle:
    
        num_consumers = cpu_count() * 2
        number_of_finished = 0
        
        consumers = [variant_consumer.VariantConsumer(lock, tasks, results, my_family, 
                        args.verbose) for i in xrange(num_consumers)]
    
        for w in consumers:
            w.start()
    
        var_parser = variant_parser.VariantParser(var_file, tasks, head.individuals, head.header, args.verbose)
    
    
        for i in xrange(num_consumers):
            tasks.put(None)
    
    
        # Print the results to a temporary file:
        while True:
            next_result = results.get()
            if type(next_result) == type('a'):
                if next_result == 'Done':
                    number_of_finished += 1
                if number_of_finished == num_consumers:
                    break
            else:
                for variant_id in next_result:
                    file_handle.write('\t'.join(next_result[variant_id].get_cmms_variant())+'\n')
            
        tasks.join()
        file_handle.close()

    if args.verbose:
        print 'Variants done!. Time to check models: ', (datetime.now() - start_time_variant_parsing)
        print ''
        print 'Start sorting the variants:'
        start_time_variant_sorting = datetime.now()

    if args.outfile:
        results = args.outfile[0]
    else:
        results = 'results.tmp'
    
    with open(results, 'w') as f:
        for line_number in head.metadata:
            f.write(head.metadata[line_number] + '\n')
        f.write('#' + '\t'.join(head.header) + '\n')
    
    var_sorter = variant_sorter.FileSort(temp_file, results)
    var_sorter.sort()
    os.remove(temp_file)
    
    if args.verbose:
        print 'Variants sorted!. Time to sort variants: ', (datetime.now() - start_time_variant_sorting)
        print ''
    
    if not args.outfile:
        with open(results, 'r') as f:
            for line in f:
                print line
        os.remove(results)
    
    
    
    # 
    # # Add info about variant file:
    # new_headers = my_variant_parser.header_lines 
    # 
    # # Add new headers:
    # 
    # new_headers.append('Inheritance_model')
    # new_headers.append('Compounds')
    # new_headers.append('Rank_score')
    # 
    # 
    # if args.verbose:
    #     print 'Checking genetic models...'
    #     print ''
    # 
    # for data in my_variant_parser.metadata:
    #     print data
    # 
    # print '#'+'\t'.join(new_headers)
    # 
    # if not args.position:
    #     all_variants = {}
    # 
    # # Check the genetic models
    # 
    # jobs=[]
    # 
    # for chrom in my_variant_parser.chrom_shelves:
    #     
    #     shelve_directory = os.path.split(my_variant_parser.chrom_shelves[chrom])[0]
    #     current_shelve = my_variant_parser.chrom_shelves[chrom]
    #     
    #     p = multiprocessing.Process(target=check_variants, args=(current_shelve, my_family, gene_annotation, args, preferred_models))
    #     jobs.append(p)
    #     p.start()
    # 
    # for job in jobs:
    #     job.join()
    # 
    # # Print all variants:
    # 
    # for chrom in my_variant_parser.chrom_shelves:
    #     
    #     variants = []        
    #     variant_db = shelve.open(my_variant_parser.chrom_shelves[chrom])
    #     
    #     for var_id in variant_db:
    #         variants.append(variant_db[var_id])
    #         
    #     for variant in sorted(variants, key=lambda genetic_variant:genetic_variant.start):
    #         pass
    #         # print '\t'.join(variant.get_cmms_variant())
    # 
    # 
    #     os.remove(my_variant_parser.chrom_shelves[chrom])
    # os.removedirs(shelve_directory)
    # 
    # # Else print by rank score:
    # # if not args.position:
    # #     for variant in sorted(all_variants.iteritems(), key=lambda (k,v): int(operator.itemgetter(-1)(v)), reverse=True):
    # #         if args.treshold:
    # #             rank_score = int(variant[-1][-1])
    # #             if rank_score >= args.treshold[0]:
    # #                 print '\t'.join(variant[1])
    # #         else:
    #             # print '\t'.join(variant[1])
    # if args.verbose:
    #     print 'Finished analysis!'
    #     print 'Time for analysis', (datetime.now() - start_time_variant_parsing)



if __name__ == '__main__':
    main()

