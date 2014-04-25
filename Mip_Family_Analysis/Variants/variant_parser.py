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
import re

from pprint import pprint as pp
from datetime import datetime

from Mip_Family_Analysis.Variants import genotype

class VariantFileParser(object):
    """docstring for VariantParser"""
    def __init__(self, variant_file, batch_queue, head, verbosity = False):
        super(VariantFileParser, self).__init__()
        self.variant_file = variant_file
        self.batch_queue = batch_queue
        self.verbosity = verbosity
        self.individuals = head.individuals
        self.header_line = head.header
    
    def parse(self):
        """Start the parsing"""        
        start_parsing = datetime.now()
        start_chrom = start_parsing
        start_twenty = start_parsing
        beginning = True
        batch = {}
        new_chrom = None
        current_chrom = None
        current_features = []
        nr_of_variants = 0
        if self.verbosity:
            print('Start parsing the variants ...\n')
        with open(self.variant_file, 'rb') as f:
            for line in f:
                if not line.startswith('#'):
                    variant_line = line.rstrip().split('\t')
                    variant, new_features = self.cmms_variant(variant_line, self.individuals)
                    if self.verbosity:
                        nr_of_variants += 1
                        new_chrom = variant['Chromosome']
                        if nr_of_variants % 20000 == 0:
                            print('%s variants parsed!' % str(nr_of_variants))
                            print('Last 20.000 took %s to parse. \n' % str(datetime.now() - start_twenty))
                            start_twenty = datetime.now()
                    # If we look at the first variant, setup boundary conditions:
                    if beginning:
                        current_features = new_features
                        beginning = False
                        # Add the variant to each of its features in a batch
                        batch = self.add_variant(batch, variant, new_features)
                        if self.verbosity:
                            current_chrom = new_chrom
                    else:
                        send = True
                    
                    # Check if we are in a space between features:
                        if len(new_features) == 0:
                            if len(current_features) == 0:
                                send = False
                    #If not check if we are in a consecutive region
                        elif len(set.intersection(set(new_features),set(current_features))) > 0:
                            send = False
                        
                        if send:
                            # If there is an intergenetic region we do not look at the compounds.
                            # The tasks are tuples like (variant_list, bool(if compounds))
                            self.batch_queue.put(batch)
                            current_features = new_features
                            batch = self.add_variant({}, variant, new_features)
                        else:
                            current_features = list(set(current_features) | set(new_features))
                            batch = self.add_variant(batch, variant, new_features) # Add variant batch
                    
                    if self.verbosity:
                        if new_chrom != current_chrom:
                            print('Chromosome %s parsed!' % current_chrom)
                            print('Time to parse chromosome %s' % str(datetime.now()-start_chrom))
                            current_chrom = new_chrom
                            start_chrom = datetime.now()
                        
        if self.verbosity:
            print('Chromosome %s parsed!' % current_chrom)
            print('Time to parse chromosome %s \n' % str(datetime.now()-start_chrom))
            print('Variants done!. Time to parse variants: %s \n' % str(datetime.now() - start_parsing))
        self.batch_queue.put(batch)
        
        return
    
    def add_variant(self, batch, variant, features):
        """Adds the variant to the proper gene(s) in the batch."""
        variant_id = [variant['Chromosome'], variant['Variant_start'], variant['Reference_allele'], variant['Alternative_allele']]
        variant_id = '_'.join(variant_id)
        if len(features) == 0:
            if len(batch) == 0:
                batch['-'] = {variant_id:variant}
            else:
                batch['-'][variant_id] = variant
        for feature in features:
            if feature in batch:
                batch[feature][variant_id] = variant
            else:
                batch[feature] = {variant_id:variant}
        return batch
    
    def get_genes(self, annotation_string, annotation_type = 'HGNC'):
        """Parse the annotation and return a list of genes that the variant belongs to. annotation_type in [HGNC, ENSEMBL]"""
        genes = []
        if annotation_type == 'ENSEMBL':
            return annotation_string.split(';')
        elif annotation_type == 'HGNC':
            if 'dist' in annotation_string:
                return genes
            hgnc_genes = {}
            for gene_string in annotation_string.split(';'):
                new_gene_string = re.sub('\(.*\)','', gene_string)
                for gene in new_gene_string.split(','):
                    hgnc_genes[gene] = ''
                return list(hgnc_genes.keys())
        return genes
    
    def cmms_variant(self, splitted_variant_line, individuals):
        """Returns a variant object in the cmms format."""
        
        variant = dict(zip(self.header_line, splitted_variant_line))
        
        # Get the genes:
        features_overlapped = self.get_genes(variant['HGNC_symbol'], 'HGNC')
        
        variant['Genotypes'] = {}
        
        for individual in individuals:
            try:
                gt_info = variant['IDN:'+individual].split(':')[1].split('=')[1]
            except (IndexError, KeyError):
                gt_info = './.'
            
            variant['Genotypes'][individual] = genotype.Genotype(GT=gt_info)
        
        return variant, features_overlapped
    


def main():
    from multiprocessing import JoinableQueue
    from Mip_Family_Analysis.Utils import header_parser
    from Mip_Family_Analysis.Utils import annotation_parser
    parser = argparse.ArgumentParser(description="Parse different kind of pedigree files.")
    parser.add_argument('variant_file', type=str, nargs=1 , help='A file with variant information.')
    parser.add_argument('-ann', '--annotation', nargs=1 , help='A file with annotations.')    
    parser.add_argument('-v', '--verbose', action="store_true", help='Increase output verbosity.')
    args = parser.parse_args()
    
    infile = args.variant_file[0]
    
    annotation_trees = {}
    if args.annotation:
        print('Parsing annotation:')
        annotation_trees = annotation_parser.AnnotationParser(args.annotation[0], 'ref_gene')
        print('Annotation parsed.')
        
    head = header_parser.HeaderParser(infile)
    file_type = 'cmms'
    variant_queue = JoinableQueue()
    start_time = datetime.now()
    my_parser = VariantFileParser(infile, variant_queue, head, annotation_trees, args.verbose)
    my_parser.parse()
    # print(('Time to parse variants: %s' % (datetime.now()-start_time)))

if __name__ == '__main__':
    main()
